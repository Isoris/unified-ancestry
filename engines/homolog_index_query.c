// =============================================================================
// homolog_index_query.c — microsecond gene→hits lookup over a homolog index
//                          built by homolog_index_build. mmap + binary search.
//
// Usage:
//   homolog_index_query --index <f.holindx> --gene <id>
//                       [--format tsv|json] [--source diamond|miniprot]
//                       [--min_identity F] [--min_bitscore F] [--limit N]
//                       [--no_header]
//
// Output JSON shape (also doc'd in
// engines/schemas/homolog_index.query.output.schema.json):
//   {"query":"GENE_ID","found":true,"n_hits":N,"hits":[{...},{...}]}
//
// Compile: gcc -O3 -march=native -o homolog_index_query homolog_index_query.c
// =============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "homolog_index.h"

typedef struct {
    int fd;
    void* base;
    size_t size;
    const IndexHeader* hdr;
    const GeneRecord* genes;
    const HitRecord* hits;
    const char* pool;
} Index;

static int idx_open(const char* path, Index* idx) {
    idx->fd = open(path, O_RDONLY);
    if (idx->fd < 0) return -1;
    struct stat st;
    if (fstat(idx->fd, &st) < 0) { close(idx->fd); return -2; }
    idx->size = (size_t)st.st_size;
    idx->base = mmap(NULL, idx->size, PROT_READ, MAP_PRIVATE, idx->fd, 0);
    if (idx->base == MAP_FAILED) { close(idx->fd); return -3; }
    idx->hdr = (const IndexHeader*)idx->base;
    if (idx->hdr->magic != HOLINDX_MAGIC)   { munmap(idx->base, idx->size); close(idx->fd); return -4; }
    if (idx->hdr->version != HOLINDX_VERSION) { munmap(idx->base, idx->size); close(idx->fd); return -5; }
    idx->genes = (const GeneRecord*)((const char*)idx->base + idx->hdr->off_gene_index);
    idx->hits  = (const HitRecord*) ((const char*)idx->base + idx->hdr->off_hits);
    idx->pool  = (const char*)      idx->base + idx->hdr->off_string_pool;
    return 0;
}

static void idx_close(Index* idx) {
    if (idx->base) munmap(idx->base, idx->size);
    if (idx->fd >= 0) close(idx->fd);
}

static const GeneRecord* idx_find(const Index* idx, const char* name) {
    int lo = 0, hi = (int)idx->hdr->n_genes - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        int cmp = strcmp(idx->pool + idx->genes[mid].name_off, name);
        if (cmp == 0) return &idx->genes[mid];
        if (cmp < 0) lo = mid + 1;
        else         hi = mid - 1;
    }
    return NULL;
}

static const char* src_type_str(uint8_t s) {
    return s == HOLINDX_SRC_DIAMOND ? "diamond" : (s == HOLINDX_SRC_MINIPROT ? "miniprot" : "unknown");
}

static void emit_tsv_header(FILE* out) {
    fprintf(out,
        "query\ttarget\tsource_label\tsource_type\tidentity\tbitscore\tevalue\t"
        "qstart\tqend\ttstart\ttend\talnlen\tn_segments\tstrand\n");
}

static void emit_tsv_hit(FILE* out, const Index* idx, const GeneRecord* g, const HitRecord* h) {
    fprintf(out, "%s\t%s\t%s\t%s\t%.4f\t%.2f\t%.3g\t%d\t%d\t%d\t%d\t%u\t%u\t%c\n",
            idx->pool + g->name_off,
            idx->pool + h->target_off,
            idx->pool + h->source_off,
            src_type_str(h->source_type),
            (double)h->identity, (double)h->bitscore, h->evalue,
            h->qstart, h->qend, h->tstart, h->tend,
            h->alnlen, (unsigned)h->n_segments,
            h->strand ? h->strand : '.');
}

// Minimal JSON string escape: escape backslash, double-quote, and control chars.
static void emit_json_str(FILE* out, const char* s) {
    fputc('"', out);
    for (; *s; s++) {
        unsigned char c = (unsigned char)*s;
        if (c == '"' || c == '\\') { fputc('\\', out); fputc(c, out); }
        else if (c < 0x20) fprintf(out, "\\u%04x", c);
        else fputc(c, out);
    }
    fputc('"', out);
}

static void emit_json_hit(FILE* out, const Index* idx, const HitRecord* h) {
    char strand_buf[2] = { h->strand ? (char)h->strand : '.', 0 };
    fputc('{', out);
    fputs("\"target\":", out);   emit_json_str(out, idx->pool + h->target_off);
    fputs(",\"source_label\":", out); emit_json_str(out, idx->pool + h->source_off);
    fprintf(out, ",\"source_type\":\"%s\"", src_type_str(h->source_type));
    fprintf(out, ",\"identity\":%.4f", (double)h->identity);
    fprintf(out, ",\"bitscore\":%.2f", (double)h->bitscore);
    fprintf(out, ",\"evalue\":%.3g", h->evalue);
    fprintf(out, ",\"qstart\":%d,\"qend\":%d,\"tstart\":%d,\"tend\":%d",
            h->qstart, h->qend, h->tstart, h->tend);
    fprintf(out, ",\"alnlen\":%u,\"n_segments\":%u", h->alnlen, (unsigned)h->n_segments);
    fputs(",\"strand\":", out);  emit_json_str(out, strand_buf);
    fputc('}', out);
}

static void print_usage(void) {
    fprintf(stderr,
        "homolog_index_query — fast gene→hits lookup over a homolog index.\n"
        "\n"
        "Usage:\n"
        "  homolog_index_query --index <f.holindx> --gene <id>\n"
        "                      [--format tsv|json] (default: tsv)\n"
        "                      [--source diamond|miniprot]\n"
        "                      [--min_identity F] [--min_bitscore F]\n"
        "                      [--limit N] [--no_header]\n");
}

int main(int argc, char** argv) {
    const char* idx_path = NULL;
    const char* gene = NULL;
    const char* format = "tsv";
    const char* source_filter = NULL;
    double min_identity = -1, min_bitscore = -INFINITY;
    int limit = -1;
    int with_header = 1;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--index") && i+1<argc)        idx_path = argv[++i];
        else if (!strcmp(argv[i], "--gene")  && i+1<argc)        gene = argv[++i];
        else if (!strcmp(argv[i], "--format") && i+1<argc)       format = argv[++i];
        else if (!strcmp(argv[i], "--source") && i+1<argc)       source_filter = argv[++i];
        else if (!strcmp(argv[i], "--min_identity") && i+1<argc) min_identity = atof(argv[++i]);
        else if (!strcmp(argv[i], "--min_bitscore") && i+1<argc) min_bitscore = atof(argv[++i]);
        else if (!strcmp(argv[i], "--limit") && i+1<argc)        limit = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--no_header"))                with_header = 0;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[hiq] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }

    if (!idx_path || !gene) { print_usage(); return 1; }

    int wanted_source = -1;
    if (source_filter) {
        if      (!strcmp(source_filter, "diamond"))  wanted_source = HOLINDX_SRC_DIAMOND;
        else if (!strcmp(source_filter, "miniprot")) wanted_source = HOLINDX_SRC_MINIPROT;
        else { fprintf(stderr, "[hiq] Unknown --source '%s'\n", source_filter); return 1; }
    }

    Index idx; memset(&idx, 0, sizeof(idx)); idx.fd = -1;
    int rc = idx_open(idx_path, &idx);
    if (rc < 0) { fprintf(stderr, "[hiq] Cannot open index %s (rc=%d)\n", idx_path, rc); return 2; }

    int is_json = !strcmp(format, "json");
    const GeneRecord* g = idx_find(&idx, gene);

    if (!g) {
        if (is_json) {
            fputs("{\"query\":", stdout);
            emit_json_str(stdout, gene);
            fputs(",\"found\":false,\"n_hits\":0,\"hits\":[]}\n", stdout);
        } else if (with_header) {
            emit_tsv_header(stdout);
        }
        idx_close(&idx);
        return 0;
    }

    if (is_json) {
        fputs("{\"query\":", stdout);
        emit_json_str(stdout, idx.pool + g->name_off);
        fprintf(stdout, ",\"found\":true,\"n_hits_total\":%u,\"hits\":[", g->n_hits);
    } else if (with_header) {
        emit_tsv_header(stdout);
    }

    int emitted = 0;
    for (uint32_t k = 0; k < g->n_hits; k++) {
        const HitRecord* h = &idx.hits[g->hits_off + k];
        if (wanted_source >= 0 && h->source_type != (uint8_t)wanted_source) continue;
        if (min_identity >= 0 && h->identity < (float)min_identity) continue;
        if (h->bitscore < (float)min_bitscore) continue;
        if (limit > 0 && emitted >= limit) break;
        if (is_json) {
            if (emitted > 0) fputc(',', stdout);
            emit_json_hit(stdout, &idx, h);
        } else {
            emit_tsv_hit(stdout, &idx, g, h);
        }
        emitted++;
    }

    if (is_json) fprintf(stdout, "],\"n_hits_emitted\":%d}\n", emitted);

    idx_close(&idx);
    return 0;
}
