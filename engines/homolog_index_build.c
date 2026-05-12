// =============================================================================
// homolog_index_build.c — build a sorted, mmap-friendly homolog index from
//                         DIAMOND tabular and/or miniprot PAF files.
//
// Goal: a single offline build → microsecond gene→hits lookups via
//       homolog_index_query.
//
// Inputs (any combination):
//   --diamond <path>[=<label>]   DIAMOND -outfmt 6 (12 cols: qseqid sseqid pident
//                                length mismatch gapopen qstart qend sstart send
//                                evalue bitscore).
//   --miniprot <path>[=<label>]  miniprot PAF (standard 12 cols + optional tags;
//                                ms:i: or AS:i: parsed as bitscore if present).
//   --manifest <file.tsv>        TSV (optional header): type<TAB>path<TAB>label
//
// Output:
//   --out <file.holindx>         Binary index (see engines/homolog_index.h).
//
// Compile: gcc -O3 -march=native -o homolog_index_build homolog_index_build.c
// =============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include <ctype.h>

#include "homolog_index.h"

// ── String pool ─────────────────────────────────────────────────────────────

typedef struct { char* buf; size_t len; size_t cap; } StrPool;

static void sp_init(StrPool* p) { p->cap = 1<<16; p->buf = (char*)malloc(p->cap); p->len = 0; }
static uint32_t sp_append(StrPool* p, const char* s) {
    size_t n = strlen(s) + 1;
    while (p->len + n > p->cap) { p->cap *= 2; p->buf = (char*)realloc(p->buf, p->cap); }
    uint32_t off = (uint32_t)p->len;
    memcpy(p->buf + p->len, s, n);
    p->len += n;
    return off;
}

// ── Hash map (FNV-1a + open addressing) string → pool offset ────────────────

typedef struct { uint32_t hash; uint32_t pool_off; } HashEntry;
typedef struct { HashEntry* entries; size_t cap; size_t n; StrPool* pool; } StrMap;

static uint32_t fnv1a(const char* s) {
    uint32_t h = 2166136261u;
    while (*s) { h ^= (uint8_t)*s++; h *= 16777619u; }
    return h;
}

static void sm_init(StrMap* m, StrPool* p) {
    m->cap = 1<<16;
    m->n = 0;
    m->pool = p;
    m->entries = (HashEntry*)malloc(m->cap * sizeof(HashEntry));
    for (size_t i = 0; i < m->cap; i++) m->entries[i].pool_off = 0xFFFFFFFFu;
}

static void sm_grow(StrMap* m) {
    size_t old_cap = m->cap;
    HashEntry* old = m->entries;
    m->cap *= 2;
    m->entries = (HashEntry*)malloc(m->cap * sizeof(HashEntry));
    for (size_t i = 0; i < m->cap; i++) m->entries[i].pool_off = 0xFFFFFFFFu;
    for (size_t i = 0; i < old_cap; i++) {
        if (old[i].pool_off != 0xFFFFFFFFu) {
            size_t mask = m->cap - 1;
            size_t pos = old[i].hash & mask;
            while (m->entries[pos].pool_off != 0xFFFFFFFFu) pos = (pos + 1) & mask;
            m->entries[pos] = old[i];
        }
    }
    free(old);
}

static uint32_t sm_intern(StrMap* m, const char* s) {
    uint32_t h = fnv1a(s);
    size_t mask = m->cap - 1;
    size_t pos = h & mask;
    while (m->entries[pos].pool_off != 0xFFFFFFFFu) {
        if (m->entries[pos].hash == h &&
            strcmp(m->pool->buf + m->entries[pos].pool_off, s) == 0) {
            return m->entries[pos].pool_off;
        }
        pos = (pos + 1) & mask;
    }
    uint32_t off = sp_append(m->pool, s);
    m->entries[pos].hash = h;
    m->entries[pos].pool_off = off;
    m->n++;
    if (m->n * 2 > m->cap) sm_grow(m);
    return off;
}

// ── In-memory hit list (built before sorting/dedup) ─────────────────────────

typedef struct {
    uint32_t query_off;
    uint32_t target_off;
    uint32_t source_off;
    int32_t  qstart, qend, tstart, tend;
    float    identity;
    float    bitscore;
    double   evalue;
    uint32_t alnlen;
    uint16_t n_segments;
    uint8_t  source_type;
    uint8_t  strand;
} BuildHit;

typedef struct { BuildHit* items; size_t n; size_t cap; } HitList;

static void hl_push(HitList* h, BuildHit b) {
    if (h->n == h->cap) {
        h->cap = h->cap ? h->cap * 2 : 1024;
        h->items = (BuildHit*)realloc(h->items, h->cap * sizeof(BuildHit));
    }
    h->items[h->n++] = b;
}

// ── Parsers ─────────────────────────────────────────────────────────────────

static int parse_diamond(const char* path, const char* source_label,
                          StrMap* m, HitList* h) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[hib] Cannot open %s\n", path); return 0; }
    uint32_t src_off = sm_intern(m, source_label);
    char* line = NULL;
    size_t cap = 0;
    ssize_t got;
    int n = 0;
    while ((got = getline(&line, &cap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        char q[256], t[256];
        double pident, evalue, bits;
        int alnlen, mism, gap, qs, qe, ss, se;
        int nf = sscanf(line, "%255s %255s %lf %d %d %d %d %d %d %d %lf %lf",
                        q, t, &pident, &alnlen, &mism, &gap,
                        &qs, &qe, &ss, &se, &evalue, &bits);
        if (nf < 12) continue;
        BuildHit b;
        memset(&b, 0, sizeof(b));
        b.query_off = sm_intern(m, q);
        b.target_off = sm_intern(m, t);
        b.source_off = src_off;
        b.qstart = qs; b.qend = qe;
        b.tstart = ss; b.tend = se;
        b.identity = (float)(pident / 100.0);
        b.bitscore = (float)bits;
        b.evalue = evalue;
        b.alnlen = (uint32_t)alnlen;
        b.n_segments = 1;
        b.source_type = HOLINDX_SRC_DIAMOND;
        b.strand = 0;
        hl_push(h, b);
        n++;
    }
    free(line);
    fclose(f);
    fprintf(stderr, "[hib] DIAMOND: %d hits from %s (label=%s)\n", n, path, source_label);
    return n;
}

static int parse_miniprot(const char* path, const char* source_label,
                          StrMap* m, HitList* h) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[hib] Cannot open %s\n", path); return 0; }
    uint32_t src_off = sm_intern(m, source_label);
    char* line = NULL;
    size_t cap = 0;
    ssize_t got;
    int n = 0;
    while ((got = getline(&line, &cap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        char q[256], t[256], strand_s[8];
        int qlen, qs, qe, tlen, ts, te, matches, alnlen, mapq;
        int nf = sscanf(line, "%255s %d %d %d %7s %255s %d %d %d %d %d %d",
                        q, &qlen, &qs, &qe, strand_s, t, &tlen, &ts, &te,
                        &matches, &alnlen, &mapq);
        if (nf < 12) continue;
        double bitscore = -1.0;
        const char* tag = strstr(line, "\tms:i:");
        if (tag) bitscore = atof(tag + 6);
        else if ((tag = strstr(line, "\tAS:i:"))) bitscore = atof(tag + 6);
        BuildHit b;
        memset(&b, 0, sizeof(b));
        b.query_off = sm_intern(m, q);
        b.target_off = sm_intern(m, t);
        b.source_off = src_off;
        b.qstart = qs; b.qend = qe;
        b.tstart = ts; b.tend = te;
        b.identity = (alnlen > 0) ? (float)matches / (float)alnlen : 0.0f;
        b.bitscore = (float)bitscore;
        b.evalue = -1.0;
        b.alnlen = (uint32_t)alnlen;
        b.n_segments = 1;
        b.source_type = HOLINDX_SRC_MINIPROT;
        b.strand = (strand_s[0] == '-') ? '-' : '+';
        hl_push(h, b);
        n++;
    }
    free(line);
    fclose(f);
    fprintf(stderr, "[hib] miniprot: %d hits from %s (label=%s)\n", n, path, source_label);
    return n;
}

// ── Sorting comparators (need a global pool pointer for strcmp by offset) ──

static StrPool* sort_pool = NULL;

static int hit_cmp(const void* a, const void* b) {
    const BuildHit* ha = (const BuildHit*)a;
    const BuildHit* hb = (const BuildHit*)b;
    int r = strcmp(sort_pool->buf + ha->query_off, sort_pool->buf + hb->query_off);
    if (r) return r;
    if (ha->bitscore > hb->bitscore) return -1;
    if (ha->bitscore < hb->bitscore) return 1;
    if (ha->identity > hb->identity) return -1;
    if (ha->identity < hb->identity) return 1;
    return 0;
}

static int gene_cmp(const void* a, const void* b) {
    const GeneRecord* ga = (const GeneRecord*)a;
    const GeneRecord* gb = (const GeneRecord*)b;
    return strcmp(sort_pool->buf + ga->name_off, sort_pool->buf + gb->name_off);
}

// ── Manifest parsing ────────────────────────────────────────────────────────

typedef struct { int kind; char path[512]; char label[128]; } SrcEntry;

static int parse_manifest(const char* path, SrcEntry** out, int* n_out) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[hib] Cannot open manifest %s\n", path); return 0; }
    int cap = 16, n = 0;
    SrcEntry* arr = (SrcEntry*)malloc(cap * sizeof(SrcEntry));
    char* line = NULL; size_t lcap = 0; ssize_t got;
    int header_done = 0;
    while ((got = getline(&line, &lcap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        char type[32], p[512], lab[128];
        int nf = sscanf(line, "%31s %511s %127s", type, p, lab);
        if (nf < 3) continue;
        if (!header_done && (!strcasecmp(type, "type") || !strcasecmp(type, "kind"))) {
            header_done = 1; continue;
        }
        header_done = 1;
        if (n == cap) { cap *= 2; arr = (SrcEntry*)realloc(arr, cap * sizeof(SrcEntry)); }
        if (!strcasecmp(type, "diamond"))       arr[n].kind = HOLINDX_SRC_DIAMOND;
        else if (!strcasecmp(type, "miniprot")) arr[n].kind = HOLINDX_SRC_MINIPROT;
        else { fprintf(stderr, "[hib] Unknown source type '%s' in manifest; skipping\n", type); continue; }
        strncpy(arr[n].path, p, sizeof(arr[n].path) - 1);
        arr[n].path[sizeof(arr[n].path) - 1] = 0;
        strncpy(arr[n].label, lab, sizeof(arr[n].label) - 1);
        arr[n].label[sizeof(arr[n].label) - 1] = 0;
        n++;
    }
    free(line);
    fclose(f);
    *out = arr;
    *n_out = n;
    return n;
}

// ── Main ────────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "homolog_index_build — build a sorted, mmap-friendly homolog index.\n"
        "\n"
        "Inputs (any combination, may repeat):\n"
        "  --diamond <path>[=<label>]   DIAMOND tabular -outfmt 6\n"
        "  --miniprot <path>[=<label>]  miniprot PAF\n"
        "  --manifest <file.tsv>        TSV: type<TAB>path<TAB>label\n"
        "\n"
        "Output:\n"
        "  --out <file.holindx>         Output binary index (required)\n"
        "\n"
        "Label defaults to basename of the path. The label is what shows up\n"
        "as `source_label` in query output (e.g. species name).\n");
}

static const char* basename_of(const char* p) {
    const char* s = strrchr(p, '/');
    return s ? s + 1 : p;
}

int main(int argc, char** argv) {
    const char* out_path = NULL;
    SrcEntry* mani = NULL; int mani_n = 0;
    SrcEntry* srcs = NULL; int n_srcs = 0, cap_srcs = 0;

    for (int i = 1; i < argc; i++) {
        if ((!strcmp(argv[i], "--diamond") || !strcmp(argv[i], "--miniprot")) && i + 1 < argc) {
            int kind = !strcmp(argv[i], "--diamond") ? HOLINDX_SRC_DIAMOND : HOLINDX_SRC_MINIPROT;
            const char* spec = argv[++i];
            const char* eq = strchr(spec, '=');
            char path[512], label[128];
            if (eq) {
                size_t pl = (size_t)(eq - spec);
                if (pl >= sizeof(path)) pl = sizeof(path) - 1;
                memcpy(path, spec, pl); path[pl] = 0;
                strncpy(label, eq + 1, sizeof(label) - 1); label[sizeof(label) - 1] = 0;
            } else {
                strncpy(path, spec, sizeof(path) - 1); path[sizeof(path) - 1] = 0;
                strncpy(label, basename_of(spec), sizeof(label) - 1); label[sizeof(label) - 1] = 0;
            }
            if (n_srcs == cap_srcs) {
                cap_srcs = cap_srcs ? cap_srcs * 2 : 16;
                srcs = (SrcEntry*)realloc(srcs, cap_srcs * sizeof(SrcEntry));
            }
            srcs[n_srcs].kind = kind;
            strncpy(srcs[n_srcs].path, path, sizeof(srcs[n_srcs].path) - 1);
            srcs[n_srcs].path[sizeof(srcs[n_srcs].path) - 1] = 0;
            strncpy(srcs[n_srcs].label, label, sizeof(srcs[n_srcs].label) - 1);
            srcs[n_srcs].label[sizeof(srcs[n_srcs].label) - 1] = 0;
            n_srcs++;
        }
        else if (!strcmp(argv[i], "--manifest") && i + 1 < argc) parse_manifest(argv[++i], &mani, &mani_n);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[hib] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }

    if (!out_path) { print_usage(); return 1; }
    if (n_srcs == 0 && mani_n == 0) {
        fprintf(stderr, "[hib] No sources provided.\n"); print_usage(); return 1;
    }

    StrPool pool; sp_init(&pool);
    sp_append(&pool, ""); // reserve offset 0 for the empty string
    StrMap map; sm_init(&map, &pool);
    sm_intern(&map, "");
    HitList hits = {0};

    for (int i = 0; i < n_srcs; i++) {
        if (srcs[i].kind == HOLINDX_SRC_DIAMOND) parse_diamond(srcs[i].path, srcs[i].label, &map, &hits);
        else                                      parse_miniprot(srcs[i].path, srcs[i].label, &map, &hits);
    }
    for (int i = 0; i < mani_n; i++) {
        if (mani[i].kind == HOLINDX_SRC_DIAMOND) parse_diamond(mani[i].path, mani[i].label, &map, &hits);
        else                                      parse_miniprot(mani[i].path, mani[i].label, &map, &hits);
    }

    fprintf(stderr, "[hib] Total: %zu hits, %zu pool bytes, %zu unique strings\n",
            hits.n, pool.len, map.n);
    if (hits.n == 0) { fprintf(stderr, "[hib] No hits parsed; aborting\n"); return 1; }

    sort_pool = &pool;
    qsort(hits.items, hits.n, sizeof(BuildHit), hit_cmp);

    // Collapse contiguous runs of the same query into GeneRecords (in hit-array order)
    GeneRecord* genes = NULL;
    size_t n_genes = 0, cap_genes = 0;
    size_t i = 0;
    while (i < hits.n) {
        size_t j = i;
        uint32_t qoff = hits.items[i].query_off;
        while (j < hits.n && hits.items[j].query_off == qoff) j++;
        if (n_genes == cap_genes) {
            cap_genes = cap_genes ? cap_genes * 2 : 1024;
            genes = (GeneRecord*)realloc(genes, cap_genes * sizeof(GeneRecord));
        }
        memset(&genes[n_genes], 0, sizeof(GeneRecord));
        genes[n_genes].name_off = qoff;
        genes[n_genes].name_len = (uint32_t)strlen(pool.buf + qoff);
        genes[n_genes].hits_off = i;
        genes[n_genes].n_hits = (uint32_t)(j - i);
        n_genes++;
        i = j;
    }
    fprintf(stderr, "[hib] %zu genes indexed\n", n_genes);

    // Sort gene records by name (alphabetical) for binary search.
    // The hits_off values remain valid because each gene still references a contiguous
    // range in the (alphabetically-by-query) hit array.
    qsort(genes, n_genes, sizeof(GeneRecord), gene_cmp);

    FILE* f = fopen(out_path, "wb");
    if (!f) { fprintf(stderr, "[hib] Cannot open %s for writing\n", out_path); return 2; }

    IndexHeader hdr;
    memset(&hdr, 0, sizeof(hdr));
    hdr.magic = HOLINDX_MAGIC;
    hdr.version = HOLINDX_VERSION;
    hdr.n_genes = (uint32_t)n_genes;
    hdr.n_hits = hits.n;
    hdr.string_pool_size = pool.len;
    hdr.header_size = sizeof(IndexHeader);
    hdr.gene_record_size = sizeof(GeneRecord);
    hdr.hit_record_size = sizeof(HitRecord);
    hdr.off_gene_index = sizeof(IndexHeader);
    hdr.off_hits = hdr.off_gene_index + (uint64_t)n_genes * sizeof(GeneRecord);
    hdr.off_string_pool = hdr.off_hits + (uint64_t)hits.n * sizeof(HitRecord);

    fwrite(&hdr, sizeof(hdr), 1, f);
    fwrite(genes, sizeof(GeneRecord), n_genes, f);

    for (size_t k = 0; k < hits.n; k++) {
        HitRecord hr;
        memset(&hr, 0, sizeof(hr));
        hr.target_off = hits.items[k].target_off;
        hr.source_off = hits.items[k].source_off;
        hr.qstart = hits.items[k].qstart;
        hr.qend = hits.items[k].qend;
        hr.tstart = hits.items[k].tstart;
        hr.tend = hits.items[k].tend;
        hr.identity = hits.items[k].identity;
        hr.bitscore = hits.items[k].bitscore;
        hr.evalue = hits.items[k].evalue;
        hr.alnlen = hits.items[k].alnlen;
        hr.n_segments = hits.items[k].n_segments;
        hr.source_type = hits.items[k].source_type;
        hr.strand = hits.items[k].strand;
        fwrite(&hr, sizeof(hr), 1, f);
    }
    fwrite(pool.buf, 1, pool.len, f);
    fclose(f);

    fprintf(stderr, "[hib] Wrote %s (%zu genes, %zu hits, %zu pool bytes)\n",
            out_path, n_genes, hits.n, pool.len);

    free(srcs); free(mani);
    free(hits.items); free(genes);
    free(pool.buf); free(map.entries);
    return 0;
}
