// =============================================================================
// homolog_atlas_server.c — tiny single-binary HTTP server over a homolog index.
//
// Mmaps a .holindx (built by homolog_index_build) once, then serves
// microsecond-latency lookups over plain HTTP. No external deps — pure C +
// libc + POSIX sockets + pthread. Drop into the atlas: the page calls
// fetch("/lookup?gene=GENE_X").
//
// Endpoints (all return application/json with CORS *):
//   GET /lookup?gene=<id>[&source=diamond|miniprot]
//              [&min_identity=F][&min_bitscore=F][&limit=N]
//   GET /health
//   GET /
//
// Concurrency: thread-per-connection (pthread_detach). Mmap'd index is
// read-only so all threads share it for free.
//
// Compile:
//   gcc -O3 -march=native -o homolog_atlas_server homolog_atlas_server.c -lpthread
// =============================================================================

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <pthread.h>
#include <signal.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "homolog_index.h"

#define SERVER_NAME "homolog_atlas_server/1"

// ── Shared mmap'd index ─────────────────────────────────────────────────────

static struct {
    int fd;
    void* base;
    size_t size;
    const IndexHeader* hdr;
    const GeneRecord* genes;
    const HitRecord* hits;
    const char* pool;
    const char* path;
} g_idx;

static volatile sig_atomic_t g_running = 1;
static void sigint_handler(int sig) { (void)sig; g_running = 0; }

static int idx_open(const char* path) {
    g_idx.fd = open(path, O_RDONLY);
    if (g_idx.fd < 0) return -1;
    struct stat st;
    if (fstat(g_idx.fd, &st) < 0) { close(g_idx.fd); return -2; }
    g_idx.size = (size_t)st.st_size;
    g_idx.base = mmap(NULL, g_idx.size, PROT_READ, MAP_PRIVATE, g_idx.fd, 0);
    if (g_idx.base == MAP_FAILED) { close(g_idx.fd); return -3; }
    g_idx.hdr = (const IndexHeader*)g_idx.base;
    if (g_idx.hdr->magic != HOLINDX_MAGIC)     { munmap(g_idx.base, g_idx.size); close(g_idx.fd); return -4; }
    if (g_idx.hdr->version != HOLINDX_VERSION) { munmap(g_idx.base, g_idx.size); close(g_idx.fd); return -5; }
    g_idx.genes = (const GeneRecord*)((const char*)g_idx.base + g_idx.hdr->off_gene_index);
    g_idx.hits  = (const HitRecord*) ((const char*)g_idx.base + g_idx.hdr->off_hits);
    g_idx.pool  = (const char*)g_idx.base + g_idx.hdr->off_string_pool;
    g_idx.path  = path;
    return 0;
}

static const GeneRecord* idx_find(const char* name) {
    int lo = 0, hi = (int)g_idx.hdr->n_genes - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        int cmp = strcmp(g_idx.pool + g_idx.genes[mid].name_off, name);
        if (cmp == 0) return &g_idx.genes[mid];
        if (cmp < 0) lo = mid + 1;
        else         hi = mid - 1;
    }
    return NULL;
}

static const char* src_str(uint8_t s) {
    return s == HOLINDX_SRC_DIAMOND ? "diamond" : (s == HOLINDX_SRC_MINIPROT ? "miniprot" : "unknown");
}

// ── Growable string buffer ──────────────────────────────────────────────────

typedef struct { char* buf; size_t len; size_t cap; } StrBuf;

static void sb_init(StrBuf* s) { s->cap = 4096; s->buf = (char*)malloc(s->cap); s->len = 0; s->buf[0] = 0; }
static void sb_free(StrBuf* s) { free(s->buf); s->buf = NULL; }
static void sb_reserve(StrBuf* s, size_t need) {
    while (s->len + need + 1 > s->cap) { s->cap *= 2; s->buf = (char*)realloc(s->buf, s->cap); }
}
static void sb_putc(StrBuf* s, char c) { sb_reserve(s, 1); s->buf[s->len++] = c; s->buf[s->len] = 0; }
static void sb_puts(StrBuf* s, const char* str) {
    size_t n = strlen(str); sb_reserve(s, n);
    memcpy(s->buf + s->len, str, n); s->len += n; s->buf[s->len] = 0;
}
static void sb_printf(StrBuf* s, const char* fmt, ...) {
    va_list ap; va_start(ap, fmt);
    va_list ap2; va_copy(ap2, ap);
    int need = vsnprintf(NULL, 0, fmt, ap);
    va_end(ap);
    if (need < 0) { va_end(ap2); return; }
    sb_reserve(s, (size_t)need);
    vsnprintf(s->buf + s->len, s->cap - s->len, fmt, ap2);
    va_end(ap2);
    s->len += (size_t)need;
}
static void sb_json_str(StrBuf* s, const char* str) {
    sb_putc(s, '"');
    for (; *str; str++) {
        unsigned char c = (unsigned char)*str;
        if (c == '"' || c == '\\') { sb_putc(s, '\\'); sb_putc(s, (char)c); }
        else if (c < 0x20)         sb_printf(s, "\\u%04x", c);
        else                       sb_putc(s, (char)c);
    }
    sb_putc(s, '"');
}

// ── URL decode ──────────────────────────────────────────────────────────────

static int hexval(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}
static void url_decode(char* s) {
    char* o = s;
    while (*s) {
        if (*s == '%' && s[1] && s[2]) {
            int hi = hexval(s[1]), lo = hexval(s[2]);
            if (hi >= 0 && lo >= 0) { *o++ = (char)((hi << 4) | lo); s += 3; continue; }
        }
        if (*s == '+') { *o++ = ' '; s++; continue; }
        *o++ = *s++;
    }
    *o = 0;
}

// ── Query string parsing ────────────────────────────────────────────────────

typedef struct {
    char gene[256];
    char source[32];
    int    has_min_identity;
    double min_identity;
    int    has_min_bitscore;
    double min_bitscore;
    int    has_limit;
    int    limit;
} QueryParams;

static void qp_parse(QueryParams* q, char* qs) {
    memset(q, 0, sizeof(*q));
    char* p = qs;
    while (p && *p) {
        char* amp = strchr(p, '&');
        if (amp) *amp = 0;
        char* eq = strchr(p, '=');
        if (eq) {
            *eq = 0;
            char* k = p;
            char* v = eq + 1;
            url_decode(v);
            if      (!strcmp(k, "gene"))         { strncpy(q->gene, v, sizeof(q->gene) - 1); }
            else if (!strcmp(k, "source"))       { strncpy(q->source, v, sizeof(q->source) - 1); }
            else if (!strcmp(k, "min_identity")) { q->has_min_identity = 1; q->min_identity = atof(v); }
            else if (!strcmp(k, "min_bitscore")) { q->has_min_bitscore = 1; q->min_bitscore = atof(v); }
            else if (!strcmp(k, "limit"))        { q->has_limit = 1; q->limit = atoi(v); }
        }
        if (!amp) break;
        p = amp + 1;
    }
}

// ── HTTP response writer ────────────────────────────────────────────────────

static const char* status_text(int code) {
    switch (code) {
        case 200: return "OK";
        case 400: return "Bad Request";
        case 404: return "Not Found";
        case 405: return "Method Not Allowed";
        case 500: return "Internal Server Error";
        default:  return "Error";
    }
}

static void send_full(int fd, const void* buf, size_t n) {
    const char* p = (const char*)buf;
    while (n > 0) {
        ssize_t w = write(fd, p, n);
        if (w <= 0) {
            if (w < 0 && errno == EINTR) continue;
            return;
        }
        p += w; n -= (size_t)w;
    }
}

static void send_response(int fd, int code, const char* ctype, const char* body, size_t blen) {
    char hdr[1024];
    int n = snprintf(hdr, sizeof(hdr),
                     "HTTP/1.1 %d %s\r\n"
                     "Server: %s\r\n"
                     "Content-Type: %s\r\n"
                     "Content-Length: %zu\r\n"
                     "Access-Control-Allow-Origin: *\r\n"
                     "Access-Control-Allow-Methods: GET, OPTIONS\r\n"
                     "Access-Control-Allow-Headers: *\r\n"
                     "Connection: close\r\n"
                     "\r\n",
                     code, status_text(code), SERVER_NAME, ctype, blen);
    if (n > 0) send_full(fd, hdr, (size_t)n);
    if (body && blen) send_full(fd, body, blen);
}

// ── Endpoint handlers ───────────────────────────────────────────────────────

static void emit_hit_json(StrBuf* s, const HitRecord* h) {
    char strand_buf[2] = { h->strand ? (char)h->strand : '.', 0 };
    sb_putc(s, '{');
    sb_puts(s, "\"target\":");        sb_json_str(s, g_idx.pool + h->target_off);
    sb_puts(s, ",\"source_label\":"); sb_json_str(s, g_idx.pool + h->source_off);
    sb_printf(s, ",\"source_type\":\"%s\"", src_str(h->source_type));
    sb_printf(s, ",\"identity\":%.4f,\"bitscore\":%.2f,\"evalue\":%.3g",
              (double)h->identity, (double)h->bitscore, h->evalue);
    sb_printf(s, ",\"qstart\":%d,\"qend\":%d,\"tstart\":%d,\"tend\":%d",
              h->qstart, h->qend, h->tstart, h->tend);
    sb_printf(s, ",\"alnlen\":%u,\"n_segments\":%u",
              h->alnlen, (unsigned)h->n_segments);
    sb_puts(s, ",\"strand\":");       sb_json_str(s, strand_buf);
    sb_putc(s, '}');
}

static void handle_lookup(int fd, QueryParams* q) {
    StrBuf body; sb_init(&body);

    if (q->gene[0] == 0) {
        sb_puts(&body, "{\"error\":\"missing required parameter: gene\"}");
        send_response(fd, 400, "application/json", body.buf, body.len);
        sb_free(&body); return;
    }

    const GeneRecord* g = idx_find(q->gene);
    if (!g) {
        sb_puts(&body, "{\"query\":"); sb_json_str(&body, q->gene);
        sb_puts(&body, ",\"found\":false,\"n_hits\":0,\"hits\":[]}");
        send_response(fd, 200, "application/json", body.buf, body.len);
        sb_free(&body); return;
    }

    int wanted_src = -1;
    if (q->source[0]) {
        if      (!strcmp(q->source, "diamond"))  wanted_src = HOLINDX_SRC_DIAMOND;
        else if (!strcmp(q->source, "miniprot")) wanted_src = HOLINDX_SRC_MINIPROT;
        else {
            sb_puts(&body, "{\"error\":\"source must be 'diamond' or 'miniprot'\"}");
            send_response(fd, 400, "application/json", body.buf, body.len);
            sb_free(&body); return;
        }
    }

    sb_puts(&body, "{\"query\":");
    sb_json_str(&body, g_idx.pool + g->name_off);
    sb_printf(&body, ",\"found\":true,\"n_hits_total\":%u,\"hits\":[", g->n_hits);

    int emitted = 0, first = 1;
    for (uint32_t k = 0; k < g->n_hits; k++) {
        const HitRecord* h = &g_idx.hits[g->hits_off + k];
        if (wanted_src >= 0 && h->source_type != (uint8_t)wanted_src) continue;
        if (q->has_min_identity && (double)h->identity < q->min_identity) continue;
        if (q->has_min_bitscore && (double)h->bitscore < q->min_bitscore) continue;
        if (q->has_limit && emitted >= q->limit) break;
        if (!first) sb_putc(&body, ',');
        emit_hit_json(&body, h);
        first = 0; emitted++;
    }
    sb_printf(&body, "],\"n_hits_emitted\":%d}", emitted);

    send_response(fd, 200, "application/json", body.buf, body.len);
    sb_free(&body);
}

static void handle_health(int fd) {
    StrBuf body; sb_init(&body);
    sb_puts(&body, "{\"ok\":true,\"service\":\"" SERVER_NAME "\",\"index_path\":");
    sb_json_str(&body, g_idx.path);
    sb_printf(&body,
              ",\"index_version\":%u,\"n_genes\":%u,\"n_hits\":%llu,"
              "\"size_mb\":%.3f,\"string_pool_bytes\":%llu}",
              g_idx.hdr->version, g_idx.hdr->n_genes,
              (unsigned long long)g_idx.hdr->n_hits,
              (double)g_idx.size / (1024.0 * 1024.0),
              (unsigned long long)g_idx.hdr->string_pool_size);
    send_response(fd, 200, "application/json", body.buf, body.len);
    sb_free(&body);
}

static void handle_root(int fd) {
    static const char* msg =
        "{\"service\":\"" SERVER_NAME "\","
        "\"endpoints\":["
        "\"GET /lookup?gene=<id>[&source=diamond|miniprot][&min_identity=F][&min_bitscore=F][&limit=N]\","
        "\"GET /health\","
        "\"GET /\""
        "]}";
    send_response(fd, 200, "application/json", msg, strlen(msg));
}

// ── Connection worker ───────────────────────────────────────────────────────

static void* handle_conn(void* arg) {
    int fd = (int)(intptr_t)arg;
    char buf[8192];
    ssize_t got = read(fd, buf, sizeof(buf) - 1);
    if (got <= 0) { close(fd); return NULL; }
    buf[got] = 0;

    char method[16] = {0}, path[1024] = {0};
    if (sscanf(buf, "%15s %1023s", method, path) < 2) {
        send_response(fd, 400, "application/json", "{\"error\":\"malformed request\"}", 30);
        close(fd); return NULL;
    }

    if (!strcmp(method, "OPTIONS")) {
        send_response(fd, 200, "application/json", "{}", 2);
        close(fd); return NULL;
    }
    if (strcmp(method, "GET") != 0) {
        send_response(fd, 405, "application/json", "{\"error\":\"only GET supported\"}", 31);
        close(fd); return NULL;
    }

    char* qs = strchr(path, '?');
    if (qs) *qs++ = 0;

    if      (!strcmp(path, "/lookup")) { QueryParams q; qp_parse(&q, qs); handle_lookup(fd, &q); }
    else if (!strcmp(path, "/health")) handle_health(fd);
    else if (!strcmp(path, "/"))       handle_root(fd);
    else                                send_response(fd, 404, "application/json", "{\"error\":\"not found\"}", 21);

    close(fd);
    return NULL;
}

// ── Main ────────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "homolog_atlas_server — tiny HTTP server over a homolog index.\n"
        "\n"
        "  --index <file.holindx>   Required.\n"
        "  --port N                 Listen port (default 8765).\n"
        "  --bind <addr>            Bind address (default 127.0.0.1; use 0.0.0.0 for any).\n"
        "\n"
        "Endpoints (all return application/json, CORS *):\n"
        "  GET /lookup?gene=<id>[&source=diamond|miniprot][&min_identity=F][&min_bitscore=F][&limit=N]\n"
        "  GET /health\n"
        "  GET /\n");
}

int main(int argc, char** argv) {
    const char* index_path = NULL;
    int port = 8765;
    const char* bind_addr = "127.0.0.1";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--index") && i + 1 < argc) index_path = argv[++i];
        else if (!strcmp(argv[i], "--port")  && i + 1 < argc) port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--bind")  && i + 1 < argc) bind_addr = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[hsrv] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }
    if (!index_path) { fprintf(stderr, "[hsrv] --index required\n"); return 1; }

    int rc = idx_open(index_path);
    if (rc < 0) { fprintf(stderr, "[hsrv] Cannot open index %s (rc=%d)\n", index_path, rc); return 2; }
    fprintf(stderr, "[hsrv] Index: %s (%u genes, %llu hits, %.2f MB)\n",
            index_path, g_idx.hdr->n_genes,
            (unsigned long long)g_idx.hdr->n_hits,
            (double)g_idx.size / (1024.0 * 1024.0));

    signal(SIGINT,  sigint_handler);
    signal(SIGTERM, sigint_handler);
    signal(SIGPIPE, SIG_IGN);

    int srv = socket(AF_INET, SOCK_STREAM, 0);
    if (srv < 0) { perror("socket"); return 3; }
    int one = 1;
    setsockopt(srv, SOL_SOCKET, SO_REUSEADDR, &one, sizeof(one));

    struct sockaddr_in addr;
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port = htons((uint16_t)port);
    if (inet_pton(AF_INET, bind_addr, &addr.sin_addr) != 1) {
        fprintf(stderr, "[hsrv] Invalid --bind: %s\n", bind_addr); return 4;
    }

    if (bind(srv, (struct sockaddr*)&addr, sizeof(addr)) < 0) { perror("bind"); return 5; }
    if (listen(srv, 64) < 0)                                   { perror("listen"); return 6; }

    fprintf(stderr, "[hsrv] Listening on http://%s:%d  (PID %d)\n", bind_addr, port, getpid());

    while (g_running) {
        struct sockaddr_in caddr;
        socklen_t clen = sizeof(caddr);
        int cfd = accept(srv, (struct sockaddr*)&caddr, &clen);
        if (cfd < 0) {
            if (errno == EINTR) continue;
            if (!g_running) break;
            perror("accept"); continue;
        }
        pthread_t tid;
        if (pthread_create(&tid, NULL, handle_conn, (void*)(intptr_t)cfd) != 0) {
            close(cfd); continue;
        }
        pthread_detach(tid);
    }

    fprintf(stderr, "[hsrv] Shutting down\n");
    close(srv);
    munmap(g_idx.base, g_idx.size);
    close(g_idx.fd);
    return 0;
}
