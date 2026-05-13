// =============================================================================
// popstats_dispatch.c — JSON activator → engine → JSON extractor dispatcher.
//
// Single binary that fronts the popstats engines for the local popstats
// server. The HTTP server hands an activator JSON to this binary on stdin
// (or via --in <file>); the dispatcher:
//
//   1. Parses the JSON activator.
//   2. Routes on the `.statistic` field to the matching per-statistic
//      handler (xpehh / iHS / outlier_scan / candidate_vs_flanks).
//   3. Maps activator JSON fields to engine CLI flags using the
//      engine_binding contract documented in
//      engines/schemas/<stat>.activator.schema.json.
//   4. Spawns the engine binary, captures TSV output to a temp file.
//   5. Parses the TSV and emits the corresponding extractor JSON to
//      stdout (or --out <file>), matching
//      engines/schemas/<stat>.extractor.schema.json.
//
// The engines themselves stay TSV-emitting plumbing. The atlas client only
// ever sees JSON — this binary is the JSON↔TSV bridge.
//
// JSON parser: minimal hand-rolled recursive descent, handles the subset
// our activator schemas use (object/array/string/number/boolean/null,
// no Unicode escapes beyond ASCII).
//
// Compile: gcc -O3 -Wall -o popstats_dispatch popstats_dispatch.c
// =============================================================================

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <errno.h>

#define DISPATCH_VERSION "popstats_dispatch_v1"

// ── Minimal JSON parser ─────────────────────────────────────────────────────

typedef enum {
    JV_NULL = 0, JV_BOOL, JV_NUMBER, JV_STRING, JV_OBJECT, JV_ARRAY
} JVType;

typedef struct JVNode JVNode;
struct JVNode {
    JVType type;
    int    b;          // BOOL
    double n;          // NUMBER
    char*  s;          // STRING
    char** keys;       // OBJECT keys
    JVNode** vals;     // OBJECT values OR ARRAY items
    int    count;      // OBJECT n_keys OR ARRAY n_items
};

static void jv_free(JVNode* v) {
    if (!v) return;
    if (v->type == JV_STRING) free(v->s);
    if (v->type == JV_OBJECT) {
        for (int i = 0; i < v->count; i++) { free(v->keys[i]); jv_free(v->vals[i]); }
        free(v->keys); free(v->vals);
    } else if (v->type == JV_ARRAY) {
        for (int i = 0; i < v->count; i++) jv_free(v->vals[i]);
        free(v->vals);
    }
    free(v);
}

typedef struct {
    const char* src;
    size_t pos;
    size_t len;
    char   err[256];
} JParser;

static void jp_skip_ws(JParser* p) {
    while (p->pos < p->len) {
        char c = p->src[p->pos];
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') p->pos++;
        else break;
    }
}

static int jp_match(JParser* p, char c) {
    jp_skip_ws(p);
    if (p->pos < p->len && p->src[p->pos] == c) { p->pos++; return 1; }
    return 0;
}

static char* jp_parse_string(JParser* p) {
    jp_skip_ws(p);
    if (p->pos >= p->len || p->src[p->pos] != '"') {
        snprintf(p->err, sizeof(p->err), "expected '\"' at %zu", p->pos); return NULL;
    }
    p->pos++;
    size_t start = p->pos;
    size_t cap = 64, n = 0;
    char* out = (char*)malloc(cap);
    while (p->pos < p->len) {
        char c = p->src[p->pos++];
        if (c == '"') { (void)start; out[n] = 0; return out; }
        if (c == '\\' && p->pos < p->len) {
            char e = p->src[p->pos++];
            char dec = 0;
            switch (e) {
                case '"': dec = '"'; break;
                case '\\': dec = '\\'; break;
                case '/': dec = '/'; break;
                case 'n': dec = '\n'; break;
                case 't': dec = '\t'; break;
                case 'r': dec = '\r'; break;
                case 'b': dec = '\b'; break;
                case 'f': dec = '\f'; break;
                case 'u': {
                    if (p->pos + 4 > p->len) { free(out); return NULL; }
                    unsigned u = 0;
                    for (int i = 0; i < 4; i++) {
                        char h = p->src[p->pos++];
                        u <<= 4;
                        if (h >= '0' && h <= '9') u |= h - '0';
                        else if (h >= 'a' && h <= 'f') u |= h - 'a' + 10;
                        else if (h >= 'A' && h <= 'F') u |= h - 'A' + 10;
                    }
                    if (u < 0x80) dec = (char)u;
                    else { /* simplistic: replace with '?' */ dec = '?'; }
                    break;
                }
                default: dec = e;
            }
            if (n + 1 >= cap) { cap *= 2; out = (char*)realloc(out, cap); }
            out[n++] = dec;
        } else {
            if (n + 1 >= cap) { cap *= 2; out = (char*)realloc(out, cap); }
            out[n++] = c;
        }
    }
    free(out);
    snprintf(p->err, sizeof(p->err), "unterminated string at %zu", p->pos);
    return NULL;
}

static JVNode* jp_parse_value(JParser* p);

static JVNode* jp_parse_object(JParser* p) {
    if (!jp_match(p, '{')) { snprintf(p->err, sizeof(p->err), "expected '{'"); return NULL; }
    JVNode* obj = (JVNode*)calloc(1, sizeof(JVNode));
    obj->type = JV_OBJECT;
    int cap = 8;
    obj->keys = (char**)malloc((size_t)cap * sizeof(char*));
    obj->vals = (JVNode**)malloc((size_t)cap * sizeof(JVNode*));
    jp_skip_ws(p);
    if (jp_match(p, '}')) return obj;
    while (1) {
        char* k = jp_parse_string(p);
        if (!k) { jv_free(obj); return NULL; }
        if (!jp_match(p, ':')) { free(k); jv_free(obj); snprintf(p->err, sizeof(p->err), "expected ':'"); return NULL; }
        JVNode* v = jp_parse_value(p);
        if (!v) { free(k); jv_free(obj); return NULL; }
        if (obj->count == cap) {
            cap *= 2;
            obj->keys = (char**)realloc(obj->keys, (size_t)cap * sizeof(char*));
            obj->vals = (JVNode**)realloc(obj->vals, (size_t)cap * sizeof(JVNode*));
        }
        obj->keys[obj->count] = k;
        obj->vals[obj->count] = v;
        obj->count++;
        if (jp_match(p, ',')) continue;
        if (jp_match(p, '}')) return obj;
        jv_free(obj);
        snprintf(p->err, sizeof(p->err), "expected ',' or '}'");
        return NULL;
    }
}

static JVNode* jp_parse_array(JParser* p) {
    if (!jp_match(p, '[')) { snprintf(p->err, sizeof(p->err), "expected '['"); return NULL; }
    JVNode* arr = (JVNode*)calloc(1, sizeof(JVNode));
    arr->type = JV_ARRAY;
    int cap = 8;
    arr->vals = (JVNode**)malloc((size_t)cap * sizeof(JVNode*));
    jp_skip_ws(p);
    if (jp_match(p, ']')) return arr;
    while (1) {
        JVNode* v = jp_parse_value(p);
        if (!v) { jv_free(arr); return NULL; }
        if (arr->count == cap) {
            cap *= 2;
            arr->vals = (JVNode**)realloc(arr->vals, (size_t)cap * sizeof(JVNode*));
        }
        arr->vals[arr->count++] = v;
        if (jp_match(p, ',')) continue;
        if (jp_match(p, ']')) return arr;
        jv_free(arr);
        snprintf(p->err, sizeof(p->err), "expected ',' or ']'");
        return NULL;
    }
}

static JVNode* jp_parse_value(JParser* p) {
    jp_skip_ws(p);
    if (p->pos >= p->len) { snprintf(p->err, sizeof(p->err), "unexpected EOF"); return NULL; }
    char c = p->src[p->pos];
    if (c == '{') return jp_parse_object(p);
    if (c == '[') return jp_parse_array(p);
    if (c == '"') {
        char* s = jp_parse_string(p);
        if (!s) return NULL;
        JVNode* v = (JVNode*)calloc(1, sizeof(JVNode));
        v->type = JV_STRING; v->s = s;
        return v;
    }
    if (c == 't' && p->pos + 4 <= p->len && !strncmp(p->src + p->pos, "true", 4)) {
        p->pos += 4;
        JVNode* v = (JVNode*)calloc(1, sizeof(JVNode));
        v->type = JV_BOOL; v->b = 1; return v;
    }
    if (c == 'f' && p->pos + 5 <= p->len && !strncmp(p->src + p->pos, "false", 5)) {
        p->pos += 5;
        JVNode* v = (JVNode*)calloc(1, sizeof(JVNode));
        v->type = JV_BOOL; v->b = 0; return v;
    }
    if (c == 'n' && p->pos + 4 <= p->len && !strncmp(p->src + p->pos, "null", 4)) {
        p->pos += 4;
        JVNode* v = (JVNode*)calloc(1, sizeof(JVNode));
        v->type = JV_NULL; return v;
    }
    // Number
    if (c == '-' || (c >= '0' && c <= '9')) {
        char* end;
        double d = strtod(p->src + p->pos, &end);
        if (end == p->src + p->pos) { snprintf(p->err, sizeof(p->err), "bad number"); return NULL; }
        p->pos = end - p->src;
        JVNode* v = (JVNode*)calloc(1, sizeof(JVNode));
        v->type = JV_NUMBER; v->n = d; return v;
    }
    snprintf(p->err, sizeof(p->err), "unexpected char '%c' at %zu", c, p->pos);
    return NULL;
}

static JVNode* jv_parse(const char* src, char* err_out, size_t err_sz) {
    JParser p = { .src = src, .pos = 0, .len = strlen(src), .err = {0} };
    JVNode* v = jp_parse_value(&p);
    if (!v) {
        if (err_out) snprintf(err_out, err_sz, "%s", p.err);
        return NULL;
    }
    return v;
}

static JVNode* jv_get(JVNode* obj, const char* key) {
    if (!obj || obj->type != JV_OBJECT) return NULL;
    for (int i = 0; i < obj->count; i++)
        if (!strcmp(obj->keys[i], key)) return obj->vals[i];
    return NULL;
}
static const char* jv_str(JVNode* v, const char* fallback) {
    return (v && v->type == JV_STRING) ? v->s : fallback;
}
static double jv_num(JVNode* v, double fallback) {
    return (v && v->type == JV_NUMBER) ? v->n : fallback;
}
static int jv_bool_(JVNode* v, int fallback) {
    return (v && v->type == JV_BOOL) ? v->b : fallback;
}
static int jv_arr_len(JVNode* v) {
    return (v && v->type == JV_ARRAY) ? v->count : 0;
}
static JVNode* jv_arr_at(JVNode* v, int i) {
    return (v && v->type == JV_ARRAY && i >= 0 && i < v->count) ? v->vals[i] : NULL;
}

// ── JSON emit helpers ──────────────────────────────────────────────────────

static void emit_jstr(FILE* f, const char* s) {
    fputc('"', f);
    for (const char* p = s ? s : ""; *p; p++) {
        unsigned char c = (unsigned char)*p;
        if (c == '"' || c == '\\') { fputc('\\', f); fputc(c, f); }
        else if (c < 0x20) fprintf(f, "\\u%04x", c);
        else fputc(c, f);
    }
    fputc('"', f);
}

static void emit_num_or_null(FILE* f, double v) {
    if (isnan(v) || !isfinite(v)) fputs("null", f);
    else fprintf(f, "%.6g", v);
}

// Parse a TSV cell value into a double; returns 0 on NA/empty/non-numeric.
static int cell_num(const char* s, double* out) {
    if (!s || !*s || !strcasecmp(s, "NA") || !strcasecmp(s, "nan")) return 0;
    if (!strcasecmp(s, "inf") || !strcasecmp(s, "-inf")) return 0;
    char* end; double v = strtod(s, &end);
    if (end == s || !isfinite(v)) return 0;
    *out = v; return 1;
}

// ── argv builder ──────────────────────────────────────────────────────────

typedef struct { char** argv; int n; int cap; } Argv;
static void argv_push(Argv* a, const char* s) {
    if (a->n == a->cap) {
        a->cap = a->cap ? a->cap * 2 : 16;
        a->argv = (char**)realloc(a->argv, (size_t)a->cap * sizeof(char*));
    }
    a->argv[a->n++] = s ? strdup(s) : NULL;
}
static void argv_push_int(Argv* a, long v)    { char buf[32]; snprintf(buf, sizeof(buf), "%ld", v); argv_push(a, buf); }
static void argv_push_double(Argv* a, double v){ char buf[32]; snprintf(buf, sizeof(buf), "%.6g", v); argv_push(a, buf); }
static void argv_free(Argv* a) {
    for (int i = 0; i < a->n; i++) free(a->argv[i]);
    free(a->argv); memset(a, 0, sizeof(*a));
}

// Append "--flag <value>" if v is non-NULL.
static void push_str_flag(Argv* a, const char* flag, const char* v) {
    if (v) { argv_push(a, flag); argv_push(a, v); }
}
static void push_str_flag_if(Argv* a, const char* flag, JVNode* v) {
    const char* s = jv_str(v, NULL);
    if (s) push_str_flag(a, flag, s);
}
static void push_num_flag_if(Argv* a, const char* flag, JVNode* v) {
    if (v && v->type == JV_NUMBER) { argv_push(a, flag); argv_push_double(a, v->n); }
}
static void push_int_flag_if(Argv* a, const char* flag, JVNode* v) {
    if (v && v->type == JV_NUMBER) { argv_push(a, flag); argv_push_int(a, (long)v->n); }
}
static void push_bool_flag_if(Argv* a, const char* flag, JVNode* v) {
    if (v && v->type == JV_BOOL && v->b) argv_push(a, flag);
}

// ── Engine execution ──────────────────────────────────────────────────────
// Spawns the engine binary with the given argv, redirects its stdout to a
// temporary file path returned in tmp_path_out. Returns the engine exit
// status, or -1 on spawn failure.

static int run_engine(const char* binary, char* const argv[], char* tmp_path_out, size_t tmp_path_sz) {
    snprintf(tmp_path_out, tmp_path_sz, "/tmp/popstats_dispatch_XXXXXX.tsv");
    int fd = mkstemps(tmp_path_out, 4);
    if (fd < 0) return -1;
    pid_t pid = fork();
    if (pid < 0) { close(fd); return -1; }
    if (pid == 0) {
        dup2(fd, STDOUT_FILENO);
        close(fd);
        execv(binary, argv);
        _exit(127);
    }
    close(fd);
    int status;
    waitpid(pid, &status, 0);
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    return -1;
}

// Same but with an output file flag (--out / --out_table / --out_regions).
// The dispatcher hands the engine a tmp path; engine writes there directly.
// Returns engine exit status.
static int run_engine_outfile(const char* binary, char* const argv[]) {
    pid_t pid = fork();
    if (pid < 0) return -1;
    if (pid == 0) {
        int dn = open("/dev/null", O_WRONLY);
        if (dn >= 0) { dup2(dn, STDOUT_FILENO); close(dn); }
        execv(binary, argv);
        _exit(127);
    }
    int status;
    waitpid(pid, &status, 0);
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    return -1;
}

// ── TSV reader → array of records → JSON emission ─────────────────────────
// Reads a TSV from `path`, returns it as parsed strings. Caller frees.

typedef struct {
    char** col_names;
    int n_cols;
    char*** rows;
    int n_rows;
    int row_cap;
    // Captured metadata key=value pairs from `# schema_version=... key=val key=val` lines.
    char** meta_keys;
    char** meta_vals;
    int    n_meta;
} TsvOut;

static char** split_tab_dup(const char* line, int* n_out) {
    int cap = 16, n = 0;
    char** out = (char**)malloc((size_t)cap * sizeof(char*));
    const char* p = line;
    while (1) {
        const char* t = strchr(p, '\t');
        size_t len = t ? (size_t)(t - p) : strlen(p);
        char* cell = (char*)malloc(len + 1);
        memcpy(cell, p, len); cell[len] = 0;
        if (n == cap) { cap *= 2; out = (char**)realloc(out, (size_t)cap * sizeof(char*)); }
        out[n++] = cell;
        if (!t) break;
        p = t + 1;
    }
    *n_out = n;
    return out;
}

static void tsv_parse_meta_line(const char* line, TsvOut* t) {
    // line starts with "# ", then "key=value key=value ..." pairs
    const char* p = line;
    while (*p == '#' || *p == ' ') p++;
    while (*p) {
        const char* eq = strchr(p, '=');
        if (!eq) break;
        const char* ke = eq;
        while (ke > p && isspace((unsigned char)ke[-1])) ke--;
        const char* vb = eq + 1;
        const char* ve = vb;
        while (*ve && !isspace((unsigned char)*ve)) ve++;
        if (eq > p && ve > vb) {
            t->meta_keys = (char**)realloc(t->meta_keys, (size_t)(t->n_meta + 1) * sizeof(char*));
            t->meta_vals = (char**)realloc(t->meta_vals, (size_t)(t->n_meta + 1) * sizeof(char*));
            t->meta_keys[t->n_meta] = strndup(p, (size_t)(ke - p));
            t->meta_vals[t->n_meta] = strndup(vb, (size_t)(ve - vb));
            t->n_meta++;
        }
        p = ve;
        while (*p && isspace((unsigned char)*p)) p++;
    }
}

static int tsv_load_file(const char* path, TsvOut* t) {
    FILE* f = fopen(path, "r");
    if (!f) return 0;
    memset(t, 0, sizeof(*t));
    char* line = NULL; size_t cap = 0; ssize_t got;
    int header_done = 0;
    while ((got = getline(&line, &cap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0) continue;
        if (line[0] == '#') { tsv_parse_meta_line(line, t); continue; }
        int nf;
        char** parts = split_tab_dup(line, &nf);
        if (!header_done) {
            t->n_cols = nf;
            t->col_names = parts;
            header_done = 1;
            continue;
        }
        if (nf != t->n_cols) {
            for (int j = 0; j < nf; j++) free(parts[j]); free(parts);
            continue;
        }
        if (t->n_rows == t->row_cap) {
            t->row_cap = t->row_cap ? t->row_cap * 2 : 256;
            t->rows = (char***)realloc(t->rows, (size_t)t->row_cap * sizeof(char**));
        }
        t->rows[t->n_rows++] = parts;
    }
    free(line); fclose(f);
    return t->n_rows;
}

static void tsv_out_free(TsvOut* t) {
    for (int j = 0; j < t->n_cols; j++) free(t->col_names[j]);
    free(t->col_names);
    for (int i = 0; i < t->n_rows; i++) {
        for (int j = 0; j < t->n_cols; j++) free(t->rows[i][j]);
        free(t->rows[i]);
    }
    free(t->rows);
    for (int i = 0; i < t->n_meta; i++) { free(t->meta_keys[i]); free(t->meta_vals[i]); }
    free(t->meta_keys); free(t->meta_vals);
    memset(t, 0, sizeof(*t));
}

static const char* meta_get(const TsvOut* t, const char* key) {
    for (int i = 0; i < t->n_meta; i++) if (!strcmp(t->meta_keys[i], key)) return t->meta_vals[i];
    return NULL;
}
static int meta_get_int(const TsvOut* t, const char* key, int fallback) {
    const char* s = meta_get(t, key);
    if (!s) return fallback;
    char* end; long v = strtol(s, &end, 10);
    return (end == s) ? fallback : (int)v;
}
static double meta_get_num(const TsvOut* t, const char* key, double fallback) {
    const char* s = meta_get(t, key);
    if (!s) return fallback;
    char* end; double v = strtod(s, &end);
    return (end == s) ? fallback : v;
}
static int meta_get_bool(const TsvOut* t, const char* key, int fallback) {
    const char* s = meta_get(t, key);
    if (!s) return fallback;
    if (!strcasecmp(s, "true") || !strcmp(s, "1")) return 1;
    if (!strcasecmp(s, "false") || !strcmp(s, "0")) return 0;
    return fallback;
}

static int tsv_col(const TsvOut* t, const char* name) {
    for (int j = 0; j < t->n_cols; j++) if (!strcmp(t->col_names[j], name)) return j;
    return -1;
}

// Emit a JSON cell value. Integers vs floats are best-effort. NA / non-numeric strings emit as strings.
static void emit_cell_typed(FILE* out, const char* cell, int as_int_hint) {
    if (!cell || !*cell || !strcasecmp(cell, "NA")) { fputs("null", out); return; }
    double d;
    if (cell_num(cell, &d)) {
        if (as_int_hint && d == floor(d) && fabs(d) < 1e15) fprintf(out, "%lld", (long long)d);
        else fprintf(out, "%.6g", d);
        return;
    }
    emit_jstr(out, cell);
}

// Heuristic: integer column if name contains 'n_', 'start', 'end', 'bp', 'length', 'count'.
static int col_is_int_hint(const char* name) {
    return (strstr(name, "n_") || strstr(name, "_n") || strstr(name, "start") ||
            strstr(name, "end") || strstr(name, "bp") || strstr(name, "length") ||
            strstr(name, "count") || !strcmp(name, "df_FIS_vsZero"));
}

// ── Per-statistic dispatchers ─────────────────────────────────────────────

static const char* engine_path_of(const char* base, const char* name, char* buf, size_t sz) {
    snprintf(buf, sz, "%s/%s", base, name);
    return buf;
}

static int dispatch_xpehh(JVNode* act, const char* engines_dir, FILE* fout) {
    JVNode* inputs = jv_get(act, "inputs");
    JVNode* scope  = jv_get(act, "scope");
    JVNode* params = jv_get(act, "params");
    JVNode* compute = jv_get(act, "compute");
    JVNode* windows = jv_get(scope, "windows");
    const char* req_id = jv_str(jv_get(act, "request_id"), NULL);

    Argv a = {0};
    char binpath[1024]; engine_path_of(engines_dir, "xpehh", binpath, sizeof(binpath));
    argv_push(&a, binpath);
    push_str_flag(&a, "--beagle",       jv_str(jv_get(inputs, "beagle"), NULL));
    push_str_flag(&a, "--sample_list",  jv_str(jv_get(inputs, "sample_list"), NULL));
    push_str_flag(&a, "--test_samples", jv_str(jv_get(inputs, "test_samples"), NULL));
    push_str_flag(&a, "--ref_samples",  jv_str(jv_get(inputs, "ref_samples"), NULL));
    push_str_flag(&a, "--chr",          jv_str(jv_get(scope, "chrom"), NULL));
    const char* wmode = jv_str(jv_get(windows, "mode"), NULL);
    if (wmode && !strcmp(wmode, "fixed_bp")) {
        long s = (long)jv_num(jv_get(windows, "size_bp"), 0);
        long step = (long)jv_num(jv_get(windows, "step_bp"), 0);
        char buf[64]; snprintf(buf, sizeof(buf), "%ld:%ld", s, step);
        argv_push(&a, "--fixed_win"); argv_push(&a, buf);
    } else if (wmode && !strcmp(wmode, "bed")) {
        push_str_flag_if(&a, "--windows", jv_get(windows, "bed_path"));
    }
    push_num_flag_if(&a, "--ehh_decay",       jv_get(params, "ehh_decay"));
    push_int_flag_if(&a, "--max_distance_bp", jv_get(params, "max_distance_bp"));
    push_int_flag_if(&a, "--min_n_test",      jv_get(params, "min_n_test"));
    push_int_flag_if(&a, "--min_n_ref",       jv_get(params, "min_n_ref"));
    push_int_flag_if(&a, "--focal_step",      jv_get(params, "focal_step"));
    push_int_flag_if(&a, "--ncores",          jv_get(compute, "ncores"));

    char tmp[256] = {0};
    snprintf(tmp, sizeof(tmp), "/tmp/popstats_dispatch_xpehh_XXXXXX");
    int fd = mkstemp(tmp);
    if (fd < 0) { argv_free(&a); return 11; }
    close(fd);
    argv_push(&a, "--out"); argv_push(&a, tmp);
    argv_push(&a, NULL); // execv terminator

    int rc = run_engine_outfile(binpath, a.argv);
    argv_free(&a);
    if (rc != 0) { unlink(tmp); return 12; }

    TsvOut t; if (!tsv_load_file(tmp, &t)) { unlink(tmp); return 13; }
    unlink(tmp);

    fputs("{", fout);
    fputs("\"schema_version\":\"xpehh_v1\",", fout);
    fputs("\"statistic\":\"xpehh\",", fout);
    if (req_id) { fputs("\"request_id\":", fout); emit_jstr(fout, req_id); fputc(',', fout); }
    fputs("\"status\":", fout);
    fputs(t.n_rows > 0 ? "\"ok\"," : "\"empty\",", fout);
    fputs("\"metadata\":{", fout);
    fputs("\"chrom\":", fout); emit_jstr(fout, jv_str(jv_get(scope, "chrom"), ""));
    fprintf(fout, ",\"n_test\":%d", meta_get_int(&t, "n_test", 0));
    fprintf(fout, ",\"n_ref\":%d", meta_get_int(&t, "n_ref", 0));
    fprintf(fout, ",\"ehh_decay\":%.6g", meta_get_num(&t, "ehh_decay", 0.05));
    fprintf(fout, ",\"max_distance_bp\":%d", meta_get_int(&t, "max_distance_bp", 1000000));
    fprintf(fout, ",\"focal_step\":%d", meta_get_int(&t, "focal_step", 1));
    fprintf(fout, ",\"n_focal_total\":%d", t.n_rows);
    int n_with = 0;
    int xpe_idx = tsv_col(&t, "xpehh_mean");
    if (xpe_idx < 0) xpe_idx = tsv_col(&t, "n_focal_with_xpehh");
    int n_with_col = tsv_col(&t, "n_focal_with_xpehh");
    if (n_with_col >= 0) {
        for (int i = 0; i < t.n_rows; i++) {
            long v; char* end;
            v = strtol(t.rows[i][n_with_col], &end, 10);
            if (end != t.rows[i][n_with_col]) n_with += (int)v;
        }
    }
    fprintf(fout, ",\"n_focal_with_xpehh\":%d", n_with);
    fprintf(fout, ",\"global_mean\":%.6g", meta_get_num(&t, "global_mean", NAN));
    fprintf(fout, ",\"global_sd\":%.6g", meta_get_num(&t, "global_sd", NAN));
    fputs(",\"approximation_note\":\"Diploid-IBS-likelihood approximation of XP-EHH from unphased BEAGLE GLs.\"", fout);
    fputs("},", fout);

    fputs("\"windows\":[", fout);
    int cols[16];
    static const char* WIN_COLS[] = {
        "chrom","window_start","window_end","n_focal","n_focal_with_xpehh",
        "iHH_test_mean","iHH_ref_mean","xpehh_mean","xpehh_max_abs",
        "norm_xpehh_mean","norm_xpehh_max_abs", NULL
    };
    int nc = 0;
    for (; WIN_COLS[nc]; nc++) cols[nc] = tsv_col(&t, WIN_COLS[nc]);
    for (int i = 0; i < t.n_rows; i++) {
        if (i) fputc(',', fout);
        fputc('{', fout);
        for (int j = 0; j < nc; j++) {
            if (j) fputc(',', fout);
            fprintf(fout, "\"%s\":", WIN_COLS[j]);
            if (cols[j] < 0) { fputs("null", fout); continue; }
            if (j == 0) emit_jstr(fout, t.rows[i][cols[j]]);
            else        emit_cell_typed(fout, t.rows[i][cols[j]], col_is_int_hint(WIN_COLS[j]));
        }
        fputc('}', fout);
    }
    fputs("]}", fout);
    fputc('\n', fout);
    tsv_out_free(&t);
    (void)xpe_idx;
    return 0;
}

static int dispatch_iHS(JVNode* act, const char* engines_dir, FILE* fout) {
    JVNode* inputs = jv_get(act, "inputs");
    JVNode* scope  = jv_get(act, "scope");
    JVNode* params = jv_get(act, "params");
    JVNode* compute = jv_get(act, "compute");
    JVNode* windows = jv_get(scope, "windows");
    const char* req_id = jv_str(jv_get(act, "request_id"), NULL);

    Argv a = {0};
    char binpath[1024]; engine_path_of(engines_dir, "iHS", binpath, sizeof(binpath));
    argv_push(&a, binpath);
    push_str_flag(&a, "--beagle",      jv_str(jv_get(inputs, "beagle"), NULL));
    push_str_flag(&a, "--sample_list", jv_str(jv_get(inputs, "sample_list"), NULL));
    push_str_flag(&a, "--cohort",      jv_str(jv_get(inputs, "cohort"), NULL));
    push_str_flag(&a, "--chr",         jv_str(jv_get(scope, "chrom"), NULL));
    const char* wmode = jv_str(jv_get(windows, "mode"), NULL);
    if (wmode && !strcmp(wmode, "fixed_bp")) {
        long s = (long)jv_num(jv_get(windows, "size_bp"), 0);
        long step = (long)jv_num(jv_get(windows, "step_bp"), 0);
        char buf[64]; snprintf(buf, sizeof(buf), "%ld:%ld", s, step);
        argv_push(&a, "--fixed_win"); argv_push(&a, buf);
    } else if (wmode && !strcmp(wmode, "bed")) {
        push_str_flag_if(&a, "--windows", jv_get(windows, "bed_path"));
    }
    push_num_flag_if(&a, "--ehh_decay",       jv_get(params, "ehh_decay"));
    push_int_flag_if(&a, "--max_distance_bp", jv_get(params, "max_distance_bp"));
    push_num_flag_if(&a, "--gl_threshold",    jv_get(params, "gl_threshold"));
    push_int_flag_if(&a, "--min_n_anc",       jv_get(params, "min_n_anc"));
    push_int_flag_if(&a, "--min_n_der",       jv_get(params, "min_n_der"));
    push_num_flag_if(&a, "--min_freq",        jv_get(params, "min_freq"));
    push_int_flag_if(&a, "--n_freq_bins",     jv_get(params, "n_freq_bins"));
    push_int_flag_if(&a, "--focal_step",      jv_get(params, "focal_step"));
    push_int_flag_if(&a, "--ncores",          jv_get(compute, "ncores"));

    char tmp[256] = "/tmp/popstats_dispatch_iHS_XXXXXX";
    int fd = mkstemp(tmp);
    if (fd < 0) { argv_free(&a); return 11; }
    close(fd);
    argv_push(&a, "--out"); argv_push(&a, tmp);
    argv_push(&a, NULL);
    int rc = run_engine_outfile(binpath, a.argv);
    argv_free(&a);
    if (rc != 0) { unlink(tmp); return 12; }

    TsvOut t; if (!tsv_load_file(tmp, &t)) { unlink(tmp); return 13; }
    unlink(tmp);

    fputs("{", fout);
    fputs("\"schema_version\":\"iHS_v1\",", fout);
    fputs("\"statistic\":\"iHS\",", fout);
    if (req_id) { fputs("\"request_id\":", fout); emit_jstr(fout, req_id); fputc(',', fout); }
    fputs("\"status\":", fout);
    fputs(t.n_rows > 0 ? "\"ok\"," : "\"empty\",", fout);
    fputs("\"metadata\":{", fout);
    fputs("\"chrom\":", fout); emit_jstr(fout, jv_str(jv_get(scope, "chrom"), ""));
    fprintf(fout, ",\"n_cohort\":%d", meta_get_int(&t, "n_cohort", 0));
    fprintf(fout, ",\"ehh_decay\":%.6g", meta_get_num(&t, "ehh_decay", 0.05));
    fprintf(fout, ",\"max_distance_bp\":%d", meta_get_int(&t, "max_distance_bp", 1000000));
    fprintf(fout, ",\"gl_threshold\":%.6g", meta_get_num(&t, "gl_threshold", 0.6));
    fprintf(fout, ",\"min_freq\":%.6g", meta_get_num(&t, "min_freq", 0.05));
    fprintf(fout, ",\"n_freq_bins\":%d", meta_get_int(&t, "n_freq_bins", 20));
    fprintf(fout, ",\"focal_step\":%d", meta_get_int(&t, "focal_step", 1));
    fprintf(fout, ",\"n_focal_total\":%d", t.n_rows);
    fputs(",\"ancestral_polarization\":\"allele_0_assumed_ancestral\"", fout);
    fputs(",\"approximation_note\":\"Diploid-IBS-likelihood approximation. Voight-strict iHS needs phasing.\"", fout);
    fputs("},\"windows\":[", fout);

    static const char* WC[] = {
        "chrom","window_start","window_end","n_focal","n_focal_with_iHS",
        "iHH_anc_mean","iHH_der_mean","iHS_mean","iHS_max_abs",
        "norm_iHS_mean","norm_iHS_max_abs","n_extreme_iHS_abs_gt_2", NULL
    };
    int cols[16]; int nc = 0;
    for (; WC[nc]; nc++) cols[nc] = tsv_col(&t, WC[nc]);
    for (int i = 0; i < t.n_rows; i++) {
        if (i) fputc(',', fout);
        fputc('{', fout);
        for (int j = 0; j < nc; j++) {
            if (j) fputc(',', fout);
            fprintf(fout, "\"%s\":", WC[j]);
            if (cols[j] < 0) { fputs("null", fout); continue; }
            if (j == 0) emit_jstr(fout, t.rows[i][cols[j]]);
            else        emit_cell_typed(fout, t.rows[i][cols[j]], col_is_int_hint(WC[j]));
        }
        fputc('}', fout);
    }
    fputs("]}", fout);
    fputc('\n', fout);
    tsv_out_free(&t);
    return 0;
}

static int dispatch_outlier_scan(JVNode* act, const char* engines_dir, FILE* fout) {
    JVNode* inputs = jv_get(act, "inputs");
    JVNode* stat = jv_get(act, "stat");
    JVNode* thr = jv_get(act, "thresholds");
    JVNode* meta = jv_get(act, "metadata");
    const char* req_id = jv_str(jv_get(act, "request_id"), NULL);

    Argv a = {0};
    char binpath[1024]; engine_path_of(engines_dir, "outlier_scan", binpath, sizeof(binpath));
    argv_push(&a, binpath);
    push_str_flag(&a, "--input",             jv_str(jv_get(inputs, "windows_tsv"), NULL));
    push_str_flag(&a, "--stronger_windows",  jv_str(jv_get(inputs, "stronger_windows_tsv"), NULL));
    push_str_flag(&a, "--genes_tsv",         jv_str(jv_get(inputs, "genes_tsv"), NULL));
    push_str_flag(&a, "--inversion_bed",     jv_str(jv_get(inputs, "inversion_bed"), NULL));
    push_str_flag(&a, "--breakpoint_bed",    jv_str(jv_get(inputs, "breakpoint_bed"), NULL));
    push_str_flag(&a, "--te_bed",            jv_str(jv_get(inputs, "te_bed"), NULL));
    push_str_flag(&a, "--stat_name",         jv_str(jv_get(stat, "name"), "stat"));
    push_str_flag(&a, "--direction",         jv_str(jv_get(stat, "direction"), "high"));
    push_str_flag(&a, "--stat_col",          jv_str(jv_get(stat, "stat_col"), NULL));
    push_str_flag(&a, "--chrom_col",         jv_str(jv_get(stat, "chrom_col"), NULL));
    push_str_flag(&a, "--start_col",         jv_str(jv_get(stat, "start_col"), NULL));
    push_str_flag(&a, "--end_col",           jv_str(jv_get(stat, "end_col"), NULL));
    push_str_flag(&a, "--nvar_col",          jv_str(jv_get(stat, "nvar_col"), NULL));
    push_num_flag_if(&a, "--threshold",          jv_get(thr, "threshold"));
    push_num_flag_if(&a, "--top_pct",            jv_get(thr, "top_pct"));
    push_int_flag_if(&a, "--merge_gap_bp",       jv_get(thr, "merge_gap_bp"));
    push_num_flag_if(&a, "--stronger_threshold", jv_get(thr, "stronger_threshold"));
    push_str_flag(&a, "--comparison_id", jv_str(jv_get(meta, "comparison_id"), NULL));
    push_str_flag(&a, "--group_a",       jv_str(jv_get(meta, "group_a"), NULL));
    push_str_flag(&a, "--group_b",       jv_str(jv_get(meta, "group_b"), NULL));
    push_str_flag(&a, "--window_mode",   jv_str(jv_get(meta, "window_mode"), NULL));
    push_int_flag_if(&a, "--genome_bp",   jv_get(meta, "genome_bp"));

    char regions_tmp[256] = "/tmp/popstats_dispatch_outlier_regions_XXXXXX";
    char summary_tmp[256] = "/tmp/popstats_dispatch_outlier_summary_XXXXXX";
    int fdr = mkstemp(regions_tmp); int fds = mkstemp(summary_tmp);
    if (fdr < 0 || fds < 0) { argv_free(&a); return 11; }
    close(fdr); close(fds);
    argv_push(&a, "--out_regions"); argv_push(&a, regions_tmp);
    argv_push(&a, "--out_summary"); argv_push(&a, summary_tmp);
    argv_push(&a, NULL);

    int rc = run_engine_outfile(binpath, a.argv);
    argv_free(&a);
    if (rc != 0) { unlink(regions_tmp); unlink(summary_tmp); return 12; }

    TsvOut t; if (!tsv_load_file(regions_tmp, &t)) { unlink(regions_tmp); unlink(summary_tmp); return 13; }
    unlink(regions_tmp);

    // Pipe the engine's JSON summary directly into our metadata.
    FILE* sf = fopen(summary_tmp, "r");
    char* summary_json = NULL; size_t summary_len = 0;
    if (sf) {
        fseek(sf, 0, SEEK_END); summary_len = (size_t)ftell(sf); fseek(sf, 0, SEEK_SET);
        summary_json = (char*)malloc(summary_len + 1);
        if (fread(summary_json, 1, summary_len, sf) != summary_len) { /* best effort */ }
        summary_json[summary_len] = 0;
        fclose(sf);
    }
    unlink(summary_tmp);

    fputs("{", fout);
    fputs("\"schema_version\":\"outlier_scan_v1\",", fout);
    fputs("\"statistic\":\"outlier_scan\",", fout);
    if (req_id) { fputs("\"request_id\":", fout); emit_jstr(fout, req_id); fputc(',', fout); }
    fputs("\"status\":", fout); fputs(t.n_rows > 0 ? "\"ok\"," : "\"empty\",", fout);
    if (summary_json) {
        // Embed the engine's summary as metadata (already JSON-object shaped)
        const char* p = summary_json;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '{') {
            fputs("\"metadata\":", fout);
            fwrite(p, 1, strlen(p), fout);
            // strip trailing newline if any
            fputc(',', fout);
        }
        free(summary_json);
    }

    static const char* RC[] = {
        "region_id","chrom","start","end","length_bp","n_windows","n_variants_total",
        "stat_mean","stat_median","stat_max","stat_min",
        "threshold_used","stronger_threshold_used","contains_strong_window","class", NULL
    };
    // outlier_scan emits <stat_name>_mean etc. — column names start with the stat name.
    const char* stat_name_used = jv_str(jv_get(stat, "name"), "stat");
    char stat_mean[64], stat_median[64], stat_max[64], stat_min[64];
    snprintf(stat_mean,   sizeof(stat_mean),   "%s_mean",   stat_name_used);
    snprintf(stat_median, sizeof(stat_median), "%s_median", stat_name_used);
    snprintf(stat_max,    sizeof(stat_max),    "%s_max",    stat_name_used);
    snprintf(stat_min,    sizeof(stat_min),    "%s_min",    stat_name_used);
    int idx_id        = tsv_col(&t, "region_id");
    int idx_chrom     = tsv_col(&t, "chrom");
    int idx_start     = tsv_col(&t, "start");
    int idx_end       = tsv_col(&t, "end");
    int idx_len       = tsv_col(&t, "length_bp");
    int idx_nwin      = tsv_col(&t, "n_windows");
    int idx_nvar      = tsv_col(&t, "n_variants_total");
    int idx_mean      = tsv_col(&t, stat_mean);
    int idx_median    = tsv_col(&t, stat_median);
    int idx_max       = tsv_col(&t, stat_max);
    int idx_min       = tsv_col(&t, stat_min);
    int idx_thr       = tsv_col(&t, "threshold_used");
    int idx_strthr    = tsv_col(&t, "stronger_threshold_used");
    int idx_cont      = tsv_col(&t, "contains_strong_window");
    int idx_class     = tsv_col(&t, "class");
    int idx_gene      = tsv_col(&t, "nearest_gene");
    int idx_oi        = tsv_col(&t, "overlaps_inversion_candidate");
    int idx_obp       = tsv_col(&t, "overlaps_breakpoint");
    int idx_ote       = tsv_col(&t, "overlaps_TE_hotspot");

    fputs("\"regions\":[", fout);
    for (int i = 0; i < t.n_rows; i++) {
        if (i) fputc(',', fout);
        fputc('{', fout);
        #define EMIT_STR(name, idx) do { fprintf(fout, "\"%s\":", name); if (idx<0) fputs("null", fout); else emit_jstr(fout, t.rows[i][idx]); } while (0)
        #define EMIT_TYPED(name, idx, as_int) do { fprintf(fout, "\"%s\":", name); if (idx<0) fputs("null", fout); else emit_cell_typed(fout, t.rows[i][idx], as_int); } while (0)
        #define EMIT_BOOL(name, idx) do { fprintf(fout, "\"%s\":", name); if (idx<0) fputs("null", fout); \
            else if (!strcasecmp(t.rows[i][idx],"TRUE")||!strcmp(t.rows[i][idx],"1")) fputs("true", fout); \
            else if (!strcasecmp(t.rows[i][idx],"FALSE")||!strcmp(t.rows[i][idx],"0")) fputs("false", fout); \
            else fputs("null", fout); } while (0)
        EMIT_STR("region_id", idx_id); fputc(',', fout);
        EMIT_STR("chrom", idx_chrom); fputc(',', fout);
        EMIT_TYPED("start", idx_start, 1); fputc(',', fout);
        EMIT_TYPED("end", idx_end, 1); fputc(',', fout);
        EMIT_TYPED("length_bp", idx_len, 1); fputc(',', fout);
        EMIT_TYPED("n_windows", idx_nwin, 1); fputc(',', fout);
        EMIT_TYPED("n_variants_total", idx_nvar, 1); fputc(',', fout);
        EMIT_TYPED("stat_mean", idx_mean, 0); fputc(',', fout);
        EMIT_TYPED("stat_median", idx_median, 0); fputc(',', fout);
        EMIT_TYPED("stat_max", idx_max, 0); fputc(',', fout);
        EMIT_TYPED("stat_min", idx_min, 0); fputc(',', fout);
        EMIT_TYPED("threshold_used", idx_thr, 0); fputc(',', fout);
        EMIT_TYPED("stronger_threshold_used", idx_strthr, 0); fputc(',', fout);
        EMIT_BOOL("contains_strong_window", idx_cont); fputc(',', fout);
        EMIT_STR("class", idx_class);
        if (idx_gene >= 0) { fputc(',', fout); EMIT_STR("nearest_gene", idx_gene); }
        if (idx_oi >= 0)   { fputc(',', fout); EMIT_BOOL("overlaps_inversion_candidate", idx_oi); }
        if (idx_obp >= 0)  { fputc(',', fout); EMIT_BOOL("overlaps_breakpoint", idx_obp); }
        if (idx_ote >= 0)  { fputc(',', fout); EMIT_BOOL("overlaps_TE_hotspot", idx_ote); }
        #undef EMIT_STR
        #undef EMIT_TYPED
        #undef EMIT_BOOL
        fputc('}', fout);
    }
    fputs("]}", fout);
    fputc('\n', fout);
    tsv_out_free(&t);
    (void)RC;
    return 0;
}

static int dispatch_candidate_vs_flanks(JVNode* act, const char* engines_dir, FILE* fout) {
    JVNode* inputs = jv_get(act, "inputs");
    JVNode* stats  = jv_get(act, "stats");
    JVNode* table_cols = jv_get(act, "table_cols");
    JVNode* params = jv_get(act, "params");
    JVNode* compute = jv_get(act, "compute");
    const char* req_id = jv_str(jv_get(act, "request_id"), NULL);

    // Build comma-joined stat_cols and stat_names
    char stat_cols_csv[2048] = "", stat_names_csv[2048] = "";
    int nstats = jv_arr_len(stats);
    int have_names = 0;
    for (int i = 0; i < nstats; i++) {
        JVNode* s = jv_arr_at(stats, i);
        const char* col = jv_str(jv_get(s, "col"), NULL);
        const char* nm  = jv_str(jv_get(s, "name"), col);
        if (i) { strncat(stat_cols_csv, ",", sizeof(stat_cols_csv) - strlen(stat_cols_csv) - 1); strncat(stat_names_csv, ",", sizeof(stat_names_csv) - strlen(stat_names_csv) - 1); }
        if (col) strncat(stat_cols_csv, col, sizeof(stat_cols_csv) - strlen(stat_cols_csv) - 1);
        if (nm)  { strncat(stat_names_csv, nm, sizeof(stat_names_csv) - strlen(stat_names_csv) - 1); have_names = 1; }
    }

    Argv a = {0};
    char binpath[1024]; engine_path_of(engines_dir, "candidate_vs_flanks", binpath, sizeof(binpath));
    argv_push(&a, binpath);
    push_str_flag(&a, "--candidates", jv_str(jv_get(inputs, "candidates_bed"), NULL));
    push_str_flag(&a, "--windows_tsv", jv_str(jv_get(inputs, "windows_tsv"), NULL));
    argv_push(&a, "--stat_cols"); argv_push(&a, stat_cols_csv);
    if (have_names) { argv_push(&a, "--stat_names"); argv_push(&a, stat_names_csv); }
    push_str_flag(&a, "--chrom_col", jv_str(jv_get(table_cols, "chrom_col"), NULL));
    push_str_flag(&a, "--start_col", jv_str(jv_get(table_cols, "start_col"), NULL));
    push_str_flag(&a, "--end_col",   jv_str(jv_get(table_cols, "end_col"), NULL));
    push_int_flag_if(&a, "--flank_bp",     jv_get(params, "flank_bp"));
    push_bool_flag_if(&a, "--exclude_flank_overlaps", jv_get(params, "exclude_flank_overlaps"));
    push_int_flag_if(&a, "--permutations", jv_get(params, "permutations"));
    push_int_flag_if(&a, "--seed",         jv_get(params, "seed"));
    push_int_flag_if(&a, "--ncores",       jv_get(compute, "ncores"));

    char table_tmp[256] = "/tmp/popstats_dispatch_cvf_table_XXXXXX";
    char sum_tmp[256] = "/tmp/popstats_dispatch_cvf_sum_XXXXXX";
    int fdt = mkstemp(table_tmp); int fds = mkstemp(sum_tmp);
    if (fdt < 0 || fds < 0) { argv_free(&a); return 11; }
    close(fdt); close(fds);
    argv_push(&a, "--out_table"); argv_push(&a, table_tmp);
    argv_push(&a, "--out_summary"); argv_push(&a, sum_tmp);
    argv_push(&a, NULL);

    int rc = run_engine_outfile(binpath, a.argv);
    argv_free(&a);
    if (rc != 0) { unlink(table_tmp); unlink(sum_tmp); return 12; }

    TsvOut t; if (!tsv_load_file(table_tmp, &t)) { unlink(table_tmp); unlink(sum_tmp); return 13; }
    unlink(table_tmp);

    FILE* sf = fopen(sum_tmp, "r");
    char* summary_json = NULL;
    if (sf) {
        fseek(sf, 0, SEEK_END); size_t slen = (size_t)ftell(sf); fseek(sf, 0, SEEK_SET);
        summary_json = (char*)malloc(slen + 1);
        if (fread(summary_json, 1, slen, sf) != slen) { /* best effort */ }
        summary_json[slen] = 0;
        fclose(sf);
    }
    unlink(sum_tmp);

    fputs("{", fout);
    fputs("\"schema_version\":\"candidate_vs_flanks_v1\",", fout);
    fputs("\"statistic\":\"candidate_vs_flanks\",", fout);
    if (req_id) { fputs("\"request_id\":", fout); emit_jstr(fout, req_id); fputc(',', fout); }
    fputs("\"status\":", fout); fputs(t.n_rows > 0 ? "\"ok\"," : "\"empty\",", fout);
    if (summary_json) {
        const char* p = summary_json;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '{') {
            fputs("\"metadata\":", fout);
            fputs(p, fout);
            fputc(',', fout);
        }
        free(summary_json);
    }

    fputs("\"tests\":[", fout);
    static const char* TC[] = {
        "candidate_id","chrom","cand_start","cand_end","length_bp","flank_bp","stat_name",
        "n_inside","n_flank","mean_inside","mean_flank","median_inside","median_flank",
        "effect_diff","inside_minus_flank_direction","W","z","wilcoxon_p","perm_p", NULL
    };
    int cols[32]; int nc = 0;
    for (; TC[nc]; nc++) cols[nc] = tsv_col(&t, TC[nc]);
    for (int i = 0; i < t.n_rows; i++) {
        if (i) fputc(',', fout);
        fputc('{', fout);
        for (int j = 0; j < nc; j++) {
            if (j) fputc(',', fout);
            fprintf(fout, "\"%s\":", TC[j]);
            if (cols[j] < 0) { fputs("null", fout); continue; }
            // String-typed columns: chrom, candidate_id, stat_name, direction
            if (j == 0 || j == 1 || j == 6 || !strcmp(TC[j], "inside_minus_flank_direction"))
                emit_jstr(fout, t.rows[i][cols[j]]);
            else
                emit_cell_typed(fout, t.rows[i][cols[j]], col_is_int_hint(TC[j]));
        }
        fputc('}', fout);
    }
    fputs("]}", fout);
    fputc('\n', fout);
    tsv_out_free(&t);
    return 0;
}

// ── Main ──────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "popstats_dispatch — JSON activator → engine → JSON extractor.\n"
        "\n"
        "  --in <f>          Read JSON activator from file (default: stdin).\n"
        "  --out <f>         Write JSON extractor to file (default: stdout).\n"
        "  --engines_dir <d> Path to the engines/ directory (default: dirname(argv[0])).\n"
        "  -h, --help        This help.\n"
        "\n"
        "Supported statistics (per popstats_server.activator.schema.json):\n"
        "  xpehh, iHS, outlier_scan, candidate_vs_flanks\n"
        "\n"
        "Errors emit a JSON envelope with status=\"error\" and an `error` field.\n");
}

static void emit_error(FILE* out, const char* req_id, const char* msg) {
    fputs("{\"schema_version\":\"" DISPATCH_VERSION "\",\"status\":\"error\",", out);
    if (req_id) { fputs("\"request_id\":", out); emit_jstr(out, req_id); fputc(',', out); }
    fputs("\"error\":", out); emit_jstr(out, msg);
    fputs("}\n", out);
}

int main(int argc, char** argv) {
    const char* in_path = NULL;
    const char* out_path = NULL;
    char engines_dir_buf[1024] = "";
    const char* engines_dir = NULL;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--in")  && i+1<argc) in_path = argv[++i];
        else if (!strcmp(argv[i], "--out") && i+1<argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "--engines_dir") && i+1<argc) engines_dir = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[dispatch] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }

    if (!engines_dir) {
        strncpy(engines_dir_buf, argv[0], sizeof(engines_dir_buf) - 1);
        char* slash = strrchr(engines_dir_buf, '/');
        if (slash) *slash = 0;
        else        strcpy(engines_dir_buf, ".");
        engines_dir = engines_dir_buf;
    }

    FILE* in = in_path ? fopen(in_path, "r") : stdin;
    if (!in) { fprintf(stderr, "[dispatch] Cannot open --in %s\n", in_path); return 2; }
    FILE* out = out_path ? fopen(out_path, "w") : stdout;
    if (!out) { fprintf(stderr, "[dispatch] Cannot open --out %s\n", out_path); return 3; }

    // Slurp the entire activator JSON.
    size_t cap = 65536, n = 0;
    char* buf = (char*)malloc(cap);
    int c;
    while ((c = fgetc(in)) != EOF) {
        if (n + 1 >= cap) { cap *= 2; buf = (char*)realloc(buf, cap); }
        buf[n++] = (char)c;
    }
    buf[n] = 0;
    if (in_path) fclose(in);

    char err_buf[256];
    JVNode* act = jv_parse(buf, err_buf, sizeof(err_buf));
    free(buf);
    if (!act) {
        emit_error(out, NULL, err_buf[0] ? err_buf : "JSON parse failed");
        if (out_path) fclose(out);
        return 4;
    }
    if (act->type != JV_OBJECT) {
        emit_error(out, NULL, "activator must be a JSON object");
        jv_free(act); if (out_path) fclose(out); return 5;
    }

    const char* stat = jv_str(jv_get(act, "statistic"), NULL);
    const char* req_id = jv_str(jv_get(act, "request_id"), NULL);
    if (!stat) {
        emit_error(out, req_id, "activator is missing 'statistic' field");
        jv_free(act); if (out_path) fclose(out); return 6;
    }
    fprintf(stderr, "[dispatch] statistic=%s engines_dir=%s\n", stat, engines_dir);

    int rc = 7;
    if      (!strcmp(stat, "xpehh"))               rc = dispatch_xpehh(act, engines_dir, out);
    else if (!strcmp(stat, "iHS"))                 rc = dispatch_iHS(act, engines_dir, out);
    else if (!strcmp(stat, "outlier_scan"))        rc = dispatch_outlier_scan(act, engines_dir, out);
    else if (!strcmp(stat, "candidate_vs_flanks")) rc = dispatch_candidate_vs_flanks(act, engines_dir, out);
    else {
        char msg[256];
        snprintf(msg, sizeof(msg), "unsupported statistic '%s' (try: xpehh, iHS, outlier_scan, candidate_vs_flanks)", stat);
        emit_error(out, req_id, msg);
        rc = 8;
    }
    if (rc != 0 && rc < 100) {
        fprintf(stderr, "[dispatch] handler returned %d for %s\n", rc, stat);
    }

    jv_free(act);
    if (out_path) fclose(out);
    return rc;
}
