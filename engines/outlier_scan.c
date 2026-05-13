// =============================================================================
// outlier_scan.c — genome-wide statistic-agnostic outlier + HDR scan.
//
// Generalization of fst_outlier_scan: same threshold → flag → gap-merge → HDR
// pipeline, but works on any numeric column of a per-window TSV
// (FST, π, dXY, dA, Tajima's D, fdM, ...). Direction-aware:
//
//   --direction high   flag rows with value >= --threshold   (FST, dXY, fdM>0, ...)
//   --direction low    flag rows with value <= --threshold   (π troughs, low recombination)
//   --direction abs    flag rows with |value| >= --threshold (fdM both ways, Tajima's D)
//
// For --top_pct, the direction picks the appropriate tail percentile.
//
// Output:
//   --out_windows   <f>   flagged input windows (passing the direction-aware filter)
//   --out_regions   <f>   merged regions, class = OUTLIER or HDR
//   --out_hdr       <f>   HDR-class regions only
//   --out_summary   <f>   JSON summary
//
// Output columns and JSON keys use a configurable --stat_name (default "stat"),
// so e.g. --stat_name fst gives fst_mean/fst_max/..., --stat_name pi gives
// pi_mean/pi_max/..., etc.
//
// Caveat phrasing: emits "empirical <stat> outliers" / "highly extreme
// regions" — NOT "selected regions". No neutral calibration is done here.
//
// Compile: gcc -O3 -march=native -Wall -o outlier_scan outlier_scan.c -lm
// =============================================================================

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>

#define SCHEMA_VERSION "outlier_scan_v1"

// ── Direction enum ─────────────────────────────────────────────────────────
typedef enum { DIR_HIGH = 0, DIR_LOW = 1, DIR_ABS = 2 } Direction;

static int parse_direction(const char* s, Direction* d) {
    if      (!strcasecmp(s, "high") || !strcasecmp(s, "upper")) { *d = DIR_HIGH; return 1; }
    else if (!strcasecmp(s, "low")  || !strcasecmp(s, "lower")) { *d = DIR_LOW;  return 1; }
    else if (!strcasecmp(s, "abs")  || !strcasecmp(s, "two_sided") || !strcasecmp(s, "two-sided")) { *d = DIR_ABS; return 1; }
    return 0;
}

static const char* direction_str(Direction d) {
    return d == DIR_HIGH ? "high" : d == DIR_LOW ? "low" : "abs";
}

// Returns 1 if value passes the outlier-flag test against threshold, given direction.
static inline int passes_outlier(double v, double threshold, Direction d) {
    if (!isfinite(v)) return 0;
    if (d == DIR_HIGH) return v >= threshold;
    if (d == DIR_LOW)  return v <= threshold;
    return fabs(v) >= threshold;
}

// ── Generic TSV reader ─────────────────────────────────────────────────────

typedef struct {
    char** col_names;
    int    n_cols;
    char*** cells;
    int    n_rows;
    int    cap;
} Tsv;

static char** split_tab(const char* line, int* n_out) {
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

static int tsv_load(const char* path, Tsv* t) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[outlier_scan] Cannot open %s\n", path); return 0; }
    char* line = NULL; size_t cap = 0; ssize_t got;
    memset(t, 0, sizeof(*t));
    int header_done = 0;
    while ((got = getline(&line, &cap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        int nf;
        char** parts = split_tab(line, &nf);
        if (!header_done) {
            t->n_cols = nf;
            t->col_names = (char**)malloc((size_t)nf * sizeof(char*));
            for (int j = 0; j < nf; j++) t->col_names[j] = parts[j];
            free(parts);
            header_done = 1;
            continue;
        }
        if (nf != t->n_cols) {
            for (int j = 0; j < nf; j++) free(parts[j]);
            free(parts);
            continue;
        }
        if (t->n_rows == t->cap) {
            t->cap = t->cap ? t->cap * 2 : 4096;
            t->cells = (char***)realloc(t->cells, (size_t)t->cap * sizeof(char**));
        }
        t->cells[t->n_rows++] = parts;
    }
    free(line); fclose(f);
    return t->n_rows;
}

static void tsv_free(Tsv* t) {
    for (int j = 0; j < t->n_cols; j++) free(t->col_names[j]);
    free(t->col_names);
    for (int i = 0; i < t->n_rows; i++) {
        for (int j = 0; j < t->n_cols; j++) free(t->cells[i][j]);
        free(t->cells[i]);
    }
    free(t->cells);
    memset(t, 0, sizeof(*t));
}

static int tsv_find_col(const Tsv* t, const char* const* aliases) {
    for (int a = 0; aliases[a]; a++)
        for (int j = 0; j < t->n_cols; j++)
            if (!strcasecmp(t->col_names[j], aliases[a])) return j;
    return -1;
}
static int tsv_col_idx(const Tsv* t, const char* name) {
    for (int j = 0; j < t->n_cols; j++)
        if (!strcasecmp(t->col_names[j], name)) return j;
    return -1;
}

static int parse_d(const char* s, double* o) {
    if (!s || !*s || !strcasecmp(s, "NA") || !strcasecmp(s, "nan")) return 0;
    if (!strcasecmp(s, "inf") || !strcasecmp(s, "-inf")) return 0;
    char* end; double v = strtod(s, &end);
    if (end == s || !isfinite(v)) return 0;
    *o = v; return 1;
}
static int parse_i(const char* s, long* o) {
    if (!s || !*s) return 0;
    char* end; long v = strtol(s, &end, 10);
    if (end == s) return 0;
    *o = v; return 1;
}

// ── Window struct (parsed from input) ──────────────────────────────────────

typedef struct {
    char* chrom;
    long  start, end;
    double value;
    long   n_variants;
    int    flagged;
} Window;

static int win_cmp(const void* a, const void* b) {
    const Window* wa = (const Window*)a; const Window* wb = (const Window*)b;
    int c = strcmp(wa->chrom, wb->chrom);
    if (c) return c;
    return (wa->start < wb->start) ? -1 : (wa->start > wb->start ? 1 : 0);
}

static int dbl_cmp_asc(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y ? 1 : 0);
}

static double percentile_value(const double* arr, int n, double q) {
    if (n == 0) return NAN;
    double* tmp = (double*)malloc((size_t)n * sizeof(double));
    memcpy(tmp, arr, (size_t)n * sizeof(double));
    qsort(tmp, (size_t)n, sizeof(double), dbl_cmp_asc);
    int idx = (int)floor(q * (n - 1));
    if (idx < 0) idx = 0;
    if (idx >= n) idx = n - 1;
    double v = tmp[idx];
    free(tmp);
    return v;
}

// ── BED loader + overlap check ─────────────────────────────────────────────

typedef struct { char* chrom; long start, end; char* name; } BedRec;

static int bed_load(const char* path, BedRec** out) {
    *out = NULL;
    if (!path) return 0;
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[outlier_scan] Cannot open BED %s\n", path); return 0; }
    int cap = 1024, n = 0;
    BedRec* arr = (BedRec*)malloc((size_t)cap * sizeof(BedRec));
    char* line = NULL; size_t lcap = 0; ssize_t got;
    while ((got = getline(&line, &lcap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        char chrom[128] = "", name[128] = "";
        long s, e;
        int nf = sscanf(line, "%127s %ld %ld %127s", chrom, &s, &e, name);
        if (nf < 3) continue;
        if (n == cap) { cap *= 2; arr = (BedRec*)realloc(arr, (size_t)cap * sizeof(BedRec)); }
        arr[n].chrom = strdup(chrom); arr[n].start = s; arr[n].end = e;
        arr[n].name = (nf >= 4) ? strdup(name) : strdup("");
        n++;
    }
    free(line); fclose(f);
    qsort(arr, (size_t)n, sizeof(BedRec), (int(*)(const void*, const void*))win_cmp);
    *out = arr;
    return n;
}

static void bed_free(BedRec* arr, int n) {
    for (int i = 0; i < n; i++) { free(arr[i].chrom); free(arr[i].name); }
    free(arr);
}

static int bed_overlaps(const BedRec* bed, int n_bed, const char* chrom, long s, long e) {
    for (int i = 0; i < n_bed; i++) {
        if (strcmp(bed[i].chrom, chrom) != 0) continue;
        if (bed[i].end < s) continue;
        if (bed[i].start > e) continue;
        return 1;
    }
    return 0;
}

// Nearest gene from a TSV.
typedef struct { char* chrom; long start, end; char* name; } GeneRec;

static int genes_load(const char* path, GeneRec** out) {
    *out = NULL;
    if (!path) return 0;
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[outlier_scan] Cannot open genes TSV %s\n", path); return 0; }
    int cap = 1024, n = 0;
    GeneRec* arr = (GeneRec*)malloc((size_t)cap * sizeof(GeneRec));
    char* line = NULL; size_t lcap = 0; ssize_t got;
    int header_seen = 0;
    while ((got = getline(&line, &lcap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        char chrom[128] = "", name[128] = "";
        long s, e;
        int nf = sscanf(line, "%127s %ld %ld %127s", chrom, &s, &e, name);
        if (nf < 4) continue;
        if (!header_seen && (!strcasecmp(chrom, "chrom") || !strcasecmp(chrom, "chromosome"))) {
            header_seen = 1; continue;
        }
        header_seen = 1;
        if (n == cap) { cap *= 2; arr = (GeneRec*)realloc(arr, (size_t)cap * sizeof(GeneRec)); }
        arr[n].chrom = strdup(chrom);
        arr[n].start = s; arr[n].end = e;
        arr[n].name = strdup(name);
        n++;
    }
    free(line); fclose(f);
    *out = arr;
    return n;
}

static void genes_free(GeneRec* arr, int n) {
    for (int i = 0; i < n; i++) { free(arr[i].chrom); free(arr[i].name); }
    free(arr);
}

static const char* nearest_gene(const GeneRec* g, int n_g, const char* chrom, long mid) {
    const char* best = NULL;
    long best_dist = -1;
    for (int i = 0; i < n_g; i++) {
        if (strcmp(g[i].chrom, chrom) != 0) continue;
        long d;
        if (mid < g[i].start)      d = g[i].start - mid;
        else if (mid > g[i].end)   d = mid - g[i].end;
        else                       d = 0;
        if (best_dist < 0 || d < best_dist) { best_dist = d; best = g[i].name; }
    }
    return best ? best : "NA";
}

// ── Region (merged outlier run) ────────────────────────────────────────────

typedef struct {
    char chrom[128];
    long start, end;
    int  n_windows;
    long n_variants_total;
    double val_sum, val_min, val_max;
    double* values;
    int  values_cap;
    int  contains_strong_window;
} Region;

static void region_push(Region* r, double v) {
    if (r->n_windows == r->values_cap) {
        r->values_cap = r->values_cap ? r->values_cap * 2 : 16;
        r->values = (double*)realloc(r->values, (size_t)r->values_cap * sizeof(double));
    }
    r->values[r->n_windows++] = v;
    r->val_sum += v;
    if (v < r->val_min) r->val_min = v;
    if (v > r->val_max) r->val_max = v;
}

typedef struct {
    Region* regs;
    int nr;
    const char* comp_id;
    const char* stat_name;
    Direction direction;
    double threshold_used;
    double stronger_threshold;
    long merge_gap_bp;
    GeneRec* genes; int n_genes;
    BedRec* inv_bed; int n_inv;
    BedRec* bp_bed; int n_bp;
    BedRec* teh_bed; int n_te;
} RegionsEmitCtx;

static void emit_regions(FILE* f, const RegionsEmitCtx* c, int hdr_only) {
    fprintf(f, "# schema_version=%s kind=%s comparison_id=%s stat=%s direction=%s "
               "threshold=%.6g stronger_threshold=%.6g merge_gap_bp=%ld\n",
            SCHEMA_VERSION, hdr_only ? "HDR_regions" : "outlier_regions",
            c->comp_id, c->stat_name, direction_str(c->direction),
            c->threshold_used, c->stronger_threshold, c->merge_gap_bp);
    fprintf(f, "region_id\tchrom\tstart\tend\tlength_bp\tn_windows\tn_variants_total\t"
               "%s_mean\t%s_median\t%s_max\t%s_min\t"
               "threshold_used\tstronger_threshold_used\tcontains_strong_window\tclass",
            c->stat_name, c->stat_name, c->stat_name, c->stat_name);
    if (c->n_genes > 0)  fprintf(f, "\tnearest_gene");
    if (c->n_inv > 0)    fprintf(f, "\toverlaps_inversion_candidate");
    if (c->n_bp > 0)     fprintf(f, "\toverlaps_breakpoint");
    if (c->n_te > 0)     fprintf(f, "\toverlaps_TE_hotspot");
    fprintf(f, "\tnotes\n");

    for (int i = 0; i < c->nr; i++) {
        const Region* r = &c->regs[i];
        if (hdr_only && !r->contains_strong_window) continue;
        double* tmp = (double*)malloc((size_t)r->n_windows * sizeof(double));
        memcpy(tmp, r->values, (size_t)r->n_windows * sizeof(double));
        qsort(tmp, (size_t)r->n_windows, sizeof(double), dbl_cmp_asc);
        double med = (r->n_windows % 2)
                   ? tmp[r->n_windows / 2]
                   : 0.5 * (tmp[r->n_windows / 2 - 1] + tmp[r->n_windows / 2]);
        free(tmp);
        double mean = r->val_sum / r->n_windows;
        long len = r->end - r->start + 1;
        const char* cls = r->contains_strong_window ? "HDR" : "OUTLIER";
        fprintf(f, "%s_%s_%ld_%ld\t%s\t%ld\t%ld\t%ld\t%d",
                c->comp_id, r->chrom, r->start, r->end,
                r->chrom, r->start, r->end, len, r->n_windows);
        if (r->n_variants_total > 0) fprintf(f, "\t%ld", r->n_variants_total);
        else                          fprintf(f, "\tNA");
        fprintf(f, "\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%s\t%s",
                mean, med, r->val_max, r->val_min,
                c->threshold_used, c->stronger_threshold,
                r->contains_strong_window ? "TRUE" : "FALSE", cls);
        if (c->n_genes > 0) {
            long mid = (r->start + r->end) / 2;
            fprintf(f, "\t%s", nearest_gene(c->genes, c->n_genes, r->chrom, mid));
        }
        if (c->n_inv > 0) fprintf(f, "\t%s",
            bed_overlaps(c->inv_bed, c->n_inv, r->chrom, r->start, r->end) ? "TRUE" : "FALSE");
        if (c->n_bp > 0)  fprintf(f, "\t%s",
            bed_overlaps(c->bp_bed, c->n_bp, r->chrom, r->start, r->end) ? "TRUE" : "FALSE");
        if (c->n_te > 0)  fprintf(f, "\t%s",
            bed_overlaps(c->teh_bed, c->n_te, r->chrom, r->start, r->end) ? "TRUE" : "FALSE");
        fprintf(f, "\t\n");
    }
}

// ── Usage ──────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "outlier_scan — genome-wide outlier + HDR scan, statistic-agnostic.\n"
        "\n"
        "Input:\n"
        "  --input <windows.tsv>          Required.\n"
        "  --chrom_col / --start_col / --end_col / --nvar_col   (auto-detected if omitted)\n"
        "  --stat_col <name>              Numeric column to scan. (Aliases: --fst_col, --value_col)\n"
        "                                  If omitted, auto-detects fst_*, pi_*, hudson_fst, etc.\n"
        "  --stat_name <s>                Output column / JSON-key prefix (default 'stat').\n"
        "                                  E.g. 'fst', 'pi', 'fdM', 'tajimaD'.\n"
        "\n"
        "Direction-aware thresholding:\n"
        "  --direction high|low|abs       (default high)  high: v>=thr   low: v<=thr   abs: |v|>=thr\n"
        "  --threshold F                  Cutoff. (Alias: --fst_threshold.)\n"
        "  --top_pct P                    Use empirical percentile from data. Direction-aware:\n"
        "                                    high → 99.x-th pct;  low → bottom pct;  abs → top |v|.\n"
        "  --merge_gap_bp B               Merge flagged windows ≤ B bp apart (default 10000).\n"
        "  --stronger_threshold F         HDR escalation threshold (defaults: high 0.30, low 0.50, abs 0.30).\n"
        "  --stronger_windows <tsv>       Optional. Stronger-window TSV (e.g. 10-kb fixed).\n"
        "                                  HDR = overlap with at least one stronger-window row\n"
        "                                  passing --stronger_threshold under same --direction.\n"
        "                                  Without it: HDR uses region max (high/abs) or min (low).\n"
        "\n"
        "Annotation (all optional):\n"
        "  --genes_tsv <f>                chrom<TAB>start<TAB>end<TAB>gene_id\n"
        "  --inversion_bed <f>            adds overlaps_inversion_candidate column\n"
        "  --breakpoint_bed <f>           adds overlaps_breakpoint column\n"
        "  --te_bed <f>                   adds overlaps_TE_hotspot column\n"
        "\n"
        "Comparison metadata (for JSON summary):\n"
        "  --comparison_id <s>            --group_a <s>  --group_b <s>\n"
        "  --window_mode <s>              free-text (e.g. 15_variant, fixed_bp_10kb)\n"
        "  --genome_bp N                  used for percent_genome_*\n"
        "\n"
        "Outputs:\n"
        "  --out_windows <f>              flagged input windows above/below/|>=| threshold\n"
        "  --out_regions <f>              merged regions, OUTLIER or HDR\n"
        "  --out_hdr     <f>              HDR class only\n"
        "  --out_summary <f>              JSON summary\n");
}

// ── Main ───────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    const char *in_path = NULL;
    const char *chrom_col_name = NULL, *start_col_name = NULL,
               *end_col_name = NULL, *stat_col_name = NULL, *nvar_col_name = NULL;
    const char *stat_name = "stat";
    Direction direction = DIR_HIGH;
    double threshold = 0.25;
    int threshold_set = 0;
    double stronger_threshold = NAN;     // direction-dependent default applied later
    double top_pct = -1;
    long merge_gap_bp = 10000;
    const char *stronger_windows_path = NULL;
    const char *genes_path = NULL, *inversion_bed = NULL,
               *breakpoint_bed = NULL, *te_bed = NULL;
    const char *out_windows = NULL, *out_regions = NULL,
               *out_hdr = NULL, *out_summary = NULL;
    const char *comp_id = "comparison", *group_a = "A", *group_b = "B",
               *window_mode = "unknown";
    long genome_bp = -1;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--input")              && i+1<argc) in_path = argv[++i];
        else if (!strcmp(argv[i], "--chrom_col")          && i+1<argc) chrom_col_name = argv[++i];
        else if (!strcmp(argv[i], "--start_col")          && i+1<argc) start_col_name = argv[++i];
        else if (!strcmp(argv[i], "--end_col")            && i+1<argc) end_col_name = argv[++i];
        else if ((!strcmp(argv[i], "--stat_col") ||
                  !strcmp(argv[i], "--fst_col") ||
                  !strcmp(argv[i], "--value_col")) && i+1<argc) stat_col_name = argv[++i];
        else if (!strcmp(argv[i], "--stat_name")          && i+1<argc) stat_name = argv[++i];
        else if (!strcmp(argv[i], "--nvar_col")           && i+1<argc) nvar_col_name = argv[++i];
        else if (!strcmp(argv[i], "--direction")          && i+1<argc) {
            if (!parse_direction(argv[++i], &direction)) {
                fprintf(stderr, "[outlier_scan] --direction must be high|low|abs\n"); return 1;
            }
        }
        else if ((!strcmp(argv[i], "--threshold") ||
                  !strcmp(argv[i], "--fst_threshold")) && i+1<argc) { threshold = atof(argv[++i]); threshold_set = 1; }
        else if (!strcmp(argv[i], "--top_pct")            && i+1<argc) top_pct = atof(argv[++i]);
        else if (!strcmp(argv[i], "--merge_gap_bp")       && i+1<argc) merge_gap_bp = atol(argv[++i]);
        else if (!strcmp(argv[i], "--stronger_threshold") && i+1<argc) stronger_threshold = atof(argv[++i]);
        else if (!strcmp(argv[i], "--stronger_windows")   && i+1<argc) stronger_windows_path = argv[++i];
        else if (!strcmp(argv[i], "--genes_tsv")          && i+1<argc) genes_path = argv[++i];
        else if (!strcmp(argv[i], "--inversion_bed")      && i+1<argc) inversion_bed = argv[++i];
        else if (!strcmp(argv[i], "--breakpoint_bed")     && i+1<argc) breakpoint_bed = argv[++i];
        else if (!strcmp(argv[i], "--te_bed")             && i+1<argc) te_bed = argv[++i];
        else if (!strcmp(argv[i], "--out_windows")        && i+1<argc) out_windows = argv[++i];
        else if (!strcmp(argv[i], "--out_regions")        && i+1<argc) out_regions = argv[++i];
        else if (!strcmp(argv[i], "--out_hdr")            && i+1<argc) out_hdr = argv[++i];
        else if (!strcmp(argv[i], "--out_summary")        && i+1<argc) out_summary = argv[++i];
        else if (!strcmp(argv[i], "--comparison_id")      && i+1<argc) comp_id = argv[++i];
        else if (!strcmp(argv[i], "--group_a")            && i+1<argc) group_a = argv[++i];
        else if (!strcmp(argv[i], "--group_b")            && i+1<argc) group_b = argv[++i];
        else if (!strcmp(argv[i], "--window_mode")        && i+1<argc) window_mode = argv[++i];
        else if (!strcmp(argv[i], "--genome_bp")          && i+1<argc) genome_bp = atol(argv[++i]);
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[outlier_scan] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }
    if (!in_path) { print_usage(); return 1; }

    // Apply direction-dependent default for stronger_threshold if not set.
    if (isnan(stronger_threshold)) {
        if      (direction == DIR_HIGH) stronger_threshold = 0.30;
        else if (direction == DIR_LOW)  stronger_threshold = (threshold_set ? threshold * 0.5 : 0.10);
        else                             stronger_threshold = 0.30;
    }

    Tsv t;
    if (!tsv_load(in_path, &t)) { fprintf(stderr, "[outlier_scan] empty input\n"); return 2; }
    fprintf(stderr, "[outlier_scan] Loaded %d rows × %d cols\n", t.n_rows, t.n_cols);

    static const char* CHROM_ALIASES[] = {"chrom", "chromosome", "chr", "contig", NULL};
    static const char* START_ALIASES[] = {"start", "win_start", "window_start", NULL};
    static const char* END_ALIASES[]   = {"end", "win_end", "window_end", NULL};
    static const char* STAT_ALIASES[]  = {"fst", "FST", "Fst", "hudson_fst", "fst_mean",
                                            "pi", "theta_pi", "theta_pi_all",
                                            "dXY", "dxy", "dA",
                                            "Tajima_D", "tajima_D", "tajimaD",
                                            "fdM", "xpehh_mean", "iHS_mean", NULL};
    static const char* NVAR_ALIASES[]  = {"n_variants", "n_sites", "n_snps", NULL};

    int chrom_idx = chrom_col_name ? tsv_col_idx(&t, chrom_col_name) : tsv_find_col(&t, CHROM_ALIASES);
    int start_idx = start_col_name ? tsv_col_idx(&t, start_col_name) : tsv_find_col(&t, START_ALIASES);
    int end_idx   = end_col_name   ? tsv_col_idx(&t, end_col_name)   : tsv_find_col(&t, END_ALIASES);
    int stat_idx  = stat_col_name  ? tsv_col_idx(&t, stat_col_name)  : tsv_find_col(&t, STAT_ALIASES);
    int nvar_idx  = nvar_col_name  ? tsv_col_idx(&t, nvar_col_name)  : tsv_find_col(&t, NVAR_ALIASES);

    if (stat_idx < 0 && !stat_col_name) {
        // fall back to first column starting with fst_/pi_/dxy_/dA_
        for (int j = 0; j < t.n_cols; j++) {
            const char* c = t.col_names[j];
            if (!strncasecmp(c, "fst_", 4) || !strncasecmp(c, "pi_", 3) ||
                !strncasecmp(c, "dxy_", 4) || !strncasecmp(c, "dA_", 3) ||
                !strncasecmp(c, "fdM_", 4)) { stat_idx = j; break; }
        }
    }
    if (chrom_idx < 0 || start_idx < 0 || end_idx < 0 || stat_idx < 0) {
        fprintf(stderr, "[outlier_scan] Missing required column. "
                "chrom_idx=%d start_idx=%d end_idx=%d stat_idx=%d\n",
                chrom_idx, start_idx, end_idx, stat_idx);
        tsv_free(&t); return 3;
    }
    fprintf(stderr, "[outlier_scan] Using cols: chrom=%s start=%s end=%s stat=%s%s%s; "
            "direction=%s stat_name=%s\n",
            t.col_names[chrom_idx], t.col_names[start_idx], t.col_names[end_idx],
            t.col_names[stat_idx],
            nvar_idx >= 0 ? " n_variants=" : "",
            nvar_idx >= 0 ? t.col_names[nvar_idx] : "",
            direction_str(direction), stat_name);

    // Parse into Window array.
    Window* w = (Window*)malloc((size_t)t.n_rows * sizeof(Window));
    int nw = 0;
    for (int i = 0; i < t.n_rows; i++) {
        double v; long s, e;
        if (!parse_d(t.cells[i][stat_idx], &v)) continue;
        if (!parse_i(t.cells[i][start_idx], &s)) continue;
        if (!parse_i(t.cells[i][end_idx], &e)) continue;
        w[nw].chrom = strdup(t.cells[i][chrom_idx]);
        w[nw].start = s; w[nw].end = e;
        w[nw].value = v;
        w[nw].n_variants = -1;
        if (nvar_idx >= 0) {
            long n; if (parse_i(t.cells[i][nvar_idx], &n)) w[nw].n_variants = n;
        }
        w[nw].flagged = 0;
        nw++;
    }
    qsort(w, (size_t)nw, sizeof(Window), win_cmp);

    // Threshold resolution (--top_pct overrides --threshold).
    double threshold_used = threshold;
    const char* threshold_source = "absolute";
    if (top_pct >= 0 && top_pct <= 100) {
        threshold_source = "top_pct";
        double* arr = (double*)malloc((size_t)nw * sizeof(double));
        int m = 0;
        for (int i = 0; i < nw; i++) {
            if (!isfinite(w[i].value)) continue;
            if (direction == DIR_ABS) arr[m++] = fabs(w[i].value);
            else                       arr[m++] = w[i].value;
        }
        double q;
        if (direction == DIR_HIGH || direction == DIR_ABS) q = 1.0 - top_pct / 100.0;
        else                                                 q = top_pct / 100.0;
        threshold_used = percentile_value(arr, m, q);
        free(arr);
        fprintf(stderr, "[outlier_scan] --top_pct %.4g (dir=%s) → threshold = %.6g\n",
                top_pct, direction_str(direction), threshold_used);
    } else {
        fprintf(stderr, "[outlier_scan] threshold = %.6g (dir=%s)\n",
                threshold_used, direction_str(direction));
    }

    int n_flagged = 0;
    for (int i = 0; i < nw; i++)
        if (passes_outlier(w[i].value, threshold_used, direction)) { w[i].flagged = 1; n_flagged++; }
    fprintf(stderr, "[outlier_scan] flagged %d / %d windows\n", n_flagged, nw);

    // Optional stronger windows for HDR.
    Window* sw = NULL;
    int nsw = 0;
    if (stronger_windows_path) {
        Tsv ts;
        if (tsv_load(stronger_windows_path, &ts)) {
            int sc = chrom_col_name ? tsv_col_idx(&ts, chrom_col_name) : tsv_find_col(&ts, CHROM_ALIASES);
            int ss = start_col_name ? tsv_col_idx(&ts, start_col_name) : tsv_find_col(&ts, START_ALIASES);
            int se = end_col_name   ? tsv_col_idx(&ts, end_col_name)   : tsv_find_col(&ts, END_ALIASES);
            int sf = stat_col_name  ? tsv_col_idx(&ts, stat_col_name)  : tsv_find_col(&ts, STAT_ALIASES);
            if (sf < 0 && !stat_col_name) {
                for (int j = 0; j < ts.n_cols; j++) {
                    const char* c = ts.col_names[j];
                    if (!strncasecmp(c, "fst_", 4) || !strncasecmp(c, "pi_", 3)) { sf = j; break; }
                }
            }
            if (sc < 0 || ss < 0 || se < 0 || sf < 0) {
                fprintf(stderr, "[outlier_scan] --stronger_windows: missing required column\n");
            } else {
                sw = (Window*)malloc((size_t)ts.n_rows * sizeof(Window));
                for (int i = 0; i < ts.n_rows; i++) {
                    double v; long s, e;
                    if (!parse_d(ts.cells[i][sf], &v)) continue;
                    if (!parse_i(ts.cells[i][ss], &s)) continue;
                    if (!parse_i(ts.cells[i][se], &e)) continue;
                    if (!passes_outlier(v, stronger_threshold, direction)) continue;
                    sw[nsw].chrom = strdup(ts.cells[i][sc]);
                    sw[nsw].start = s; sw[nsw].end = e; sw[nsw].value = v;
                    sw[nsw].n_variants = -1; sw[nsw].flagged = 1;
                    nsw++;
                }
                fprintf(stderr, "[outlier_scan] stronger windows passing %.4g (%s): %d\n",
                        stronger_threshold, direction_str(direction), nsw);
            }
            tsv_free(&ts);
        }
    }

    // Merge flagged windows.
    Region* regs = (Region*)calloc((size_t)nw + 1, sizeof(Region));
    int nr = 0;
    for (int i = 0; i < nw; i++) {
        if (!w[i].flagged) continue;
        if (nr > 0
            && strcmp(regs[nr-1].chrom, w[i].chrom) == 0
            && (w[i].start - regs[nr-1].end) <= merge_gap_bp) {
            if (w[i].end > regs[nr-1].end) regs[nr-1].end = w[i].end;
            if (w[i].n_variants >= 0) regs[nr-1].n_variants_total += w[i].n_variants;
            region_push(&regs[nr-1], w[i].value);
        } else {
            strncpy(regs[nr].chrom, w[i].chrom, sizeof(regs[nr].chrom) - 1);
            regs[nr].chrom[sizeof(regs[nr].chrom) - 1] = 0;
            regs[nr].start = w[i].start; regs[nr].end = w[i].end;
            regs[nr].n_windows = 0;
            regs[nr].n_variants_total = (w[i].n_variants >= 0) ? w[i].n_variants : 0;
            regs[nr].val_sum = 0;
            regs[nr].val_min = INFINITY; regs[nr].val_max = -INFINITY;
            regs[nr].contains_strong_window = 0;
            region_push(&regs[nr], w[i].value);
            nr++;
        }
    }
    fprintf(stderr, "[outlier_scan] merged into %d outlier regions\n", nr);

    // HDR check (direction-aware fallback uses min for DIR_LOW, max for HIGH/ABS).
    for (int i = 0; i < nr; i++) {
        if (sw && nsw > 0) {
            for (int j = 0; j < nsw; j++) {
                if (strcmp(sw[j].chrom, regs[i].chrom) != 0) continue;
                if (sw[j].end < regs[i].start) continue;
                if (sw[j].start > regs[i].end) continue;
                regs[i].contains_strong_window = 1;
                break;
            }
        } else {
            double probe = (direction == DIR_LOW) ? regs[i].val_min : regs[i].val_max;
            regs[i].contains_strong_window = passes_outlier(probe, stronger_threshold, direction);
        }
    }

    int n_hdr = 0;
    long total_outlier_bp = 0, total_hdr_bp = 0;
    for (int i = 0; i < nr; i++) {
        long len = regs[i].end - regs[i].start + 1;
        total_outlier_bp += len;
        if (regs[i].contains_strong_window) { n_hdr++; total_hdr_bp += len; }
    }

    // Annotation loads
    BedRec *inv_bed = NULL, *bp_bed = NULL, *teh_bed = NULL;
    int n_inv = 0, n_bp = 0, n_te = 0;
    if (inversion_bed)  n_inv = bed_load(inversion_bed, &inv_bed);
    if (breakpoint_bed) n_bp  = bed_load(breakpoint_bed, &bp_bed);
    if (te_bed)          n_te = bed_load(te_bed, &teh_bed);

    GeneRec *genes = NULL; int n_genes = 0;
    if (genes_path) n_genes = genes_load(genes_path, &genes);

    // ── Output: flagged windows ──
    if (out_windows) {
        FILE* f = fopen(out_windows, "w");
        if (f) {
            fprintf(f, "# schema_version=%s kind=outlier_windows comparison_id=%s stat=%s direction=%s "
                       "threshold=%.6g\n", SCHEMA_VERSION, comp_id, stat_name,
                       direction_str(direction), threshold_used);
            fprintf(f, "chrom\tstart\tend\t%s\tn_variants\tflagged\n", stat_name);
            for (int i = 0; i < nw; i++) {
                if (!w[i].flagged) continue;
                fprintf(f, "%s\t%ld\t%ld\t%.6g", w[i].chrom, w[i].start, w[i].end, w[i].value);
                if (w[i].n_variants >= 0) fprintf(f, "\t%ld", w[i].n_variants); else fprintf(f, "\tNA");
                fprintf(f, "\t1\n");
            }
            fclose(f);
        }
    }

    // ── Regions ──
    RegionsEmitCtx ctx = {
        .regs = regs, .nr = nr, .comp_id = comp_id, .stat_name = stat_name,
        .direction = direction,
        .threshold_used = threshold_used, .stronger_threshold = stronger_threshold,
        .merge_gap_bp = merge_gap_bp,
        .genes = genes, .n_genes = n_genes,
        .inv_bed = inv_bed, .n_inv = n_inv,
        .bp_bed = bp_bed,   .n_bp = n_bp,
        .teh_bed = teh_bed, .n_te = n_te,
    };
    if (out_regions) { FILE* f = fopen(out_regions, "w"); if (f) { emit_regions(f, &ctx, 0); fclose(f); } }
    if (out_hdr)     { FILE* f = fopen(out_hdr, "w");     if (f) { emit_regions(f, &ctx, 1); fclose(f); } }

    // ── JSON summary ──
    if (out_summary) {
        FILE* f = fopen(out_summary, "w");
        if (f) {
            double pct_out = (genome_bp > 0) ? 100.0 * total_outlier_bp / genome_bp : NAN;
            double pct_hdr = (genome_bp > 0) ? 100.0 * total_hdr_bp     / genome_bp : NAN;
            fprintf(f, "{\n");
            fprintf(f, "  \"schema_version\": \"%s\",\n", SCHEMA_VERSION);
            fprintf(f, "  \"comparison_id\": \"%s\",\n", comp_id);
            fprintf(f, "  \"group_a\": \"%s\",\n", group_a);
            fprintf(f, "  \"group_b\": \"%s\",\n", group_b);
            fprintf(f, "  \"window_mode\": \"%s\",\n", window_mode);
            fprintf(f, "  \"stat_name\": \"%s\",\n", stat_name);
            fprintf(f, "  \"direction\": \"%s\",\n", direction_str(direction));
            fprintf(f, "  \"threshold\": %.6g,\n", threshold_used);
            fprintf(f, "  \"threshold_source\": \"%s\",\n", threshold_source);
            if (top_pct >= 0) fprintf(f, "  \"top_pct\": %.6g,\n", top_pct);
            fprintf(f, "  \"stronger_threshold\": %.6g,\n", stronger_threshold);
            fprintf(f, "  \"merge_gap_bp\": %ld,\n", merge_gap_bp);
            fprintf(f, "  \"stronger_windows_used\": %s,\n", (sw && nsw > 0) ? "true" : "false");
            fprintf(f, "  \"n_input_windows\": %d,\n", nw);
            fprintf(f, "  \"n_outlier_windows\": %d,\n", n_flagged);
            fprintf(f, "  \"n_outlier_regions\": %d,\n", nr);
            fprintf(f, "  \"n_HDR_regions\": %d,\n", n_hdr);
            fprintf(f, "  \"total_outlier_bp\": %ld,\n", total_outlier_bp);
            fprintf(f, "  \"total_HDR_bp\": %ld,\n", total_hdr_bp);
            if (genome_bp > 0) {
                fprintf(f, "  \"genome_bp\": %ld,\n", genome_bp);
                fprintf(f, "  \"percent_genome_outlier\": %.6g,\n", pct_out);
                fprintf(f, "  \"percent_genome_HDR\": %.6g,\n", pct_hdr);
            }
            fprintf(f, "  \"note\": \"Empirical %s outliers; not a neutral-calibrated selection scan.\"\n",
                    stat_name);
            fprintf(f, "}\n");
            fclose(f);
        }
    }

    // Cleanup
    for (int i = 0; i < nw; i++) free(w[i].chrom);
    free(w);
    if (sw) { for (int i = 0; i < nsw; i++) free(sw[i].chrom); free(sw); }
    for (int i = 0; i < nr; i++) free(regs[i].values);
    free(regs);
    if (inv_bed) bed_free(inv_bed, n_inv);
    if (bp_bed)  bed_free(bp_bed, n_bp);
    if (teh_bed) bed_free(teh_bed, n_te);
    if (genes)   genes_free(genes, n_genes);
    tsv_free(&t);
    return 0;
}
