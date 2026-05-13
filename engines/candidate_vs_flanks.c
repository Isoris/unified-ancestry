// =============================================================================
// candidate_vs_flanks.c — inversion-candidate support track from any per-window
//                        statistic TSV. Sibling to outlier_scan.c:
//
//   outlier_scan          : genome-wide scan — "where are the special regions?"
//   candidate_vs_flanks   : candidate-support — "is THIS candidate different
//                                                from its local flanks?"
//
// For each candidate interval and each requested stat column, computes:
//   - mean / median inside candidate vs in flanking windows
//   - Mann-Whitney U (inside vs flanks), tie-corrected z, normal-approx p
//   - Optional permutation test (random-relabel inside/flank, count
//     |U_perm - mu| >= |U_obs - mu|), +1 smoothing
//   - inside-minus-flank effect direction
//
// Inputs:
//   --candidates <bed>            chrom start end [candidate_id]
//   --windows_tsv <tsv>           per-window stat table (same column conventions
//                                  as outlier_scan: chrom/start/end + one or more
//                                  numeric stat columns)
//   --stat_cols a,b,c             column(s) to test
//   --stat_names f,p,d            optional display names (same length); default = col names
//   --chrom_col / --start_col / --end_col   (auto-detected if omitted)
//   --flank_bp N                  per-side flank size (default 300000)
//   --exclude_flank_overlaps      drop flank windows that fall in OTHER candidates
//   --permutations N              perm reps (default 1000; 0 disables)
//   --seed K                      RNG seed
//   --ncores N                    OpenMP across (candidate × stat) tasks
//   --out_table   <f>             per-(candidate,stat) TSV
//   --out_summary <f>             JSON summary
//
// Compile: gcc -O3 -march=native -fopenmp -o candidate_vs_flanks candidate_vs_flanks.c -lm
// =============================================================================

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define SCHEMA_VERSION "candidate_vs_flanks_v1"

// ── Tsv ─────────────────────────────────────────────────────────────────────

typedef struct {
    char** col_names;
    int n_cols;
    char*** cells;
    int n_rows;
    int cap;
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
    if (!f) { fprintf(stderr, "[cvf] Cannot open %s\n", path); return 0; }
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

static char** split_csv(const char* s, int* n_out) {
    int cap = 8, n = 0;
    char** out = (char**)malloc((size_t)cap * sizeof(char*));
    const char* p = s;
    while (1) {
        const char* c = strchr(p, ',');
        size_t len = c ? (size_t)(c - p) : strlen(p);
        char* tok = (char*)malloc(len + 1);
        memcpy(tok, p, len); tok[len] = 0;
        if (n == cap) { cap *= 2; out = (char**)realloc(out, (size_t)cap * sizeof(char*)); }
        out[n++] = tok;
        if (!c) break;
        p = c + 1;
    }
    *n_out = n;
    return out;
}

// ── BED candidates ─────────────────────────────────────────────────────────

typedef struct {
    char chrom[128];
    long start, end;
    char id[128];
} Candidate;

static int candidates_load(const char* path, Candidate** out) {
    *out = NULL;
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[cvf] Cannot open candidates %s\n", path); return 0; }
    int cap = 256, n = 0;
    Candidate* arr = (Candidate*)malloc((size_t)cap * sizeof(Candidate));
    char* line = NULL; size_t lcap = 0; ssize_t got;
    while ((got = getline(&line, &lcap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        char ch[128] = "", id[128] = "";
        long s, e;
        int nf = sscanf(line, "%127s %ld %ld %127s", ch, &s, &e, id);
        if (nf < 3) continue;
        if (n == cap) { cap *= 2; arr = (Candidate*)realloc(arr, (size_t)cap * sizeof(Candidate)); }
        strncpy(arr[n].chrom, ch, sizeof(arr[n].chrom) - 1); arr[n].chrom[sizeof(arr[n].chrom) - 1] = 0;
        arr[n].start = s; arr[n].end = e;
        if (nf >= 4) {
            strncpy(arr[n].id, id, sizeof(arr[n].id) - 1); arr[n].id[sizeof(arr[n].id) - 1] = 0;
        } else {
            snprintf(arr[n].id, sizeof(arr[n].id), "cand_%d", n + 1);
        }
        n++;
    }
    free(line); fclose(f);
    *out = arr;
    return n;
}

// True if [s, e] on chrom overlaps any candidate.
static int candidate_overlaps_any(const Candidate* arr, int n, const char* chrom,
                                   long s, long e, int skip_idx) {
    for (int i = 0; i < n; i++) {
        if (i == skip_idx) continue;
        if (strcmp(arr[i].chrom, chrom) != 0) continue;
        if (arr[i].end < s) continue;
        if (arr[i].start > e) continue;
        return 1;
    }
    return 0;
}

// ── Wilcoxon + permutation ──────────────────────────────────────────────────

typedef struct { double v; int group; } VG;

static int vg_cmp(const void* a, const void* b) {
    const VG* va = (const VG*)a, *vb = (const VG*)b;
    if (va->v < vb->v) return -1;
    if (va->v > vb->v) return  1;
    return 0;
}

static void wilcoxon(VG* arr, int n, int n0, int n1,
                     double* W_out, double* z_out, double* p_out) {
    qsort(arr, (size_t)n, sizeof(VG), vg_cmp);
    double sum_rank0 = 0;
    double sum_tie = 0;
    int i = 0;
    while (i < n) {
        int k = i + 1;
        while (k < n && arr[k].v == arr[i].v) k++;
        double avg_rank = (double)(i + 1 + k) / 2.0;
        int tcount = k - i;
        if (tcount > 1) sum_tie += (double)tcount * tcount * tcount - (double)tcount;
        for (int j = i; j < k; j++) if (arr[j].group == 0) sum_rank0 += avg_rank;
        i = k;
    }
    double U = sum_rank0 - (double)n0 * (n0 + 1) / 2.0;
    double mu_U = (double)n0 * n1 / 2.0;
    double tie_term = (n > 1) ? sum_tie / ((double)n * (n - 1)) : 0;
    double var_U = ((double)n0 * n1 / 12.0) * ((double)(n + 1) - tie_term);
    *W_out = U;
    if (var_U > 0) {
        double z = (U - mu_U) / sqrt(var_U);
        *z_out = z;
        *p_out = erfc(fabs(z) * M_SQRT1_2);
    } else {
        *z_out = NAN; *p_out = NAN;
    }
}

static void perm_test(const VG* arr, int n, int n0, int n1,
                      double U_obs, int n_perm, unsigned int* rng,
                      double* p_out) {
    if (n_perm <= 0 || n0 <= 0 || n1 <= 0) { *p_out = NAN; return; }
    VG* sorted = (VG*)malloc((size_t)n * sizeof(VG));
    memcpy(sorted, arr, (size_t)n * sizeof(VG));
    qsort(sorted, (size_t)n, sizeof(VG), vg_cmp);
    double* rank = (double*)malloc((size_t)n * sizeof(double));
    int i = 0;
    while (i < n) {
        int k = i + 1;
        while (k < n && sorted[k].v == sorted[i].v) k++;
        double ar = (double)(i + 1 + k) / 2.0;
        for (int j = i; j < k; j++) rank[j] = ar;
        i = k;
    }
    double mu_U = (double)n0 * n1 / 2.0;
    double abs_obs = fabs(U_obs - mu_U);
    int* labels = (int*)malloc((size_t)n * sizeof(int));
    for (int j = 0; j < n; j++) labels[j] = (j < n0) ? 0 : 1;
    int ge = 0;
    for (int b = 0; b < n_perm; b++) {
        for (int j = n - 1; j > 0; j--) {
            int r = (int)(rand_r(rng) % (unsigned int)(j + 1));
            int tmp = labels[j]; labels[j] = labels[r]; labels[r] = tmp;
        }
        double sr = 0;
        for (int j = 0; j < n; j++) if (labels[j] == 0) sr += rank[j];
        double U = sr - (double)n0 * (n0 + 1) / 2.0;
        if (fabs(U - mu_U) >= abs_obs) ge++;
    }
    *p_out = (double)(ge + 1) / (double)(n_perm + 1);
    free(labels); free(rank); free(sorted);
}

static int dbl_cmp(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y ? 1 : 0);
}

static double mean_of(const double* v, int n) {
    if (n == 0) return NAN;
    double s = 0;
    for (int i = 0; i < n; i++) s += v[i];
    return s / n;
}
static double median_of(double* v, int n) {
    if (n == 0) return NAN;
    qsort(v, (size_t)n, sizeof(double), dbl_cmp);
    return (n % 2) ? v[n/2] : 0.5 * (v[n/2 - 1] + v[n/2]);
}

// ── Result struct ──────────────────────────────────────────────────────────

typedef struct {
    char* cand_id;
    char* chrom;
    long  cand_start, cand_end;
    char* stat_name;
    int   n_inside, n_flank;
    double mean_inside, mean_flank;
    double median_inside, median_flank;
    double effect_diff;
    double W, z, wilcoxon_p, perm_p;
    const char* direction;     // "inside_higher", "inside_lower", "same", "NA"
} CandResult;

// ── Main ───────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "candidate_vs_flanks — inversion-candidate support via Wilcoxon vs local flanks.\n"
        "\n"
        "  --candidates <bed>            chrom start end [candidate_id]\n"
        "  --windows_tsv <tsv>           per-window stat table\n"
        "  --stat_cols a,b,c             column(s) to test\n"
        "  --stat_names f,p,d            optional display names (same length)\n"
        "  --chrom_col / --start_col / --end_col   (auto-detected if omitted)\n"
        "  --flank_bp N                  per-side flank size (default 300000)\n"
        "  --exclude_flank_overlaps      drop flank windows in other candidates\n"
        "  --permutations N              default 1000; 0 disables\n"
        "  --seed K                      RNG seed\n"
        "  --ncores N                    OpenMP across (candidate × stat) tasks\n"
        "  --out_table   <f>             per-(candidate,stat) TSV\n"
        "  --out_summary <f>             JSON summary\n"
        "\n"
        "Output columns (schema_version=" SCHEMA_VERSION "):\n"
        "  candidate_id chrom cand_start cand_end length_bp flank_bp stat_name\n"
        "  n_inside n_flank mean_inside mean_flank median_inside median_flank\n"
        "  effect_diff inside_minus_flank_direction W z wilcoxon_p perm_p\n");
}

int main(int argc, char** argv) {
    const char *cand_path = NULL, *win_path = NULL;
    const char *stat_cols_spec = NULL, *stat_names_spec = NULL;
    const char *chrom_col_name = NULL, *start_col_name = NULL, *end_col_name = NULL;
    long flank_bp = 300000;
    int exclude_flank_overlaps = 0;
    int n_perm = 1000;
    unsigned int seed = (unsigned int)time(NULL);
    int seed_set = 0;
    int ncores = 1;
    const char *out_table = NULL, *out_summary = NULL;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--candidates")   && i+1<argc) cand_path = argv[++i];
        else if (!strcmp(argv[i], "--windows_tsv")  && i+1<argc) win_path  = argv[++i];
        else if (!strcmp(argv[i], "--stat_cols")    && i+1<argc) stat_cols_spec = argv[++i];
        else if (!strcmp(argv[i], "--stat_names")   && i+1<argc) stat_names_spec = argv[++i];
        else if (!strcmp(argv[i], "--chrom_col")    && i+1<argc) chrom_col_name = argv[++i];
        else if (!strcmp(argv[i], "--start_col")    && i+1<argc) start_col_name = argv[++i];
        else if (!strcmp(argv[i], "--end_col")      && i+1<argc) end_col_name   = argv[++i];
        else if (!strcmp(argv[i], "--flank_bp")     && i+1<argc) flank_bp = atol(argv[++i]);
        else if (!strcmp(argv[i], "--exclude_flank_overlaps"))   exclude_flank_overlaps = 1;
        else if (!strcmp(argv[i], "--permutations") && i+1<argc) n_perm = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed")         && i+1<argc) { seed = (unsigned int)atoi(argv[++i]); seed_set = 1; }
        else if (!strcmp(argv[i], "--ncores")       && i+1<argc) ncores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out_table")    && i+1<argc) out_table = argv[++i];
        else if (!strcmp(argv[i], "--out_summary")  && i+1<argc) out_summary = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[cvf] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }
    if (!cand_path || !win_path || !stat_cols_spec) { print_usage(); return 1; }

    #ifdef _OPENMP
    omp_set_num_threads(ncores > 0 ? ncores : 1);
    #else
    (void)ncores;
    #endif

    Candidate* cands = NULL;
    int n_cands = candidates_load(cand_path, &cands);
    fprintf(stderr, "[cvf] %d candidates loaded\n", n_cands);
    if (n_cands == 0) return 2;

    Tsv t;
    if (!tsv_load(win_path, &t)) { fprintf(stderr, "[cvf] empty windows TSV\n"); return 3; }
    fprintf(stderr, "[cvf] %d window rows × %d cols\n", t.n_rows, t.n_cols);

    static const char* CHROM_ALIASES[] = {"chrom", "chromosome", "chr", "contig", NULL};
    static const char* START_ALIASES[] = {"start", "win_start", "window_start", NULL};
    static const char* END_ALIASES[]   = {"end", "win_end", "window_end", NULL};

    int chrom_idx = chrom_col_name ? tsv_col_idx(&t, chrom_col_name) : tsv_find_col(&t, CHROM_ALIASES);
    int start_idx = start_col_name ? tsv_col_idx(&t, start_col_name) : tsv_find_col(&t, START_ALIASES);
    int end_idx   = end_col_name   ? tsv_col_idx(&t, end_col_name)   : tsv_find_col(&t, END_ALIASES);
    if (chrom_idx < 0 || start_idx < 0 || end_idx < 0) {
        fprintf(stderr, "[cvf] Missing chrom/start/end columns\n");
        tsv_free(&t); free(cands); return 4;
    }

    int nstat;
    char** stat_cols = split_csv(stat_cols_spec, &nstat);
    int* stat_idx = (int*)malloc((size_t)nstat * sizeof(int));
    for (int i = 0; i < nstat; i++) {
        stat_idx[i] = tsv_col_idx(&t, stat_cols[i]);
        if (stat_idx[i] < 0) {
            fprintf(stderr, "[cvf] stat_col '%s' not found in TSV\n", stat_cols[i]);
            return 5;
        }
    }
    char** stat_names = NULL;
    int nname = 0;
    if (stat_names_spec) {
        stat_names = split_csv(stat_names_spec, &nname);
        if (nname != nstat) {
            fprintf(stderr, "[cvf] --stat_names must have same length as --stat_cols\n");
            return 6;
        }
    }
    fprintf(stderr, "[cvf] %d candidates × %d stats = %d tasks; flank_bp=%ld permutations=%d seed=%u%s\n",
            n_cands, nstat, n_cands * nstat, flank_bp, n_perm, seed,
            seed_set ? "" : " (time-based)");

    int n_tasks = n_cands * nstat;
    CandResult* results = (CandResult*)calloc((size_t)n_tasks, sizeof(CandResult));

    #pragma omp parallel for schedule(dynamic)
    for (int task = 0; task < n_tasks; task++) {
        int c_idx = task / nstat;
        int s_idx = task % nstat;
        unsigned int local_seed = seed ^ (unsigned int)(task * 2654435761u);
        const Candidate* C = &cands[c_idx];
        int vi = stat_idx[s_idx];

        VG* arr = (VG*)malloc((size_t)t.n_rows * sizeof(VG));
        double* inside_vals = (double*)malloc((size_t)t.n_rows * sizeof(double));
        double* flank_vals  = (double*)malloc((size_t)t.n_rows * sizeof(double));
        int n_inside = 0, n_flank = 0, n_total = 0;

        long flank_lo_start = C->start - flank_bp;
        long flank_lo_end   = C->start - 1;
        long flank_hi_start = C->end + 1;
        long flank_hi_end   = C->end + flank_bp;

        for (int i = 0; i < t.n_rows; i++) {
            if (strcmp(t.cells[i][chrom_idx], C->chrom) != 0) continue;
            long ws, we;
            if (!parse_i(t.cells[i][start_idx], &ws)) continue;
            if (!parse_i(t.cells[i][end_idx], &we))   continue;
            long mid = (ws + we) / 2;
            double v;
            if (!parse_d(t.cells[i][vi], &v)) continue;
            int is_inside = (mid >= C->start && mid <= C->end);
            int is_flank  = ((mid >= flank_lo_start && mid <= flank_lo_end) ||
                              (mid >= flank_hi_start && mid <= flank_hi_end));
            if (!is_inside && !is_flank) continue;
            if (is_flank && exclude_flank_overlaps &&
                candidate_overlaps_any(cands, n_cands, C->chrom, ws, we, c_idx))
                continue;
            if (is_inside) {
                inside_vals[n_inside++] = v;
                arr[n_total].v = v; arr[n_total].group = 0; n_total++;
            } else {
                flank_vals[n_flank++] = v;
                arr[n_total].v = v; arr[n_total].group = 1; n_total++;
            }
        }

        CandResult* R = &results[task];
        R->cand_id = strdup(C->id);
        R->chrom = strdup(C->chrom);
        R->cand_start = C->start;
        R->cand_end = C->end;
        R->stat_name = strdup(stat_names ? stat_names[s_idx] : stat_cols[s_idx]);
        R->n_inside = n_inside; R->n_flank = n_flank;
        R->mean_inside = mean_of(inside_vals, n_inside);
        R->mean_flank  = mean_of(flank_vals, n_flank);
        R->median_inside = median_of(inside_vals, n_inside);
        R->median_flank  = median_of(flank_vals,  n_flank);
        R->effect_diff   = R->mean_inside - R->mean_flank;
        R->W = R->z = R->wilcoxon_p = R->perm_p = NAN;
        R->direction = "NA";

        if (n_inside >= 2 && n_flank >= 2) {
            wilcoxon(arr, n_total, n_inside, n_flank, &R->W, &R->z, &R->wilcoxon_p);
            if (n_perm > 0) perm_test(arr, n_total, n_inside, n_flank, R->W, n_perm, &local_seed, &R->perm_p);
            if (R->mean_inside > R->mean_flank)      R->direction = "inside_higher";
            else if (R->mean_inside < R->mean_flank) R->direction = "inside_lower";
            else                                      R->direction = "same";
        }
        free(arr); free(inside_vals); free(flank_vals);
    }

    // ── Output: per-(candidate, stat) TSV ──
    FILE* fout = out_table ? fopen(out_table, "w") : stdout;
    if (!fout) { fprintf(stderr, "[cvf] Cannot open --out_table %s\n", out_table); return 7; }
    fprintf(fout, "# schema_version=%s flank_bp=%ld permutations=%d\n",
            SCHEMA_VERSION, flank_bp, n_perm);
    fprintf(fout, "candidate_id\tchrom\tcand_start\tcand_end\tlength_bp\tflank_bp\tstat_name\t"
                   "n_inside\tn_flank\tmean_inside\tmean_flank\tmedian_inside\tmedian_flank\t"
                   "effect_diff\tinside_minus_flank_direction\tW\tz\twilcoxon_p\tperm_p\n");
    for (int i = 0; i < n_tasks; i++) {
        CandResult* R = &results[i];
        fprintf(fout, "%s\t%s\t%ld\t%ld\t%ld\t%ld\t%s\t%d\t%d",
                R->cand_id, R->chrom, R->cand_start, R->cand_end,
                R->cand_end - R->cand_start + 1, flank_bp,
                R->stat_name, R->n_inside, R->n_flank);
        #define EM(v) do { if (isnan(v)) fprintf(fout, "\tNA"); else fprintf(fout, "\t%.6g", v); } while (0)
        EM(R->mean_inside); EM(R->mean_flank);
        EM(R->median_inside); EM(R->median_flank);
        EM(R->effect_diff);
        fprintf(fout, "\t%s", R->direction);
        EM(R->W); EM(R->z); EM(R->wilcoxon_p); EM(R->perm_p);
        #undef EM
        fputc('\n', fout);
    }
    if (out_table) fclose(fout);

    // ── JSON summary ──
    if (out_summary) {
        FILE* f = fopen(out_summary, "w");
        if (f) {
            int n_sig = 0;
            for (int i = 0; i < n_tasks; i++)
                if (isfinite(results[i].wilcoxon_p) && results[i].wilcoxon_p < 0.05) n_sig++;
            fprintf(f, "{\n");
            fprintf(f, "  \"schema_version\": \"%s\",\n", SCHEMA_VERSION);
            fprintf(f, "  \"n_candidates\": %d,\n", n_cands);
            fprintf(f, "  \"n_stats\": %d,\n", nstat);
            fprintf(f, "  \"n_tests\": %d,\n", n_tasks);
            fprintf(f, "  \"flank_bp\": %ld,\n", flank_bp);
            fprintf(f, "  \"exclude_flank_overlaps\": %s,\n", exclude_flank_overlaps ? "true" : "false");
            fprintf(f, "  \"permutations\": %d,\n", n_perm);
            fprintf(f, "  \"n_sig_wilcoxon_p_lt_0.05\": %d,\n", n_sig);
            fprintf(f, "  \"note\": \"Inside-vs-flanks Wilcoxon. Not a calibrated selection test; "
                       "use multiple-testing correction across candidates × stats.\"\n");
            fprintf(f, "}\n");
            fclose(f);
        }
    }

    for (int i = 0; i < n_tasks; i++) {
        free(results[i].cand_id); free(results[i].chrom); free(results[i].stat_name);
    }
    free(results);
    for (int i = 0; i < nstat; i++) free(stat_cols[i]);
    free(stat_cols);
    if (stat_names) { for (int i = 0; i < nname; i++) free(stat_names[i]); free(stat_names); }
    free(stat_idx);
    free(cands);
    tsv_free(&t);
    return 0;
}
