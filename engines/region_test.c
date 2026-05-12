// =============================================================================
// region_test.c — per-region vs rest-of-genome Wilcoxon + permutation tests
//                  on arbitrary value columns of a TSV.
//
// Driver: Bonhomme et al. methods excerpt — "Wilcoxon tests or resampling
// procedures (1000 tests) to compare diversity estimates (π, πN, πS) between
// candidate regions and the rest of the genome."
//
// Generic enough to test any numeric column. Typical use: feed it the output
// of pi_NS (one row per locus × group) and ask which inversion regions show
// significantly different πN/πS, π0/π4, etc., from rest of genome.
//
// Input:
//   --input <tsv>            TSV with header. NA/Inf cells dropped per test.
//   --region_col <name>      Partition column (default: region).
//   --value_cols a,b,c       One or more numeric columns to test (required).
//   --rest_label <name>      Value of region_col denoting "the rest". If not
//                            given, the rest is "any row whose region_col is
//                            different from the region being tested".
//   --skip_regions a,b       Regions to skip entirely (not tested, not in rest).
//   --permutations N         Permutation test reps (default 1000; 0 disables).
//   --seed K                 RNG seed (default time(NULL)).
//   --ncores N               OpenMP across (region × value_col) tests.
//   --out <f>                Output TSV.
//
// Output (one row per region × value_col):
//   region value_col n_in n_out mean_in mean_out median_in median_out
//   W z wilcoxon_p perm_p effect_diff
//
// Compile: gcc -O3 -march=native -fopenmp -o region_test region_test.c -lm
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

#define SCHEMA_VERSION "region_test_v1"

// ── TSV in-memory representation ─────────────────────────────────────────────

typedef struct {
    char** col_names;
    int    n_cols;
    char*** cells;       // [row][col]
    int    n_rows;
    int    cap;
} Tsv;

static void tsv_free(Tsv* t) {
    if (!t) return;
    for (int j = 0; j < t->n_cols; j++) free(t->col_names[j]);
    free(t->col_names);
    for (int i = 0; i < t->n_rows; i++) {
        for (int j = 0; j < t->n_cols; j++) free(t->cells[i][j]);
        free(t->cells[i]);
    }
    free(t->cells);
    memset(t, 0, sizeof(*t));
}

static int tsv_col_idx(const Tsv* t, const char* name) {
    for (int j = 0; j < t->n_cols; j++) if (!strcmp(t->col_names[j], name)) return j;
    return -1;
}

static char** split_tab(const char* line, int* n_out) {
    int cap = 16, n = 0;
    char** out = (char**)malloc((size_t)cap * sizeof(char*));
    const char* p = line;
    while (1) {
        const char* tab = strchr(p, '\t');
        size_t len = tab ? (size_t)(tab - p) : strlen(p);
        char* cell = (char*)malloc(len + 1);
        memcpy(cell, p, len); cell[len] = 0;
        if (n == cap) { cap *= 2; out = (char**)realloc(out, (size_t)cap * sizeof(char*)); }
        out[n++] = cell;
        if (!tab) break;
        p = tab + 1;
    }
    *n_out = n;
    return out;
}

static int tsv_load(const char* path, Tsv* t) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[region_test] Cannot open %s\n", path); return 0; }
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
            fprintf(stderr, "[region_test] Skipping row with %d cols (expected %d)\n", nf, t->n_cols);
            for (int j = 0; j < nf; j++) free(parts[j]);
            free(parts);
            continue;
        }
        if (t->n_rows == t->cap) {
            t->cap = t->cap ? t->cap * 2 : 1024;
            t->cells = (char***)realloc(t->cells, (size_t)t->cap * sizeof(char**));
        }
        t->cells[t->n_rows++] = parts;
    }
    free(line);
    fclose(f);
    return t->n_rows;
}

static int parse_numeric(const char* s, double* out) {
    if (!s || s[0] == 0) return 0;
    if (!strcasecmp(s, "NA") || !strcasecmp(s, "nan") || !strcasecmp(s, "na")) return 0;
    if (!strcasecmp(s, "inf") || !strcasecmp(s, "infinity")) { *out = INFINITY; return 0; }
    if (!strcasecmp(s, "-inf"))                              { *out = -INFINITY; return 0; }
    char* end;
    double v = strtod(s, &end);
    if (end == s) return 0;
    if (!isfinite(v)) return 0;
    *out = v;
    return 1;
}

// ── Comma-split helper ───────────────────────────────────────────────────────

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

// ── Wilcoxon rank-sum (tie-corrected normal approximation) ───────────────────

typedef struct { double v; int group; } VG;

static int vg_cmp(const void* a, const void* b) {
    const VG* va = (const VG*)a, *vb = (const VG*)b;
    if (va->v < vb->v) return -1;
    if (va->v > vb->v) return  1;
    return 0;
}

// Returns W (rank-sum for group 0), z, p (two-sided). Inputs are interleaved
// values + binary group labels (0 or 1).
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
        *z_out = NAN;
        *p_out = NAN;
    }
}

// ── Permutation test: shuffle group labels, recompute |U - mu_U|, count >= obs.

static void perm_test(const VG* arr, int n, int n0, int n1,
                      double U_obs, int n_perm, unsigned int* rng,
                      double* p_out) {
    if (n_perm <= 0 || n0 <= 0 || n1 <= 0) { *p_out = NAN; return; }
    // Sort once; we'll permute group labels without re-sorting (ranks come from sorted index).
    VG* sorted = (VG*)malloc((size_t)n * sizeof(VG));
    memcpy(sorted, arr, (size_t)n * sizeof(VG));
    qsort(sorted, (size_t)n, sizeof(VG), vg_cmp);
    // Precompute rank with tie averaging.
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

    // Index pool for permutations: 0..n-1 → first n0 are group 0.
    int* perm_labels = (int*)malloc((size_t)n * sizeof(int));
    for (int j = 0; j < n; j++) perm_labels[j] = (j < n0) ? 0 : 1;

    int ge = 0;
    for (int b = 0; b < n_perm; b++) {
        // Fisher-Yates on perm_labels.
        for (int j = n - 1; j > 0; j--) {
            int r = (int)(rand_r(rng) % (unsigned int)(j + 1));
            int tmp = perm_labels[j]; perm_labels[j] = perm_labels[r]; perm_labels[r] = tmp;
        }
        // Compute U for group 0 under this permutation.
        double sr = 0;
        for (int j = 0; j < n; j++) if (perm_labels[j] == 0) sr += rank[j];
        double U = sr - (double)n0 * (n0 + 1) / 2.0;
        if (fabs(U - mu_U) >= abs_obs) ge++;
    }
    *p_out = (double)(ge + 1) / (double)(n_perm + 1);  // +1 smoothing (standard)

    free(perm_labels); free(rank); free(sorted);
}

static double mean_of(const double* v, int n) {
    if (n == 0) return NAN;
    double s = 0;
    for (int i = 0; i < n; i++) s += v[i];
    return s / n;
}

static int dbl_cmp(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y ? 1 : 0);
}

static double median_of(double* v, int n) {
    if (n == 0) return NAN;
    qsort(v, (size_t)n, sizeof(double), dbl_cmp);
    return (n % 2) ? v[n/2] : 0.5 * (v[n/2 - 1] + v[n/2]);
}

// ── Main ──────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "region_test — Wilcoxon + permutation: per-region vs rest, for any TSV value cols.\n"
        "\n"
        "  --input <tsv>            Required. TSV with header.\n"
        "  --region_col <name>      Partition column (default: region).\n"
        "  --value_cols a,b,c       Numeric cols to test (required, comma-separated).\n"
        "  --rest_label <name>      Optional. Value meaning 'rest of genome' explicitly;\n"
        "                            otherwise rest = all rows not in the region being tested.\n"
        "  --skip_regions a,b       Optional. Region labels to ignore entirely.\n"
        "  --permutations N         Perm reps (default 1000; 0 = skip perm test).\n"
        "  --seed K                 RNG seed (default time(NULL)).\n"
        "  --ncores N               OpenMP across (region × value_col) tests.\n"
        "  --out <f>                Output TSV (default stdout).\n");
}

typedef struct {
    char* region;
    char* value_col;
    int   n_in, n_out;
    double mean_in, mean_out, median_in, median_out;
    double W, z, wilcoxon_p, perm_p;
} TestResult;

int main(int argc, char** argv) {
    const char *in_path = NULL, *out_path = NULL;
    const char *region_col = "region";
    const char *value_cols_spec = NULL;
    const char *rest_label = NULL;
    const char *skip_regions_spec = NULL;
    int n_perm = 1000, ncores = 1;
    unsigned int seed = (unsigned int)time(NULL);
    int seed_set = 0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--input")        && i+1<argc) in_path = argv[++i];
        else if (!strcmp(argv[i], "--region_col")   && i+1<argc) region_col = argv[++i];
        else if (!strcmp(argv[i], "--value_cols")   && i+1<argc) value_cols_spec = argv[++i];
        else if (!strcmp(argv[i], "--rest_label")   && i+1<argc) rest_label = argv[++i];
        else if (!strcmp(argv[i], "--skip_regions") && i+1<argc) skip_regions_spec = argv[++i];
        else if (!strcmp(argv[i], "--permutations") && i+1<argc) n_perm = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed")         && i+1<argc) { seed = (unsigned int)atoi(argv[++i]); seed_set = 1; }
        else if (!strcmp(argv[i], "--ncores")       && i+1<argc) ncores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out")          && i+1<argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[region_test] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }

    if (!in_path || !value_cols_spec) { print_usage(); return 1; }

    #ifdef _OPENMP
    omp_set_num_threads(ncores > 0 ? ncores : 1);
    #else
    (void)ncores;
    #endif

    Tsv t;
    if (!tsv_load(in_path, &t)) { fprintf(stderr, "[region_test] empty input\n"); return 2; }
    fprintf(stderr, "[region_test] Loaded %d rows × %d cols\n", t.n_rows, t.n_cols);

    int reg_idx = tsv_col_idx(&t, region_col);
    if (reg_idx < 0) { fprintf(stderr, "[region_test] Region col '%s' not found\n", region_col); tsv_free(&t); return 3; }

    int nv;
    char** value_cols = split_csv(value_cols_spec, &nv);
    int* value_idx = (int*)malloc((size_t)nv * sizeof(int));
    for (int i = 0; i < nv; i++) {
        value_idx[i] = tsv_col_idx(&t, value_cols[i]);
        if (value_idx[i] < 0) {
            fprintf(stderr, "[region_test] Value col '%s' not found\n", value_cols[i]);
            tsv_free(&t); return 4;
        }
    }

    int n_skip = 0;
    char** skip = NULL;
    if (skip_regions_spec) skip = split_csv(skip_regions_spec, &n_skip);

    // Collect unique region labels.
    char** regions = NULL;
    int n_regions = 0, cap_regions = 0;
    for (int i = 0; i < t.n_rows; i++) {
        const char* r = t.cells[i][reg_idx];
        if (rest_label && !strcmp(r, rest_label)) continue;
        int skip_it = 0;
        for (int k = 0; k < n_skip; k++) if (!strcmp(r, skip[k])) { skip_it = 1; break; }
        if (skip_it) continue;
        int seen = 0;
        for (int k = 0; k < n_regions; k++) if (!strcmp(regions[k], r)) { seen = 1; break; }
        if (!seen) {
            if (n_regions == cap_regions) {
                cap_regions = cap_regions ? cap_regions * 2 : 64;
                regions = (char**)realloc(regions, (size_t)cap_regions * sizeof(char*));
            }
            regions[n_regions++] = strdup(r);
        }
    }
    fprintf(stderr, "[region_test] %d region labels × %d value cols = %d tests\n",
            n_regions, nv, n_regions * nv);
    if (n_perm > 0) fprintf(stderr, "[region_test] permutations=%d seed=%u%s\n",
                             n_perm, seed, seed_set ? "" : " (time-based)");

    int n_tests = n_regions * nv;
    TestResult* results = (TestResult*)calloc((size_t)n_tests, sizeof(TestResult));

    #pragma omp parallel for schedule(dynamic)
    for (int t_idx = 0; t_idx < n_tests; t_idx++) {
        int r_idx = t_idx / nv;
        int v_idx = t_idx % nv;
        const char* region = regions[r_idx];
        int vi = value_idx[v_idx];
        unsigned int local_seed = seed ^ (unsigned int)(t_idx * 2654435761u);

        VG* arr = (VG*)malloc((size_t)t.n_rows * sizeof(VG));
        double* vals_in  = (double*)malloc((size_t)t.n_rows * sizeof(double));
        double* vals_out = (double*)malloc((size_t)t.n_rows * sizeof(double));
        int n = 0, n_in = 0, n_out = 0;

        for (int i = 0; i < t.n_rows; i++) {
            const char* rg = t.cells[i][reg_idx];
            int in_test = (!strcmp(rg, region)) ? 1 : 0;
            if (!in_test) {
                if (rest_label && strcmp(rg, rest_label) != 0) continue;
                int skip_it = 0;
                for (int k = 0; k < n_skip; k++) if (!strcmp(rg, skip[k])) { skip_it = 1; break; }
                if (skip_it) continue;
            }
            double v;
            if (!parse_numeric(t.cells[i][vi], &v)) continue;
            arr[n].v = v;
            arr[n].group = in_test ? 0 : 1;
            n++;
            if (in_test) vals_in[n_in++] = v;
            else         vals_out[n_out++] = v;
        }

        TestResult* R = &results[t_idx];
        R->region = strdup(region);
        R->value_col = strdup(value_cols[v_idx]);
        R->n_in = n_in; R->n_out = n_out;
        R->mean_in = mean_of(vals_in, n_in);
        R->mean_out = mean_of(vals_out, n_out);
        R->median_in = median_of(vals_in, n_in);
        R->median_out = median_of(vals_out, n_out);
        R->W = R->z = R->wilcoxon_p = R->perm_p = NAN;

        if (n_in >= 2 && n_out >= 2) {
            wilcoxon(arr, n, n_in, n_out, &R->W, &R->z, &R->wilcoxon_p);
            if (n_perm > 0) perm_test(arr, n, n_in, n_out, R->W, n_perm, &local_seed, &R->perm_p);
        }

        free(arr); free(vals_in); free(vals_out);
    }

    FILE* fout = out_path ? fopen(out_path, "w") : stdout;
    if (!fout) { fprintf(stderr, "[region_test] Cannot open --out %s\n", out_path); return 5; }

    fprintf(fout, "# schema_version=%s permutations=%d\n", SCHEMA_VERSION, n_perm);
    fprintf(fout, "region\tvalue_col\tn_in\tn_out\tmean_in\tmean_out\tmedian_in\tmedian_out\t"
                  "effect_diff\tW\tz\twilcoxon_p\tperm_p\n");

    for (int i = 0; i < n_tests; i++) {
        TestResult* R = &results[i];
        fprintf(fout, "%s\t%s\t%d\t%d", R->region, R->value_col, R->n_in, R->n_out);
        #define EM(v) do { if (isnan(v)) fprintf(fout, "\tNA"); else fprintf(fout, "\t%.6g", v); } while (0)
        EM(R->mean_in); EM(R->mean_out); EM(R->median_in); EM(R->median_out);
        EM(R->mean_in - R->mean_out);
        EM(R->W); EM(R->z); EM(R->wilcoxon_p); EM(R->perm_p);
        #undef EM
        fputc('\n', fout);
        free(R->region); free(R->value_col);
    }
    if (out_path) fclose(fout);

    for (int i = 0; i < n_regions; i++) free(regions[i]);
    free(regions); free(results); free(value_idx);
    for (int i = 0; i < nv; i++) free(value_cols[i]);
    free(value_cols);
    if (skip) { for (int i = 0; i < n_skip; i++) free(skip[i]); free(skip); }
    tsv_free(&t);
    return 0;
}
