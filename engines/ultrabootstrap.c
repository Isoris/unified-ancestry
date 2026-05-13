// =============================================================================
// ultrabootstrap.c — generic bootstrap CI for arbitrary TSV value columns.
//
// Modes:
//   row bootstrap (default): resample TSV rows with replacement.
//   block bootstrap (--block_col):  resample whole groups of rows that share
//                                    a block label as a unit (keeps within-block
//                                    correlation intact — e.g. codons within a
//                                    locus, sites within a chromosome).
//
// Statistics: mean / median / sum / sd. Multiple --statistic mean,median to get
// several at once.
//
// Optional grouping (--group_col): emit one CI per (group, column, statistic).
//
// For multi-FASTA inputs: pre-process the alignment into a TSV (one row per
// codon or site with the values you care about; --block_col locus_id) and
// then run ultrabootstrap. The C side stays purely numeric.
//
// Compile: gcc -O3 -march=native -fopenmp -o ultrabootstrap ultrabootstrap.c -lm
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

#define SCHEMA_VERSION "ultrabootstrap_v1"

// ── TSV loader (lightweight) ─────────────────────────────────────────────────

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

static int tsv_col_idx(const Tsv* t, const char* name) {
    for (int j = 0; j < t->n_cols; j++) if (!strcmp(t->col_names[j], name)) return j;
    return -1;
}

static int tsv_load(const char* path, Tsv* t) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[ultraboot] Cannot open %s\n", path); return 0; }
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
            t->cap = t->cap ? t->cap * 2 : 1024;
            t->cells = (char***)realloc(t->cells, (size_t)t->cap * sizeof(char**));
        }
        t->cells[t->n_rows++] = parts;
    }
    free(line); fclose(f);
    return t->n_rows;
}

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

static int parse_numeric(const char* s, double* out) {
    if (!s || s[0] == 0) return 0;
    if (!strcasecmp(s, "NA") || !strcasecmp(s, "nan")) return 0;
    char* end;
    double v = strtod(s, &end);
    if (end == s || !isfinite(v)) return 0;
    *out = v;
    return 1;
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

// ── Statistics ───────────────────────────────────────────────────────────────

typedef enum { STAT_MEAN, STAT_MEDIAN, STAT_SUM, STAT_SD } StatKind;

static int parse_stat(const char* s, StatKind* k) {
    if      (!strcasecmp(s, "mean"))   { *k = STAT_MEAN;   return 1; }
    else if (!strcasecmp(s, "median")) { *k = STAT_MEDIAN; return 1; }
    else if (!strcasecmp(s, "sum"))    { *k = STAT_SUM;    return 1; }
    else if (!strcasecmp(s, "sd") || !strcasecmp(s, "stdev")) { *k = STAT_SD; return 1; }
    return 0;
}

static const char* stat_name(StatKind k) {
    return k == STAT_MEAN ? "mean" : k == STAT_MEDIAN ? "median" :
           k == STAT_SUM ? "sum" : "sd";
}

static int dbl_cmp(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y ? 1 : 0);
}

static double compute_stat(double* arr, int n, StatKind k) {
    if (n == 0) return NAN;
    if (k == STAT_SUM) {
        double s = 0;
        for (int i = 0; i < n; i++) s += arr[i];
        return s;
    }
    if (k == STAT_MEAN) {
        double s = 0;
        for (int i = 0; i < n; i++) s += arr[i];
        return s / n;
    }
    if (k == STAT_SD) {
        if (n < 2) return NAN;
        double m = 0;
        for (int i = 0; i < n; i++) m += arr[i];
        m /= n;
        double v = 0;
        for (int i = 0; i < n; i++) { double d = arr[i] - m; v += d * d; }
        return sqrt(v / (n - 1));
    }
    // median: needs sorted copy
    qsort(arr, (size_t)n, sizeof(double), dbl_cmp);
    return (n % 2) ? arr[n/2] : 0.5 * (arr[n/2 - 1] + arr[n/2]);
}

static double percentile(double* arr, int n, double p) {
    int m = 0;
    for (int i = 0; i < n; i++) if (isfinite(arr[i])) arr[m++] = arr[i];
    if (m == 0) return NAN;
    qsort(arr, (size_t)m, sizeof(double), dbl_cmp);
    int idx = (int)floor(p * (m - 1) + 0.5);
    if (idx < 0) idx = 0;
    if (idx >= m) idx = m - 1;
    return arr[idx];
}

// ── Bootstrap (row or block) ─────────────────────────────────────────────────

// Block index: each row belongs to a block (0..n_blocks-1); rows[] lists
// the row indices in each block.
typedef struct {
    int* row_start;   // row_start[b] = first row index in block b (in compressed list)
    int* row_list;    // packed list of row indices, length = total rows
    int  n_blocks;
} BlockIdx;

static void build_block_idx(const Tsv* t, int block_col, BlockIdx* out) {
    // Collect unique block labels.
    char** labels = NULL;
    int n_lab = 0, cap_lab = 0;
    int* block_of_row = (int*)malloc((size_t)t->n_rows * sizeof(int));
    for (int i = 0; i < t->n_rows; i++) {
        const char* lab = t->cells[i][block_col];
        int found = -1;
        for (int k = 0; k < n_lab; k++) if (!strcmp(labels[k], lab)) { found = k; break; }
        if (found < 0) {
            if (n_lab == cap_lab) {
                cap_lab = cap_lab ? cap_lab * 2 : 64;
                labels = (char**)realloc(labels, (size_t)cap_lab * sizeof(char*));
            }
            labels[n_lab] = strdup(lab);
            found = n_lab++;
        }
        block_of_row[i] = found;
    }
    // Build counts per block, then offsets.
    int* counts = (int*)calloc((size_t)n_lab, sizeof(int));
    for (int i = 0; i < t->n_rows; i++) counts[block_of_row[i]]++;
    out->row_start = (int*)malloc((size_t)(n_lab + 1) * sizeof(int));
    out->row_start[0] = 0;
    for (int b = 0; b < n_lab; b++) out->row_start[b + 1] = out->row_start[b] + counts[b];
    out->row_list = (int*)malloc((size_t)t->n_rows * sizeof(int));
    int* fill = (int*)calloc((size_t)n_lab, sizeof(int));
    for (int i = 0; i < t->n_rows; i++) {
        int b = block_of_row[i];
        out->row_list[out->row_start[b] + fill[b]++] = i;
    }
    out->n_blocks = n_lab;
    for (int k = 0; k < n_lab; k++) free(labels[k]);
    free(labels);
    free(counts); free(fill); free(block_of_row);
}

static void free_block_idx(BlockIdx* b) {
    free(b->row_start); free(b->row_list); memset(b, 0, sizeof(*b));
}

// ── Main ──────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "ultrabootstrap — generic resample CI for arbitrary TSV value columns.\n"
        "\n"
        "  --input <tsv>            Required. TSV with header.\n"
        "  --value_cols a,b,c       Required. Numeric columns.\n"
        "  --statistic mean[,median,sum,sd]\n"
        "                           Statistics to bootstrap (default mean).\n"
        "  --group_col <name>       Optional. Emit one CI per (group, col, stat).\n"
        "  --block_col <name>       Optional. Block bootstrap: resample blocks of\n"
        "                            rows with the same block label as a unit\n"
        "                            (keeps within-block correlation). Without\n"
        "                            --block_col, rows are resampled individually.\n"
        "  --n_boot N               Bootstrap reps (default 1000).\n"
        "  --ci_lo F                Lower percentile (default 0.025).\n"
        "  --ci_hi F                Upper percentile (default 0.975).\n"
        "  --seed K                 RNG seed (default time(NULL)).\n"
        "  --ncores N               OpenMP across (group × col × stat) tasks.\n"
        "  --out <f>                Output TSV (default stdout).\n"
        "\n"
        "Output (schema_version=" SCHEMA_VERSION "):\n"
        "  group value_col statistic n observed boot_mean boot_sd boot_lo boot_hi\n"
        "\n"
        "Tip — multi-FASTA: pre-process the alignment into a TSV with one row\n"
        "per codon (or site) carrying the numeric contributions you want to\n"
        "bootstrap, with locus_id (or chrom) as --block_col. ultrabootstrap\n"
        "stays purely numeric.\n");
}

typedef struct {
    char* group;
    char* col;
    StatKind stat;
    int    n;
    double observed;
    double boot_mean;
    double boot_sd;
    double boot_lo;
    double boot_hi;
} BootResult;

int main(int argc, char** argv) {
    const char *in_path = NULL, *out_path = NULL;
    const char *value_cols_spec = NULL, *stat_spec = "mean";
    const char *group_col_name = NULL, *block_col_name = NULL;
    int n_boot = 1000, ncores = 1;
    double ci_lo = 0.025, ci_hi = 0.975;
    unsigned int seed = (unsigned int)time(NULL);
    int seed_set = 0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--input")      && i+1<argc) in_path = argv[++i];
        else if (!strcmp(argv[i], "--value_cols") && i+1<argc) value_cols_spec = argv[++i];
        else if (!strcmp(argv[i], "--statistic")  && i+1<argc) stat_spec = argv[++i];
        else if (!strcmp(argv[i], "--group_col")  && i+1<argc) group_col_name = argv[++i];
        else if (!strcmp(argv[i], "--block_col")  && i+1<argc) block_col_name = argv[++i];
        else if (!strcmp(argv[i], "--n_boot")     && i+1<argc) n_boot = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--ci_lo")      && i+1<argc) ci_lo = atof(argv[++i]);
        else if (!strcmp(argv[i], "--ci_hi")      && i+1<argc) ci_hi = atof(argv[++i]);
        else if (!strcmp(argv[i], "--seed")       && i+1<argc) { seed = (unsigned int)atoi(argv[++i]); seed_set = 1; }
        else if (!strcmp(argv[i], "--ncores")     && i+1<argc) ncores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out")        && i+1<argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[ultraboot] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }

    if (!in_path || !value_cols_spec) { print_usage(); return 1; }

    #ifdef _OPENMP
    omp_set_num_threads(ncores > 0 ? ncores : 1);
    #else
    (void)ncores;
    #endif

    Tsv t;
    if (!tsv_load(in_path, &t)) { fprintf(stderr, "[ultraboot] empty input\n"); return 2; }
    fprintf(stderr, "[ultraboot] Loaded %d rows × %d cols\n", t.n_rows, t.n_cols);

    int nv, nstat;
    char** value_cols = split_csv(value_cols_spec, &nv);
    char** stat_strs  = split_csv(stat_spec, &nstat);
    int* value_idx = (int*)malloc((size_t)nv * sizeof(int));
    for (int i = 0; i < nv; i++) {
        value_idx[i] = tsv_col_idx(&t, value_cols[i]);
        if (value_idx[i] < 0) { fprintf(stderr, "[ultraboot] col '%s' not found\n", value_cols[i]); return 3; }
    }
    StatKind* stats = (StatKind*)malloc((size_t)nstat * sizeof(StatKind));
    for (int i = 0; i < nstat; i++) {
        if (!parse_stat(stat_strs[i], &stats[i])) {
            fprintf(stderr, "[ultraboot] unknown statistic '%s'\n", stat_strs[i]); return 4;
        }
    }

    int group_col = -1;
    if (group_col_name) {
        group_col = tsv_col_idx(&t, group_col_name);
        if (group_col < 0) { fprintf(stderr, "[ultraboot] --group_col '%s' not found\n", group_col_name); return 5; }
    }
    int block_col = -1;
    if (block_col_name) {
        block_col = tsv_col_idx(&t, block_col_name);
        if (block_col < 0) { fprintf(stderr, "[ultraboot] --block_col '%s' not found\n", block_col_name); return 6; }
    }

    // Collect unique groups (one task per group).
    char** groups = NULL;
    int n_groups = 0, cap_g = 0;
    if (group_col < 0) {
        groups = (char**)malloc(sizeof(char*));
        groups[0] = strdup("ALL");
        n_groups = 1;
    } else {
        for (int i = 0; i < t.n_rows; i++) {
            const char* g = t.cells[i][group_col];
            int seen = 0;
            for (int k = 0; k < n_groups; k++) if (!strcmp(groups[k], g)) { seen = 1; break; }
            if (!seen) {
                if (n_groups == cap_g) {
                    cap_g = cap_g ? cap_g * 2 : 32;
                    groups = (char**)realloc(groups, (size_t)cap_g * sizeof(char*));
                }
                groups[n_groups++] = strdup(g);
            }
        }
    }
    fprintf(stderr, "[ultraboot] %d groups × %d cols × %d stats = %d tasks; "
                    "n_boot=%d seed=%u%s%s\n",
            n_groups, nv, nstat, n_groups * nv * nstat, n_boot, seed,
            seed_set ? "" : " (time-based)",
            block_col >= 0 ? ", block bootstrap" : ", row bootstrap");

    BlockIdx blk;
    if (block_col >= 0) {
        build_block_idx(&t, block_col, &blk);
        fprintf(stderr, "[ultraboot] %d blocks\n", blk.n_blocks);
    }

    int n_tasks = n_groups * nv * nstat;
    BootResult* results = (BootResult*)calloc((size_t)n_tasks, sizeof(BootResult));

    #pragma omp parallel for schedule(dynamic)
    for (int task = 0; task < n_tasks; task++) {
        int g_idx = task / (nv * nstat);
        int rem = task % (nv * nstat);
        int v_idx = rem / nstat;
        int s_idx = rem % nstat;

        const char* group = groups[g_idx];
        int vi = value_idx[v_idx];
        StatKind sk = stats[s_idx];
        unsigned int local_seed = seed ^ (unsigned int)(task * 2654435761u);

        // Collect numeric values for this group (and corresponding row→block map if block mode).
        double* vals = (double*)malloc((size_t)t.n_rows * sizeof(double));
        int n_vals = 0;
        // Per-block list of row indices (only used in block mode)
        int* group_blocks = NULL;
        int* group_block_row_start = NULL;
        int* group_block_row_list = NULL;
        int n_group_blocks = 0;
        if (block_col < 0) {
            for (int i = 0; i < t.n_rows; i++) {
                if (group_col >= 0 && strcmp(t.cells[i][group_col], group) != 0) continue;
                double v;
                if (parse_numeric(t.cells[i][vi], &v)) vals[n_vals++] = v;
            }
        } else {
            // Restrict to blocks that have at least one row in this group; within those
            // blocks, collect all rows (block bootstrap resamples whole blocks).
            int* tmp_blocks = (int*)malloc((size_t)blk.n_blocks * sizeof(int));
            for (int b = 0; b < blk.n_blocks; b++) {
                int in_group = 0;
                for (int k = blk.row_start[b]; k < blk.row_start[b + 1]; k++) {
                    int row = blk.row_list[k];
                    if (group_col < 0 || !strcmp(t.cells[row][group_col], group)) { in_group = 1; break; }
                }
                if (in_group) tmp_blocks[n_group_blocks++] = b;
            }
            group_blocks = tmp_blocks;
            group_block_row_start = (int*)malloc((size_t)(n_group_blocks + 1) * sizeof(int));
            group_block_row_start[0] = 0;
            int total = 0;
            for (int gb = 0; gb < n_group_blocks; gb++) {
                int b = group_blocks[gb];
                int cnt = 0;
                for (int k = blk.row_start[b]; k < blk.row_start[b + 1]; k++) {
                    int row = blk.row_list[k];
                    if (group_col >= 0 && strcmp(t.cells[row][group_col], group) != 0) continue;
                    cnt++;
                }
                total += cnt;
                group_block_row_start[gb + 1] = total;
            }
            group_block_row_list = (int*)malloc((size_t)total * sizeof(int));
            int fill = 0;
            for (int gb = 0; gb < n_group_blocks; gb++) {
                int b = group_blocks[gb];
                for (int k = blk.row_start[b]; k < blk.row_start[b + 1]; k++) {
                    int row = blk.row_list[k];
                    if (group_col >= 0 && strcmp(t.cells[row][group_col], group) != 0) continue;
                    group_block_row_list[fill++] = row;
                }
            }
            for (int i = 0; i < total; i++) {
                double v;
                if (parse_numeric(t.cells[group_block_row_list[i]][vi], &v)) vals[n_vals++] = v;
            }
        }

        BootResult* R = &results[task];
        R->group = strdup(group);
        R->col = strdup(value_cols[v_idx]);
        R->stat = sk;
        R->n = n_vals;
        R->observed = R->boot_mean = R->boot_sd = R->boot_lo = R->boot_hi = NAN;

        if (n_vals < 2 || n_boot <= 0) {
            if (n_vals > 0) {
                double* tmp = (double*)malloc((size_t)n_vals * sizeof(double));
                memcpy(tmp, vals, (size_t)n_vals * sizeof(double));
                R->observed = compute_stat(tmp, n_vals, sk);
                free(tmp);
            }
            free(vals);
            free(group_blocks); free(group_block_row_start); free(group_block_row_list);
            continue;
        }

        // Observed
        {
            double* tmp = (double*)malloc((size_t)n_vals * sizeof(double));
            memcpy(tmp, vals, (size_t)n_vals * sizeof(double));
            R->observed = compute_stat(tmp, n_vals, sk);
            free(tmp);
        }

        // Bootstrap
        double* reps = (double*)malloc((size_t)n_boot * sizeof(double));
        double* tmp_vals = (double*)malloc((size_t)n_vals * sizeof(double));

        for (int b = 0; b < n_boot; b++) {
            int rep_n = 0;
            if (block_col < 0) {
                for (int k = 0; k < n_vals; k++) {
                    int idx = (int)(rand_r(&local_seed) % (unsigned int)n_vals);
                    tmp_vals[rep_n++] = vals[idx];
                }
            } else {
                // Resample blocks with replacement (each block contributes its rows).
                for (int gb = 0; gb < n_group_blocks; gb++) {
                    int pick = (int)(rand_r(&local_seed) % (unsigned int)n_group_blocks);
                    int start = group_block_row_start[pick];
                    int end   = group_block_row_start[pick + 1];
                    for (int r = start; r < end; r++) {
                        double v;
                        if (parse_numeric(t.cells[group_block_row_list[r]][vi], &v))
                            tmp_vals[rep_n++] = v;
                        if (rep_n >= n_vals) break;
                    }
                    if (rep_n >= n_vals) break;
                }
            }
            reps[b] = compute_stat(tmp_vals, rep_n, sk);
        }

        // Summary
        double s = 0; int eff = 0;
        for (int b = 0; b < n_boot; b++) if (isfinite(reps[b])) { s += reps[b]; eff++; }
        R->boot_mean = (eff > 0) ? s / eff : NAN;
        double v = 0;
        for (int b = 0; b < n_boot; b++) if (isfinite(reps[b])) {
            double d = reps[b] - R->boot_mean; v += d * d;
        }
        R->boot_sd = (eff > 1) ? sqrt(v / (eff - 1)) : NAN;
        R->boot_lo = percentile(reps, n_boot, ci_lo);
        R->boot_hi = percentile(reps, n_boot, ci_hi);

        free(reps); free(tmp_vals); free(vals);
        free(group_blocks); free(group_block_row_start); free(group_block_row_list);
    }

    FILE* fout = out_path ? fopen(out_path, "w") : stdout;
    if (!fout) { fprintf(stderr, "[ultraboot] Cannot open --out %s\n", out_path); return 7; }
    fprintf(fout, "# schema_version=%s n_boot=%d ci=[%.4g,%.4g] mode=%s\n",
            SCHEMA_VERSION, n_boot, ci_lo, ci_hi, block_col >= 0 ? "block" : "row");
    fprintf(fout, "group\tvalue_col\tstatistic\tn\tobserved\tboot_mean\tboot_sd\tboot_lo\tboot_hi\n");

    for (int i = 0; i < n_tasks; i++) {
        BootResult* R = &results[i];
        fprintf(fout, "%s\t%s\t%s\t%d", R->group, R->col, stat_name(R->stat), R->n);
        #define EM(v) do { if (isnan(v)) fprintf(fout, "\tNA"); else fprintf(fout, "\t%.6g", v); } while (0)
        EM(R->observed); EM(R->boot_mean); EM(R->boot_sd); EM(R->boot_lo); EM(R->boot_hi);
        #undef EM
        fputc('\n', fout);
        free(R->group); free(R->col);
    }
    if (out_path) fclose(fout);

    if (block_col >= 0) free_block_idx(&blk);
    for (int i = 0; i < n_groups; i++) free(groups[i]);
    free(groups);
    for (int i = 0; i < nv; i++) free(value_cols[i]);
    free(value_cols);
    for (int i = 0; i < nstat; i++) free(stat_strs[i]);
    free(stat_strs);
    free(value_idx); free(stats); free(results);
    tsv_free(&t);
    return 0;
}
