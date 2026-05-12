// =============================================================================
// xpehh.c — per-window XP-EHH track from BEAGLE genotype likelihoods.
//
// ── What this binary computes ───────────────────────────────────────────────
//
// XP-EHH (Sabeti et al. 2007) is a haplotype-based selection scan: at each
// focal SNP it asks whether the test cohort carries longer extended haplotype
// homozygosity (EHH) than a reference cohort. Positive XP-EHH = test cohort
// has longer haplotypes (consistent with a recent sweep / inversion-suppressed
// region in the test cohort).
//
// Sabeti-strict XP-EHH needs PHASED haplotypes. This binary computes a
// defensible approximation directly from BEAGLE GLs without phasing:
//
//   For each within-cohort pair (i, j) and each site k:
//     IBS_prob_ij(k) = gl0_i*gl0_j + gl1_i*gl1_j + gl2_i*gl2_j   ∈ [0, 1]
//                   = P(samples i and j carry the same diploid genotype at k)
//
//   EHH(focal, d) = mean over pairs of  ∏_{k along walk} IBS_prob_ij(k)
//                = "diploid IBS-based EHH"
//
//   iHH(focal) = ∫ EHH(focal, d) dd        (trapezoidal, in bp)
//                walked outward from focal in both directions, stopping at
//                EHH < --ehh_decay or distance > --max_distance_bp.
//
//   XP-EHH(focal) = ln( iHH_test(focal) / iHH_ref(focal) )
//
// Caveats (documented per the spec):
//   - This is NOT Sabeti-strict XP-EHH. True XP-EHH requires phasing; this
//     binary uses diploid genotype-likelihood IBS instead, which is more
//     conservative (heterozygous-heterozygous pairs hit IBS_prob ≈ 0.5
//     even when haplotypes match).
//   - Output sign convention matches Sabeti: positive = test cohort has
//     longer EHH than reference.
//   - Norm step (z-score against autosomal background) follows Sabeti
//     recommendation. Whole-input mean/sd used.
//
// Per-window aggregation: each focal SNP belongs to the window containing
// its position. Aggregates: mean / max-|·| of XP-EHH, mean of normalized.
//
// ── Inputs ──────────────────────────────────────────────────────────────────
//   --beagle <f.beagle.gz>     BEAGLE GL file (same format as other engines)
//   --sample_list <f>          Sample IDs (same order as BEAGLE)
//   --test_samples <f>         One sample ID per line — test cohort
//   --ref_samples <f>          One sample ID per line — reference cohort
//   --chr <name>               Filter to a single chromosome.
//   --windows <bed>            BED of windows for aggregation.
//   --fixed_win W:S            Fixed-bp windows of W bp, step S.
//   --ehh_decay F              Stop walk when EHH < F (default 0.05).
//   --max_distance_bp N        Maximum walk distance per arm (default 1e6).
//   --min_n_test N             Skip if test cohort has < N samples (default 20).
//   --min_n_ref N              Same for ref (default 20).
//   --focal_step K             Compute XP-EHH at every Kth SNP (default 1).
//   --out <f>                  Output TSV (default stdout).
//   --ncores N                 OpenMP across focal SNPs.
//
// ── Outputs ─────────────────────────────────────────────────────────────────
//   Per-window TSV (schema_version=xpehh_v1):
//     chrom window_start window_end n_focal n_focal_with_xpehh
//     iHH_test_mean iHH_ref_mean
//     xpehh_mean xpehh_max_abs norm_xpehh_mean norm_xpehh_max_abs
//
//   Header line carries test_cohort_size / ref_cohort_size / global_mean /
//   global_sd for transparency.
//
// Compile: gcc -O3 -march=native -fopenmp -o xpehh xpehh.c -lz -lm
// =============================================================================

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>
#include <zlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_IND 512
#define SCHEMA_VERSION "xpehh_v1"

// ── Gz reader (same as other engines) ───────────────────────────────────────

typedef struct {
    gzFile fp;
    char buf[1<<16];
    char lo[1<<18];
    int lo_len;
} Gz;

static int gz_open(Gz* g, const char* path) {
    g->fp = gzopen(path, "rb");
    if (!g->fp) return 0;
    gzbuffer(g->fp, 1<<18);
    g->lo_len = 0;
    return 1;
}

static int gz_getline(Gz* g, char* out, int maxlen) {
    out[0] = 0;
    while (1) {
        for (int i = 0; i < g->lo_len; i++) {
            if (g->lo[i] == '\n') {
                int len = i;
                if (len > 0 && g->lo[len-1] == '\r') len--;
                if (len >= maxlen) len = maxlen - 1;
                memcpy(out, g->lo, len);
                out[len] = 0;
                int rem = g->lo_len - i - 1;
                if (rem > 0) memmove(g->lo, g->lo + i + 1, rem);
                g->lo_len = rem;
                return 1;
            }
        }
        int n = gzread(g->fp, g->buf, sizeof(g->buf) - 1);
        if (n <= 0) {
            if (g->lo_len > 0) {
                int len = g->lo_len;
                if (len >= maxlen) len = maxlen - 1;
                memcpy(out, g->lo, len);
                out[len] = 0;
                g->lo_len = 0;
                return 1;
            }
            return 0;
        }
        if (g->lo_len + n > (int)sizeof(g->lo) - 1) n = sizeof(g->lo) - 1 - g->lo_len;
        memcpy(g->lo + g->lo_len, g->buf, n);
        g->lo_len += n;
    }
}

static void gz_close(Gz* g) { if (g->fp) gzclose(g->fp); }

// ── Sample list / cohort loaders ────────────────────────────────────────────

static int load_samples(const char* path, char names[][64]) {
    FILE* f = fopen(path, "r");
    if (!f) return 0;
    int n = 0;
    char line[512];
    while (fgets(line, sizeof(line), f) && n < MAX_IND) {
        char* p = line;
        while (*p && (*p == ' ' || *p == '\t')) p++;
        int len = strlen(p);
        while (len > 0 && (p[len-1] == '\n' || p[len-1] == '\r' || p[len-1] == ' ')) len--;
        if (len == 0) continue;
        char* slash = strrchr(p, '/');
        if (slash) p = slash + 1;
        char tmp[256];
        int cl = len < 255 ? len : 255;
        memcpy(tmp, p, cl); tmp[cl] = 0;
        char* ext;
        if      ((ext = strstr(tmp, ".sorted.markdup.bam"))) *ext = 0;
        else if ((ext = strstr(tmp, ".markdup.bam")))         *ext = 0;
        else if ((ext = strstr(tmp, ".sorted.bam")))          *ext = 0;
        else if ((ext = strstr(tmp, ".bam")))                  *ext = 0;
        else if ((ext = strstr(tmp, ".cram")))                 *ext = 0;
        int nl = strlen(tmp);
        if (nl > 63) nl = 63;
        memcpy(names[n], tmp, nl); names[n][nl] = 0;
        n++;
    }
    fclose(f);
    return n;
}

static int load_cohort(const char* path, int* indices, char names[][64], int n_all) {
    FILE* f = fopen(path, "r");
    if (!f) return 0;
    int n = 0;
    char line[512];
    while (fgets(line, sizeof(line), f) && n < MAX_IND) {
        int len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) len--;
        line[len] = 0;
        if (len == 0) continue;
        char* q = line;
        char* slash = strrchr(q, '/');
        if (slash) q = slash + 1;
        char tmp[256];
        int cl = strlen(q); if (cl > 255) cl = 255;
        memcpy(tmp, q, cl); tmp[cl] = 0;
        char* ext;
        if ((ext = strstr(tmp, ".sorted.markdup.bam"))) *ext = 0;
        else if ((ext = strstr(tmp, ".bam"))) *ext = 0;
        for (int i = 0; i < n_all; i++) {
            if (!strcmp(names[i], tmp)) { indices[n++] = i; break; }
        }
    }
    fclose(f);
    return n;
}

// ── Per-site GL storage ─────────────────────────────────────────────────────
// 3 floats per (site, sample). Total: 12 * MAX_IND bytes per site.
// At MAX_IND=512, MAX_SITES=5M: ~30 GB. Reduce MAX_SITES at compile time if
// tighter footprint is needed.

typedef struct {
    int pos;
    float gl[MAX_IND][3];   // gl[i][0]=P(AA), [1]=P(AB), [2]=P(BB), normalized.
} Site;

// Grow-on-demand Site array. *cap_out reflects allocated capacity (in sites).
static int load_beagle_gl(const char* path, const char* filter_chr,
                          Site** sites_out, int* cap_out,
                          int n_ind, char* chrom_out, int chrom_out_sz) {
    Gz gz;
    if (!gz_open(&gz, path)) { fprintf(stderr, "[xpehh] Cannot open %s\n", path); return 0; }
    char line[1<<20];
    gz_getline(&gz, line, sizeof(line)); // skip header

    int cap = 65536;
    Site* sites = (Site*)malloc((size_t)cap * sizeof(Site));
    if (!sites) { fprintf(stderr, "[xpehh] initial malloc failed\n"); gz_close(&gz); return 0; }

    int n = 0;
    while (gz_getline(&gz, line, sizeof(line))) {
        if (n == cap) {
            int new_cap = cap * 2;
            Site* tmp = (Site*)realloc(sites, (size_t)new_cap * sizeof(Site));
            if (!tmp) { fprintf(stderr, "[xpehh] realloc to %d sites failed (cap %d)\n", new_cap, cap); break; }
            sites = tmp; cap = new_cap;
        }
        char* tab1 = strchr(line, '\t');
        if (!tab1) continue;
        *tab1 = 0;
        char* us = strrchr(line, '_');
        if (!us) continue;
        *us = 0;
        const char* chr = line;
        int pos = atoi(us + 1);
        if (filter_chr && strcmp(chr, filter_chr) != 0) { *tab1 = '\t'; continue; }
        if (chrom_out && chrom_out[0] == 0)
            strncpy(chrom_out, chr, chrom_out_sz - 1);

        char* p = tab1 + 1;
        char* tab2 = strchr(p, '\t');       if (!tab2) continue;
        char* tab3 = strchr(tab2 + 1, '\t'); if (!tab3) continue;
        p = tab3 + 1;

        sites[n].pos = pos;
        for (int i = 0; i < n_ind; i++) {
            double gl0 = 0, gl1 = 0, gl2 = 0;
            char* endp;
            gl0 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;
            gl1 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;
            gl2 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;
            double s = gl0 + gl1 + gl2;
            if (s > 1e-15) { gl0 /= s; gl1 /= s; gl2 /= s; }
            sites[n].gl[i][0] = (float)gl0;
            sites[n].gl[i][1] = (float)gl1;
            sites[n].gl[i][2] = (float)gl2;
        }
        n++;
    }
    gz_close(&gz);
    *sites_out = sites;
    *cap_out = cap;
    fprintf(stderr, "[xpehh] Loaded %d sites for %s (cap=%d)\n", n,
            filter_chr ? filter_chr : "all", cap);
    return n;
}

// ── Per-cohort iHH walker from a focal SNP ──────────────────────────────────
//
// Walks outward from `focal` (direction = +1 right or -1 left). Maintains a
// running per-pair product of IBS_prob across visited sites. EHH(d) = mean of
// per-pair running products. iHH = trapezoidal integral of EHH over bp.
//
// Stops when EHH < ehh_decay OR walk distance > max_distance_bp OR end-of-array.

static double walk_one_arm(const Site* sites, int n_sites, int focal, int direction,
                            const int* idx, int n,
                            double ehh_decay, int max_distance_bp,
                            int max_pair_buf, double* running_buf) {
    if (n < 2) return 0;
    int n_pairs = n * (n - 1) / 2;
    if (n_pairs > max_pair_buf) return 0; // shouldn't happen if caller sized right

    // Initialize running products at the focal site itself.
    int p = 0;
    for (int a = 0; a < n; a++) {
        int sa = idx[a];
        float ga0 = sites[focal].gl[sa][0];
        float ga1 = sites[focal].gl[sa][1];
        float ga2 = sites[focal].gl[sa][2];
        for (int b = a + 1; b < n; b++) {
            int sb = idx[b];
            float gb0 = sites[focal].gl[sb][0];
            float gb1 = sites[focal].gl[sb][1];
            float gb2 = sites[focal].gl[sb][2];
            running_buf[p++] = (double)(ga0*gb0 + ga1*gb1 + ga2*gb2);
        }
    }

    double ehh_prev = 0;
    for (int q = 0; q < n_pairs; q++) ehh_prev += running_buf[q];
    ehh_prev /= n_pairs;

    int prev_pos = sites[focal].pos;
    int focal_pos = sites[focal].pos;
    double ihh = 0;

    int k = focal + direction;
    while (k >= 0 && k < n_sites) {
        int dist = abs(sites[k].pos - focal_pos);
        if (dist > max_distance_bp) break;
        if (sites[k].pos == prev_pos) { k += direction; continue; }

        // Update running products at site k.
        p = 0;
        for (int a = 0; a < n; a++) {
            int sa = idx[a];
            float ga0 = sites[k].gl[sa][0];
            float ga1 = sites[k].gl[sa][1];
            float ga2 = sites[k].gl[sa][2];
            for (int b = a + 1; b < n; b++) {
                int sb = idx[b];
                float gb0 = sites[k].gl[sb][0];
                float gb1 = sites[k].gl[sb][1];
                float gb2 = sites[k].gl[sb][2];
                double ibs = (double)(ga0*gb0 + ga1*gb1 + ga2*gb2);
                running_buf[p++] *= ibs;
            }
        }

        double ehh = 0;
        for (int q = 0; q < n_pairs; q++) ehh += running_buf[q];
        ehh /= n_pairs;

        int dpos = abs(sites[k].pos - prev_pos);
        ihh += 0.5 * (ehh_prev + ehh) * dpos;
        ehh_prev = ehh;
        prev_pos = sites[k].pos;

        if (ehh < ehh_decay) break;
        k += direction;
    }
    return ihh;
}

// ── Per-window aggregation containers ──────────────────────────────────────

typedef struct { int start, end, n_focal; } WinAgg;

typedef struct {
    int pos;
    double iHH_test, iHH_ref, xpehh, norm_xpehh;
    int    ok;
} FocalResult;

// ── Main ────────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "xpehh — per-window XP-EHH from BEAGLE GLs (diploid-IBS approximation).\n"
        "\n"
        "  --beagle <f.beagle.gz>     Required.\n"
        "  --sample_list <f>          Sample IDs in BEAGLE order.\n"
        "  --test_samples <f>         Test-cohort sample list.\n"
        "  --ref_samples <f>          Reference-cohort sample list.\n"
        "  --chr <name>               Filter to one chromosome.\n"
        "  --windows <bed>            BED of windows (chrom start end [id]).\n"
        "  --fixed_win W:S            Fixed-bp windows of W, step S.\n"
        "  --ehh_decay F              Stop walk when EHH < F (default 0.05).\n"
        "  --max_distance_bp N        Per-arm max walk distance (default 1e6).\n"
        "  --min_n_test N             Skip if n_test < N (default 20).\n"
        "  --min_n_ref N              Skip if n_ref < N (default 20).\n"
        "  --focal_step K             Compute every Kth SNP (default 1).\n"
        "  --out <f>                  Output TSV.\n"
        "  --ncores N                 OpenMP across focal SNPs.\n"
        "\n"
        "Output (schema_version=" SCHEMA_VERSION "):\n"
        "  chrom window_start window_end n_focal n_focal_with_xpehh\n"
        "  iHH_test_mean iHH_ref_mean xpehh_mean xpehh_max_abs\n"
        "  norm_xpehh_mean norm_xpehh_max_abs\n"
        "\n"
        "Caveat: this is a diploid-IBS-likelihood approximation of XP-EHH from\n"
        "unphased BEAGLE GLs. Sabeti-strict XP-EHH needs phasing. Sign convention\n"
        "matches Sabeti: positive = test cohort has longer EHH.\n");
}

int main(int argc, char** argv) {
    const char *beagle = NULL, *sample_path = NULL,
               *test_path = NULL, *ref_path = NULL,
               *filter_chr = NULL, *win_path = NULL, *out_path = NULL;
    int fixed_win = 0, fixed_step = 0;
    double ehh_decay = 0.05;
    int max_distance_bp = 1000000;
    int min_n_test = 20, min_n_ref = 20;
    int focal_step = 1;
    int ncores = 1;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--beagle")          && i+1<argc) beagle = argv[++i];
        else if (!strcmp(argv[i], "--sample_list")     && i+1<argc) sample_path = argv[++i];
        else if (!strcmp(argv[i], "--test_samples")    && i+1<argc) test_path = argv[++i];
        else if (!strcmp(argv[i], "--ref_samples")     && i+1<argc) ref_path = argv[++i];
        else if (!strcmp(argv[i], "--chr")             && i+1<argc) filter_chr = argv[++i];
        else if (!strcmp(argv[i], "--windows")         && i+1<argc) win_path = argv[++i];
        else if (!strcmp(argv[i], "--fixed_win")       && i+1<argc) sscanf(argv[++i], "%d:%d", &fixed_win, &fixed_step);
        else if (!strcmp(argv[i], "--ehh_decay")       && i+1<argc) ehh_decay = atof(argv[++i]);
        else if (!strcmp(argv[i], "--max_distance_bp") && i+1<argc) max_distance_bp = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--min_n_test")      && i+1<argc) min_n_test = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--min_n_ref")       && i+1<argc) min_n_ref = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--focal_step")      && i+1<argc) focal_step = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out")             && i+1<argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "--ncores")          && i+1<argc) ncores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[xpehh] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }
    if (!beagle || !sample_path || !test_path || !ref_path) { print_usage(); return 1; }

    #ifdef _OPENMP
    omp_set_num_threads(ncores > 0 ? ncores : 1);
    #else
    (void)ncores;
    #endif

    static char names[MAX_IND][64];
    int n_ind = load_samples(sample_path, names);
    fprintf(stderr, "[xpehh] %d samples\n", n_ind);
    if (n_ind == 0) return 2;

    int* test_idx = (int*)malloc(MAX_IND * sizeof(int));
    int* ref_idx  = (int*)malloc(MAX_IND * sizeof(int));
    int n_test = load_cohort(test_path, test_idx, names, n_ind);
    int n_ref  = load_cohort(ref_path,  ref_idx,  names, n_ind);
    fprintf(stderr, "[xpehh] test=%d ref=%d\n", n_test, n_ref);
    if (n_test < min_n_test || n_ref < min_n_ref) {
        fprintf(stderr, "[xpehh] Cohort too small (min_n_test=%d min_n_ref=%d)\n",
                min_n_test, min_n_ref);
        return 3;
    }

    Site* sites = NULL;
    int site_cap = 0;
    char chrom[32] = "";
    int n_sites = load_beagle_gl(beagle, filter_chr, &sites, &site_cap, n_ind, chrom, sizeof(chrom));
    if (n_sites == 0) { free(sites); return 0; }

    fprintf(stderr, "[xpehh] Computing XP-EHH at every %d SNP (decay %.4g, max_dist %d bp)\n",
            focal_step, ehh_decay, max_distance_bp);

    int n_pairs_test = n_test * (n_test - 1) / 2;
    int n_pairs_ref  = n_ref  * (n_ref  - 1) / 2;
    int max_pairs = n_pairs_test > n_pairs_ref ? n_pairs_test : n_pairs_ref;

    int n_focal = (n_sites + focal_step - 1) / focal_step;
    FocalResult* results = (FocalResult*)calloc((size_t)n_focal, sizeof(FocalResult));

    #pragma omp parallel
    {
        double* run_test = (double*)malloc((size_t)n_pairs_test * sizeof(double));
        double* run_ref  = (double*)malloc((size_t)n_pairs_ref  * sizeof(double));

        #pragma omp for schedule(dynamic, 64)
        for (int fi = 0; fi < n_focal; fi++) {
            int f = fi * focal_step;
            if (f >= n_sites) continue;
            results[fi].pos = sites[f].pos;

            double iHH_test_R = walk_one_arm(sites, n_sites, f, +1, test_idx, n_test,
                                              ehh_decay, max_distance_bp, max_pairs, run_test);
            double iHH_test_L = walk_one_arm(sites, n_sites, f, -1, test_idx, n_test,
                                              ehh_decay, max_distance_bp, max_pairs, run_test);
            double iHH_ref_R  = walk_one_arm(sites, n_sites, f, +1, ref_idx, n_ref,
                                              ehh_decay, max_distance_bp, max_pairs, run_ref);
            double iHH_ref_L  = walk_one_arm(sites, n_sites, f, -1, ref_idx, n_ref,
                                              ehh_decay, max_distance_bp, max_pairs, run_ref);
            results[fi].iHH_test = iHH_test_R + iHH_test_L;
            results[fi].iHH_ref  = iHH_ref_R  + iHH_ref_L;

            if (results[fi].iHH_test > 0 && results[fi].iHH_ref > 0) {
                results[fi].xpehh = log(results[fi].iHH_test) - log(results[fi].iHH_ref);
                results[fi].ok = 1;
            } else {
                results[fi].xpehh = NAN;
                results[fi].ok = 0;
            }
        }
        free(run_test); free(run_ref);
    }

    // Genome-wide z-normalization
    double sum_x = 0, sum_x2 = 0;
    int n_ok = 0;
    for (int i = 0; i < n_focal; i++) {
        if (!results[i].ok) continue;
        sum_x += results[i].xpehh;
        sum_x2 += results[i].xpehh * results[i].xpehh;
        n_ok++;
    }
    double global_mean = (n_ok > 0) ? sum_x / n_ok : 0;
    double var = (n_ok > 1) ? (sum_x2 - n_ok * global_mean * global_mean) / (n_ok - 1) : 0;
    double global_sd = (var > 0) ? sqrt(var) : 1.0;
    for (int i = 0; i < n_focal; i++) {
        if (!results[i].ok) { results[i].norm_xpehh = NAN; continue; }
        results[i].norm_xpehh = (results[i].xpehh - global_mean) / global_sd;
    }
    fprintf(stderr, "[xpehh] n_focal=%d n_with_xpehh=%d global μ=%.4g σ=%.4g\n",
            n_focal, n_ok, global_mean, global_sd);

    // Build window list
    typedef struct { int start, end; char id[64]; } Window;
    Window* wins = NULL;
    int n_wins = 0;
    if (win_path) {
        FILE* f = fopen(win_path, "r");
        if (!f) { fprintf(stderr, "[xpehh] Cannot open windows %s\n", win_path); return 5; }
        int cap = 100000;
        wins = (Window*)malloc((size_t)cap * sizeof(Window));
        char line[512];
        while (fgets(line, sizeof(line), f) && n_wins < cap) {
            char ch[64], id[64] = "";
            int s, e;
            int nf = sscanf(line, "%63s %d %d %63s", ch, &s, &e, id);
            if (nf < 3) continue;
            if (filter_chr && strcmp(ch, filter_chr) != 0) continue;
            wins[n_wins].start = s; wins[n_wins].end = e;
            if (id[0]) { strncpy(wins[n_wins].id, id, 63); wins[n_wins].id[63] = 0; }
            else snprintf(wins[n_wins].id, sizeof(wins[n_wins].id), "w%d", n_wins+1);
            n_wins++;
        }
        fclose(f);
    } else if (fixed_win > 0) {
        int max_pos = sites[n_sites - 1].pos;
        int cap = (max_pos / fixed_step) + 2;
        wins = (Window*)malloc((size_t)cap * sizeof(Window));
        for (int s = 0; s <= max_pos && n_wins < cap; s += fixed_step) {
            wins[n_wins].start = s; wins[n_wins].end = s + fixed_win - 1;
            snprintf(wins[n_wins].id, sizeof(wins[n_wins].id), "w%d", n_wins + 1);
            n_wins++;
        }
    } else {
        n_wins = 1;
        wins = (Window*)malloc(sizeof(Window));
        wins[0].start = sites[0].pos;
        wins[0].end = sites[n_sites-1].pos;
        snprintf(wins[0].id, sizeof(wins[0].id), "%s_all", chrom[0] ? chrom : "all");
    }
    fprintf(stderr, "[xpehh] %d windows\n", n_wins);

    // Aggregate per window
    FILE* fout = out_path ? fopen(out_path, "w") : stdout;
    if (!fout) { fprintf(stderr, "[xpehh] Cannot open --out %s\n", out_path); return 6; }

    fprintf(fout, "# schema_version=%s n_test=%d n_ref=%d ehh_decay=%.4g "
                   "max_distance_bp=%d focal_step=%d global_mean=%.6g global_sd=%.6g\n",
            SCHEMA_VERSION, n_test, n_ref, ehh_decay, max_distance_bp, focal_step,
            global_mean, global_sd);
    fprintf(fout, "chrom\twindow_start\twindow_end\tn_focal\tn_focal_with_xpehh\t"
                   "iHH_test_mean\tiHH_ref_mean\txpehh_mean\txpehh_max_abs\t"
                   "norm_xpehh_mean\tnorm_xpehh_max_abs\n");

    for (int w = 0; w < n_wins; w++) {
        int ws = wins[w].start, we = wins[w].end;
        int n_in = 0, n_with = 0;
        double sum_iHH_t = 0, sum_iHH_r = 0;
        double sum_xp = 0, sum_xp_norm = 0;
        double max_abs = 0, max_abs_norm = 0;
        for (int i = 0; i < n_focal; i++) {
            if (results[i].pos < ws || results[i].pos > we) continue;
            n_in++;
            if (!results[i].ok) continue;
            n_with++;
            sum_iHH_t += results[i].iHH_test;
            sum_iHH_r += results[i].iHH_ref;
            sum_xp += results[i].xpehh;
            sum_xp_norm += results[i].norm_xpehh;
            if (fabs(results[i].xpehh) > max_abs) max_abs = fabs(results[i].xpehh);
            if (fabs(results[i].norm_xpehh) > max_abs_norm) max_abs_norm = fabs(results[i].norm_xpehh);
        }

        fprintf(fout, "%s\t%d\t%d\t%d\t%d", chrom[0] ? chrom : ".", ws, we, n_in, n_with);
        #define EM(v) do { \
            if (isnan(v) || (n_with == 0 && !isfinite(v))) fprintf(fout, "\tNA"); \
            else fprintf(fout, "\t%.6g", v); \
        } while (0)
        if (n_with > 0) {
            EM(sum_iHH_t / n_with);
            EM(sum_iHH_r / n_with);
            EM(sum_xp / n_with);
            EM(max_abs);
            EM(sum_xp_norm / n_with);
            EM(max_abs_norm);
        } else {
            fputs("\tNA\tNA\tNA\tNA\tNA\tNA", fout);
        }
        #undef EM
        fputc('\n', fout);
    }

    if (out_path) fclose(fout);
    free(wins); free(results);
    free(test_idx); free(ref_idx);
    free(sites);
    return 0;
}
