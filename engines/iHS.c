// =============================================================================
// iHS.c — per-window iHS (intra-cohort EHH-based selection scan) from BEAGLE GLs.
//
// ── What this binary computes ───────────────────────────────────────────────
//
// iHS (Voight, Kudaravalli, Wen, Pritchard 2006) is a single-cohort selection
// scan: at each focal SNP it asks whether the DERIVED-allele homozygotes share
// longer extended haplotype homozygosity than the ANCESTRAL-allele homozygotes.
// Negative iHS = derived allele has longer EHH → consistent with recent
// positive selection on the derived allele.
//
// Voight-strict iHS needs PHASED haplotypes. This binary computes a
// defensible approximation directly from BEAGLE GLs without phasing, using
// the same diploid-IBS-likelihood EHH machinery as xpehh.c. Within-cohort,
// samples are partitioned at each focal SNP into:
//
//   ANC sub-cohort: argmax-genotype == HOM_REF at focal AND max GL ≥ --gl_threshold
//   DER sub-cohort: argmax-genotype == HOM_ALT at focal AND max GL ≥ --gl_threshold
//   (heterozygotes and low-confidence sites are dropped at the focal)
//
//   IBS_prob_ij(k) = gl0_i·gl0_j + gl1_i·gl1_j + gl2_i·gl2_j   ∈ [0, 1]
//   EHH_sub(focal, d) = mean over within-subcohort pairs of running ∏ IBS_prob.
//   iHH_sub(focal)    = ∫ EHH_sub(d) dd   (trapezoidal in bp, walked both arms,
//                        stopping at EHH < --ehh_decay or > --max_distance_bp).
//
//   iHS_raw(focal)  = ln( iHH_anc(focal) / iHH_der(focal) )
//
// Frequency-binned standardization (Voight 2006 §Methods):
//   Bin focals by derived-allele frequency into K bins.
//   Within each bin: compute mean(iHS_raw) and sd(iHS_raw).
//   iHS_norm(focal) = (iHS_raw − bin_mean) / bin_sd.
//
// Per-window aggregation: each focal SNP belongs to the window containing its
// position. Aggregates: mean and max-|·| of iHS_raw and iHS_norm.
//
// Ancestral polarization:
//   By default, BEAGLE allele 0 is treated as ancestral. To polarize against
//   an outgroup-derived ancestral state, pre-process the BEAGLE so allele 0 =
//   ancestral at every site (out of scope for this binary).
//
// Caveats (documented):
//   - Not Voight-strict (no phasing). Heterozygous-heterozygous pairs hit
//     IBS_prob ≈ 0.5 even when haplotypes match → conservative EHH.
//   - Sign convention matches Voight: negative iHS = derived has longer EHH
//     → potential recent positive selection on derived.
//   - Sub-cohort sizes vary per focal; at low-derived-frequency focals the DER
//     sub-cohort may be empty (skipped). --min_n_anc / --min_n_der gate this.
//
// ── Inputs ──────────────────────────────────────────────────────────────────
//   --beagle <f.beagle.gz>     BEAGLE GL file (allele 0 treated as ancestral)
//   --sample_list <f>          Sample IDs in BEAGLE order
//   --cohort <f>               One sample ID per line — the cohort to scan
//                              (default: all samples in --sample_list)
//   --chr <name>               Filter to a single chromosome
//   --windows <bed>            BED of windows (chrom start end [id])
//   --fixed_win W:S            Fixed-bp windows
//   --ehh_decay F              Stop walk when EHH < F (default 0.05)
//   --max_distance_bp N        Per-arm max walk distance (default 1e6)
//   --gl_threshold F           Min max-GL for confident genotype (default 0.6)
//   --min_n_anc N              Skip focal if anc sub-cohort < N (default 8)
//   --min_n_der N              Skip focal if der sub-cohort < N (default 8)
//   --min_freq F               Skip focal if MAF < F (default 0.05)
//   --n_freq_bins K            Number of derived-allele-freq bins for norm (default 20)
//   --focal_step K             Compute every Kth SNP (default 1)
//   --out <f>                  Output TSV
//   --ncores N                 OpenMP across focal SNPs
//
// ── Outputs ─────────────────────────────────────────────────────────────────
//   Per-window TSV (schema_version=iHS_v1):
//     chrom window_start window_end n_focal n_focal_with_iHS
//     iHH_anc_mean iHH_der_mean iHS_mean iHS_max_abs
//     norm_iHS_mean norm_iHS_max_abs n_extreme_iHS_abs_gt_2
//
// Compile: gcc -O3 -march=native -fopenmp -o iHS iHS.c -lz -lm
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
#define SCHEMA_VERSION "iHS_v1"

// ── Gz reader (same as xpehh.c / region_popstats.c) ─────────────────────────

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

// ── Sample loaders ──────────────────────────────────────────────────────────

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

// ── Site storage (same model as xpehh.c, grow-on-demand) ───────────────────

typedef struct {
    int pos;
    float gl[MAX_IND][3];
} Site;

static int load_beagle_gl(const char* path, const char* filter_chr,
                          Site** sites_out, int* cap_out,
                          int n_ind, char* chrom_out, int chrom_out_sz) {
    Gz gz;
    if (!gz_open(&gz, path)) { fprintf(stderr, "[iHS] Cannot open %s\n", path); return 0; }
    char line[1<<20];
    gz_getline(&gz, line, sizeof(line));

    int cap = 65536;
    Site* sites = (Site*)malloc((size_t)cap * sizeof(Site));
    if (!sites) { fprintf(stderr, "[iHS] initial malloc failed\n"); gz_close(&gz); return 0; }

    int n = 0;
    while (gz_getline(&gz, line, sizeof(line))) {
        if (n == cap) {
            int new_cap = cap * 2;
            Site* tmp = (Site*)realloc(sites, (size_t)new_cap * sizeof(Site));
            if (!tmp) { fprintf(stderr, "[iHS] realloc to %d sites failed\n", new_cap); break; }
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
        if (filter_chr && strcmp(chr, filter_chr) != 0) continue;
        if (chrom_out && chrom_out[0] == 0) {
            size_t cn = strlen(chr);
            if (cn >= (size_t)chrom_out_sz) cn = chrom_out_sz - 1;
            memcpy(chrom_out, chr, cn);
            chrom_out[cn] = 0;
        }
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
    fprintf(stderr, "[iHS] Loaded %d sites for %s (cap=%d)\n", n,
            filter_chr ? filter_chr : "all", cap);
    return n;
}

// ── EHH walking arm (identical to xpehh's; computes one iHH per cohort) ─────

static double walk_one_arm(const Site* sites, int n_sites, int focal, int direction,
                            const int* idx, int n,
                            double ehh_decay, int max_distance_bp,
                            double* running_buf) {
    if (n < 2) return 0;
    int n_pairs = n * (n - 1) / 2;

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

    int focal_pos = sites[focal].pos;
    int prev_pos = focal_pos;
    double ihh = 0;
    int k = focal + direction;
    while (k >= 0 && k < n_sites) {
        int dist = abs(sites[k].pos - focal_pos);
        if (dist > max_distance_bp) break;
        if (sites[k].pos == prev_pos) { k += direction; continue; }

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

// ── Focal results ───────────────────────────────────────────────────────────

typedef struct {
    int pos;
    int    n_anc, n_der;
    double der_freq;
    double iHH_anc, iHH_der;
    double iHS, iHS_norm;
    int    freq_bin;
    int    ok;
} FocalResult;

// ── Main ────────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "iHS — per-window iHS from BEAGLE GLs (Voight 2006 framing, diploid-IBS approx).\n"
        "\n"
        "  --beagle <f.beagle.gz>     Required.\n"
        "  --sample_list <f>          Sample IDs in BEAGLE order.\n"
        "  --cohort <f>               One sample ID per line (default: all samples).\n"
        "  --chr <name>               Filter to one chromosome.\n"
        "  --windows <bed>            BED of windows.\n"
        "  --fixed_win W:S            Fixed-bp windows.\n"
        "  --ehh_decay F              Stop walk when EHH < F (default 0.05).\n"
        "  --max_distance_bp N        Per-arm walk cap (default 1e6).\n"
        "  --gl_threshold F           Min max-GL for confident hom call (default 0.6).\n"
        "  --min_n_anc N              Skip focal if anc subcohort < N (default 8).\n"
        "  --min_n_der N              Skip focal if der subcohort < N (default 8).\n"
        "  --min_freq F               Skip focal if MAF < F (default 0.05).\n"
        "  --n_freq_bins K            Derived-allele-freq bins for norm (default 20).\n"
        "  --focal_step K             Compute every Kth SNP (default 1).\n"
        "  --out <f>                  Output TSV.\n"
        "  --ncores N                 OpenMP across focal SNPs.\n"
        "\n"
        "Output (schema_version=" SCHEMA_VERSION "):\n"
        "  chrom window_start window_end n_focal n_focal_with_iHS\n"
        "  iHH_anc_mean iHH_der_mean iHS_mean iHS_max_abs\n"
        "  norm_iHS_mean norm_iHS_max_abs n_extreme_iHS_abs_gt_2\n"
        "\n"
        "Caveat: diploid-IBS approximation. Voight-strict iHS needs phasing.\n"
        "        Allele 0 in BEAGLE is treated as ancestral.\n"
        "        Sign convention: negative = derived has longer EHH → potential\n"
        "        recent positive selection on derived.\n");
}

int main(int argc, char** argv) {
    const char *beagle = NULL, *sample_path = NULL, *cohort_path = NULL,
               *filter_chr = NULL, *win_path = NULL, *out_path = NULL;
    int fixed_win = 0, fixed_step = 0;
    double ehh_decay = 0.05;
    int max_distance_bp = 1000000;
    double gl_threshold = 0.6;
    int min_n_anc = 8, min_n_der = 8;
    double min_freq = 0.05;
    int n_freq_bins = 20;
    int focal_step = 1;
    int ncores = 1;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--beagle")          && i+1<argc) beagle = argv[++i];
        else if (!strcmp(argv[i], "--sample_list")     && i+1<argc) sample_path = argv[++i];
        else if (!strcmp(argv[i], "--cohort")          && i+1<argc) cohort_path = argv[++i];
        else if (!strcmp(argv[i], "--chr")             && i+1<argc) filter_chr = argv[++i];
        else if (!strcmp(argv[i], "--windows")         && i+1<argc) win_path = argv[++i];
        else if (!strcmp(argv[i], "--fixed_win")       && i+1<argc) sscanf(argv[++i], "%d:%d", &fixed_win, &fixed_step);
        else if (!strcmp(argv[i], "--ehh_decay")       && i+1<argc) ehh_decay = atof(argv[++i]);
        else if (!strcmp(argv[i], "--max_distance_bp") && i+1<argc) max_distance_bp = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--gl_threshold")    && i+1<argc) gl_threshold = atof(argv[++i]);
        else if (!strcmp(argv[i], "--min_n_anc")       && i+1<argc) min_n_anc = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--min_n_der")       && i+1<argc) min_n_der = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--min_freq")        && i+1<argc) min_freq = atof(argv[++i]);
        else if (!strcmp(argv[i], "--n_freq_bins")     && i+1<argc) n_freq_bins = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--focal_step")      && i+1<argc) focal_step = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out")             && i+1<argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "--ncores")          && i+1<argc) ncores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[iHS] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }
    if (!beagle || !sample_path) { print_usage(); return 1; }

    #ifdef _OPENMP
    omp_set_num_threads(ncores > 0 ? ncores : 1);
    #else
    (void)ncores;
    #endif

    static char names[MAX_IND][64];
    int n_ind = load_samples(sample_path, names);
    fprintf(stderr, "[iHS] %d samples\n", n_ind);
    if (n_ind == 0) return 2;

    int* cohort_idx = (int*)malloc(MAX_IND * sizeof(int));
    int n_cohort;
    if (cohort_path) {
        n_cohort = load_cohort(cohort_path, cohort_idx, names, n_ind);
    } else {
        n_cohort = n_ind;
        for (int i = 0; i < n_ind; i++) cohort_idx[i] = i;
    }
    fprintf(stderr, "[iHS] cohort=%d\n", n_cohort);

    Site* sites = NULL;
    int site_cap = 0;
    char chrom[32] = "";
    int n_sites = load_beagle_gl(beagle, filter_chr, &sites, &site_cap, n_ind, chrom, sizeof(chrom));
    if (n_sites == 0) { free(sites); free(cohort_idx); return 0; }

    fprintf(stderr, "[iHS] decay=%.4g max_dist=%d gl_thresh=%.4g min_freq=%.4g\n",
            ehh_decay, max_distance_bp, gl_threshold, min_freq);

    int n_focal = (n_sites + focal_step - 1) / focal_step;
    FocalResult* results = (FocalResult*)calloc((size_t)n_focal, sizeof(FocalResult));

    #pragma omp parallel
    {
        int* anc_idx = (int*)malloc((size_t)n_cohort * sizeof(int));
        int* der_idx = (int*)malloc((size_t)n_cohort * sizeof(int));
        int max_pairs = n_cohort * (n_cohort - 1) / 2;
        double* run_buf = (double*)malloc((size_t)max_pairs * sizeof(double));

        #pragma omp for schedule(dynamic, 64)
        for (int fi = 0; fi < n_focal; fi++) {
            int f = fi * focal_step;
            if (f >= n_sites) continue;
            results[fi].pos = sites[f].pos;

            // Partition the cohort at this focal.
            int n_anc = 0, n_der = 0;
            double sum_dos = 0;
            for (int k = 0; k < n_cohort; k++) {
                int s = cohort_idx[k];
                float g0 = sites[f].gl[s][0];
                float g1 = sites[f].gl[s][1];
                float g2 = sites[f].gl[s][2];
                sum_dos += g1 + 2.0 * g2;
                float maxg = g0; int arg = 0;
                if (g1 > maxg) { maxg = g1; arg = 1; }
                if (g2 > maxg) { maxg = g2; arg = 2; }
                if (maxg < gl_threshold) continue;
                if (arg == 0)      anc_idx[n_anc++] = s;
                else if (arg == 2) der_idx[n_der++] = s;
                // arg == 1 (heterozygote): dropped at focal
            }
            double der_freq = sum_dos / (2.0 * n_cohort);
            double maf = der_freq < 0.5 ? der_freq : 1.0 - der_freq;
            results[fi].n_anc = n_anc;
            results[fi].n_der = n_der;
            results[fi].der_freq = der_freq;
            results[fi].ok = 0;

            if (maf < min_freq || n_anc < min_n_anc || n_der < min_n_der) {
                results[fi].iHS = results[fi].iHS_norm = NAN;
                continue;
            }

            double iHH_anc_R = walk_one_arm(sites, n_sites, f, +1, anc_idx, n_anc,
                                             ehh_decay, max_distance_bp, run_buf);
            double iHH_anc_L = walk_one_arm(sites, n_sites, f, -1, anc_idx, n_anc,
                                             ehh_decay, max_distance_bp, run_buf);
            double iHH_der_R = walk_one_arm(sites, n_sites, f, +1, der_idx, n_der,
                                             ehh_decay, max_distance_bp, run_buf);
            double iHH_der_L = walk_one_arm(sites, n_sites, f, -1, der_idx, n_der,
                                             ehh_decay, max_distance_bp, run_buf);

            results[fi].iHH_anc = iHH_anc_R + iHH_anc_L;
            results[fi].iHH_der = iHH_der_R + iHH_der_L;
            if (results[fi].iHH_anc > 0 && results[fi].iHH_der > 0) {
                results[fi].iHS = log(results[fi].iHH_anc) - log(results[fi].iHH_der);
                results[fi].ok = 1;
            } else {
                results[fi].iHS = NAN;
            }
        }
        free(anc_idx); free(der_idx); free(run_buf);
    }

    // ── Frequency-binned normalization ──
    // Bin derived allele freq into n_freq_bins equal-width bins on [0,1].
    double* bin_sum  = (double*)calloc((size_t)n_freq_bins, sizeof(double));
    double* bin_sum2 = (double*)calloc((size_t)n_freq_bins, sizeof(double));
    int*    bin_n    = (int*)calloc((size_t)n_freq_bins, sizeof(int));
    for (int i = 0; i < n_focal; i++) {
        if (!results[i].ok) { results[i].freq_bin = -1; continue; }
        int b = (int)(results[i].der_freq * n_freq_bins);
        if (b < 0) b = 0;
        if (b >= n_freq_bins) b = n_freq_bins - 1;
        results[i].freq_bin = b;
        bin_sum[b]  += results[i].iHS;
        bin_sum2[b] += results[i].iHS * results[i].iHS;
        bin_n[b]++;
    }
    double* bin_mean = (double*)malloc((size_t)n_freq_bins * sizeof(double));
    double* bin_sd   = (double*)malloc((size_t)n_freq_bins * sizeof(double));
    for (int b = 0; b < n_freq_bins; b++) {
        bin_mean[b] = (bin_n[b] > 0) ? bin_sum[b] / bin_n[b] : 0;
        double var = (bin_n[b] > 1)
                     ? (bin_sum2[b] - bin_n[b] * bin_mean[b] * bin_mean[b]) / (bin_n[b] - 1) : 0;
        bin_sd[b] = (var > 0) ? sqrt(var) : 1.0;
    }
    int n_ok = 0;
    for (int i = 0; i < n_focal; i++) {
        if (!results[i].ok) { results[i].iHS_norm = NAN; continue; }
        int b = results[i].freq_bin;
        results[i].iHS_norm = (results[i].iHS - bin_mean[b]) / bin_sd[b];
        n_ok++;
    }
    fprintf(stderr, "[iHS] n_focal=%d n_with_iHS=%d (%d freq bins)\n",
            n_focal, n_ok, n_freq_bins);

    // Windows
    typedef struct { int start, end; char id[64]; } Window;
    Window* wins = NULL;
    int n_wins = 0;
    if (win_path) {
        FILE* f = fopen(win_path, "r");
        if (!f) { fprintf(stderr, "[iHS] Cannot open windows %s\n", win_path); return 5; }
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
    fprintf(stderr, "[iHS] %d windows\n", n_wins);

    FILE* fout = out_path ? fopen(out_path, "w") : stdout;
    if (!fout) { fprintf(stderr, "[iHS] Cannot open --out %s\n", out_path); return 6; }
    fprintf(fout, "# schema_version=%s n_cohort=%d ehh_decay=%.4g max_distance_bp=%d "
                   "gl_threshold=%.4g min_freq=%.4g n_freq_bins=%d focal_step=%d\n",
            SCHEMA_VERSION, n_cohort, ehh_decay, max_distance_bp,
            gl_threshold, min_freq, n_freq_bins, focal_step);
    fprintf(fout, "chrom\twindow_start\twindow_end\tn_focal\tn_focal_with_iHS\t"
                   "iHH_anc_mean\tiHH_der_mean\tiHS_mean\tiHS_max_abs\t"
                   "norm_iHS_mean\tnorm_iHS_max_abs\tn_extreme_iHS_abs_gt_2\n");

    for (int w = 0; w < n_wins; w++) {
        int ws = wins[w].start, we = wins[w].end;
        int n_in = 0, n_with = 0, n_extreme = 0;
        double s_anc = 0, s_der = 0, s_ihs = 0, s_norm = 0;
        double ma_ihs = 0, ma_norm = 0;
        for (int i = 0; i < n_focal; i++) {
            if (results[i].pos < ws || results[i].pos > we) continue;
            n_in++;
            if (!results[i].ok) continue;
            n_with++;
            s_anc  += results[i].iHH_anc;
            s_der  += results[i].iHH_der;
            s_ihs  += results[i].iHS;
            s_norm += results[i].iHS_norm;
            if (fabs(results[i].iHS) > ma_ihs)            ma_ihs  = fabs(results[i].iHS);
            if (fabs(results[i].iHS_norm) > ma_norm)      ma_norm = fabs(results[i].iHS_norm);
            if (fabs(results[i].iHS_norm) > 2.0)          n_extreme++;
        }
        fprintf(fout, "%s\t%d\t%d\t%d\t%d", chrom[0] ? chrom : ".", ws, we, n_in, n_with);
        #define EM(v) do { if (isnan(v)) fprintf(fout, "\tNA"); else fprintf(fout, "\t%.6g", (double)(v)); } while (0)
        if (n_with > 0) {
            EM(s_anc / n_with);
            EM(s_der / n_with);
            EM(s_ihs / n_with);
            EM(ma_ihs);
            EM(s_norm / n_with);
            EM(ma_norm);
        } else {
            fputs("\tNA\tNA\tNA\tNA\tNA\tNA", fout);
        }
        fprintf(fout, "\t%d\n", n_extreme);
        #undef EM
    }

    if (out_path) fclose(fout);
    free(wins); free(results);
    free(bin_sum); free(bin_sum2); free(bin_n); free(bin_mean); free(bin_sd);
    free(cohort_idx);
    free(sites);
    return 0;
}
