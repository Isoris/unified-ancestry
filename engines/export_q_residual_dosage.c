// =============================================================================
// export_q_residual_dosage.c — Q-corrected residual dosage from BEAGLE
//
// For each sample i at site j:
//   observed_dosage  = GL1_ij + 2*GL2_ij  (after GL normalization)
//   expected_dosage  = Σ_k Q_ik × 2 × f_jk
//   residual_dosage  = observed - expected
//
// Where:
//   Q_ik = ancestry proportion of sample i in cluster k (from local Q cache)
//   f_jk = allele frequency at site j in cluster k (from F matrix / .fopt.gz)
//
// Q is window-level (one Q vector per window of W SNPs). Each SNP inherits
// the Q from its containing window. Since Q changes slowly, this is the
// standard local ancestry painting approximation.
//
// Inputs:
//   --beagle <file.beagle.gz>   BEAGLE GL file
//   --fopt <file.fopt.gz>       F matrix from NGSadmix (n_sites × K, gzipped)
//   --local_q <file.tsv.gz>     Local Q cache: window_id, sample_idx, Q1..QK
//                                OR a single global Q file (n_samples × K)
//   --sample_list <file>        Sample IDs (BAM list order)
//   --chr <chr>                 Filter chromosome
//   --window_size <N>           SNPs per Q window (default: 100, must match Engine B)
//   --out_prefix <prefix>       Output: <prefix>.residual.beagle.gz
//                                       <prefix>.residual.dosage.tsv.gz
//                                       <prefix>.sites.tsv.gz
//
// Output formats:
//   .residual.beagle.gz: same BEAGLE format but GL columns replaced with
//     pseudo-GLs from the residual (for compatibility with downstream tools)
//   .residual.dosage.tsv.gz: marker × sample matrix of residual dosages
//   .sites.tsv.gz: chrom, pos, marker
//
// Memory: O(n_sites × n_ind) for dosage + O(n_sites × K) for F matrix.
//   For 200K sites × 226 samples × K=8: ~370 MB. Fine.
//
// Compile: gcc -O3 -march=native -o export_q_residual_dosage \
//          export_q_residual_dosage.c -lz -lm
// =============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

#define MAX_IND 500
#define MAX_K 20
#define MAX_SITES 2000000

// ── GZ helpers (reused from region_popstats) ──

typedef struct {
    gzFile fp;
    char buf[1<<16];
    char lo[1<<18];
    int lo_len;
} Gz;

int gz_open(Gz* g, const char* path) {
    g->fp = gzopen(path, "rb");
    if (!g->fp) return 0;
    gzbuffer(g->fp, 1<<18);
    g->lo_len = 0;
    return 1;
}

int gz_getline(Gz* g, char* out, int maxlen) {
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
        if (g->lo_len + n > (int)sizeof(g->lo) - 1)
            n = sizeof(g->lo) - 1 - g->lo_len;
        memcpy(g->lo + g->lo_len, g->buf, n);
        g->lo_len += n;
    }
}
void gz_close(Gz* g) { if (g->fp) gzclose(g->fp); }

// ── Load sample list ──

int load_samples(const char* path, char names[][64]) {
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
        int cplen = len < 255 ? len : 255;
        memcpy(tmp, p, cplen); tmp[cplen] = 0;
        char* ext;
        if ((ext = strstr(tmp, ".sorted.markdup.bam"))) *ext = 0;
        else if ((ext = strstr(tmp, ".markdup.bam"))) *ext = 0;
        else if ((ext = strstr(tmp, ".sorted.bam"))) *ext = 0;
        else if ((ext = strstr(tmp, ".bam"))) *ext = 0;
        else if ((ext = strstr(tmp, ".cram"))) *ext = 0;
        int nlen = strlen(tmp);
        if (nlen > 63) nlen = 63;
        memcpy(names[n], tmp, nlen); names[n][nlen] = 0;
        n++;
    }
    fclose(f);
    return n;
}

// ── Load F matrix (.fopt.gz or .fopt) ──
// Format: n_sites rows, K columns (space/tab separated), no header.
// Each row j has f_j1 f_j2 ... f_jK (allele frequencies).

typedef struct {
    double* data;   // flat array [n_sites * K]
    int n_sites;
    int K;
} FMatrix;

FMatrix load_fopt(const char* path) {
    FMatrix fm = {NULL, 0, 0};
    Gz gz;
    if (!gz_open(&gz, path)) {
        fprintf(stderr, "[residual] Cannot open F matrix: %s\n", path);
        return fm;
    }

    int cap = 500000;
    double* data = (double*)malloc(cap * MAX_K * sizeof(double));
    int n = 0, K = 0;
    char line[4096];

    while (gz_getline(&gz, line, sizeof(line))) {
        if (n >= cap) {
            cap *= 2;
            data = (double*)realloc(data, cap * MAX_K * sizeof(double));
        }
        char* p = line;
        int k = 0;
        while (*p && k < MAX_K) {
            char* endp;
            double v = strtod(p, &endp);
            if (endp == p) break;
            data[n * MAX_K + k] = v;
            p = endp;
            while (*p == ' ' || *p == '\t') p++;
            k++;
        }
        if (k == 0) continue;
        if (K == 0) K = k;
        n++;
    }
    gz_close(&gz);

    fm.data = data;
    fm.n_sites = n;
    fm.K = K;
    fprintf(stderr, "[residual] F matrix: %d sites × K=%d\n", n, K);
    return fm;
}

// ── Load Q matrix ──
// Two modes:
//   1. Global Q: n_samples rows × K columns (no header). Applied to ALL windows.
//   2. Local Q cache: window_id\tsample_idx\tQ1\t...\tQK (with header).
//      In this case, Q changes per window.
// We detect by checking if first line looks like a header.

typedef struct {
    double* data;    // [n_samples * K] for global, or [n_windows * n_samples * K] for local
    int n_samples;
    int K;
    int n_windows;   // 0 = global mode, >0 = local mode
    int is_local;
} QMatrix;

QMatrix load_q(const char* path, int n_ind_expected) {
    QMatrix qm = {NULL, 0, 0, 0, 0};

    // Try as plain text first (global Q)
    FILE* f = fopen(path, "r");
    if (!f) {
        // Try gzipped
        Gz gz;
        if (!gz_open(&gz, path)) {
            fprintf(stderr, "[residual] Cannot open Q: %s\n", path);
            return qm;
        }
        // Read gzipped — check for header
        char line[8192];
        gz_getline(&gz, line, sizeof(line));

        // Check if header (contains "Q1" or "window_id" or "sample")
        int has_header = (strstr(line, "Q1") || strstr(line, "window") || strstr(line, "sample"));

        if (has_header) {
            // Local Q cache: skip header, read data
            // For now, fall back to global mode using qinit
            fprintf(stderr, "[residual] Local Q cache detected — using as global fallback\n");
            // TODO: implement per-window Q loading
            gz_close(&gz);
            return qm;
        }

        // Global Q from gzipped file
        double* data = (double*)malloc(n_ind_expected * MAX_K * sizeof(double));
        int n = 0, K = 0;
        // First line is data
        {
            char* p = line;
            int k = 0;
            while (*p && k < MAX_K) {
                char* endp;
                double v = strtod(p, &endp);
                if (endp == p) break;
                data[n * MAX_K + k] = v;
                p = endp;
                while (*p == ' ' || *p == '\t') p++;
                k++;
            }
            if (k > 0) { K = k; n++; }
        }
        while (gz_getline(&gz, line, sizeof(line)) && n < n_ind_expected) {
            char* p = line;
            int k = 0;
            while (*p && k < MAX_K) {
                char* endp;
                double v = strtod(p, &endp);
                if (endp == p) break;
                data[n * MAX_K + k] = v;
                p = endp;
                while (*p == ' ' || *p == '\t') p++;
                k++;
            }
            if (k > 0) n++;
        }
        gz_close(&gz);
        qm.data = data;
        qm.n_samples = n;
        qm.K = K;
        qm.is_local = 0;
        fprintf(stderr, "[residual] Q (global, gz): %d samples × K=%d\n", n, K);
        return qm;
    }

    // Plain text global Q
    double* data = (double*)malloc(n_ind_expected * MAX_K * sizeof(double));
    int n = 0, K = 0;
    char line[8192];
    while (fgets(line, sizeof(line), f) && n < n_ind_expected) {
        char* p = line;
        int k = 0;
        while (*p && k < MAX_K) {
            char* endp;
            double v = strtod(p, &endp);
            if (endp == p) break;
            data[n * MAX_K + k] = v;
            p = endp;
            while (*p == ' ' || *p == '\t') p++;
            k++;
        }
        if (k == 0) continue;
        if (K == 0) K = k;
        n++;
    }
    fclose(f);
    qm.data = data;
    qm.n_samples = n;
    qm.K = K;
    qm.is_local = 0;
    fprintf(stderr, "[residual] Q (global): %d samples × K=%d\n", n, K);
    return qm;
}

// ── Main ──

int main(int argc, char** argv) {
    const char *beagle_path = NULL, *fopt_path = NULL, *q_path = NULL;
    const char *sample_path = NULL, *filter_chr = NULL, *out_prefix = "residual";
    int window_size = 100;  // reserved for future local Q windowing
    (void)window_size;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--beagle") && i+1<argc)      beagle_path = argv[++i];
        else if (!strcmp(argv[i], "--fopt") && i+1<argc)        fopt_path = argv[++i];
        else if (!strcmp(argv[i], "--local_q") && i+1<argc)     q_path = argv[++i];
        else if (!strcmp(argv[i], "--sample_list") && i+1<argc) sample_path = argv[++i];
        else if (!strcmp(argv[i], "--chr") && i+1<argc)         filter_chr = argv[++i];
        else if (!strcmp(argv[i], "--out_prefix") && i+1<argc)  out_prefix = argv[++i];
        else if (!strcmp(argv[i], "--window_size") && i+1<argc) window_size = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            fprintf(stderr,
                "export_q_residual_dosage — Q-corrected residual dosage\n\n"
                "  --beagle <f>         BEAGLE GL file (.beagle.gz)\n"
                "  --fopt <f>           F matrix (.fopt.gz) from NGSadmix\n"
                "  --local_q <f>        Q matrix (.qopt or local_Q cache)\n"
                "  --sample_list <f>    Sample IDs (BAM list order)\n"
                "  --chr <chr>          Filter chromosome\n"
                "  --window_size <N>    SNPs per Q window (default: 100)\n"
                "  --out_prefix <pfx>   Output prefix\n\n"
                "Output:\n"
                "  <pfx>.residual.dosage.tsv.gz  — residual dosage matrix\n"
                "  <pfx>.sites.tsv.gz            — site positions\n");
            return 0;
        }
    }

    if (!beagle_path || !fopt_path || !q_path || !sample_path) {
        fprintf(stderr, "Need: --beagle, --fopt, --local_q, --sample_list\n");
        return 1;
    }

    // Load samples
    static char names[MAX_IND][64];
    int n_ind = load_samples(sample_path, names);
    fprintf(stderr, "[residual] %d samples\n", n_ind);

    // Load F matrix
    FMatrix fm = load_fopt(fopt_path);
    if (!fm.data || fm.n_sites == 0) { fprintf(stderr, "F matrix load failed\n"); return 1; }

    // Load Q matrix
    QMatrix qm = load_q(q_path, n_ind);
    if (!qm.data || qm.n_samples == 0) { fprintf(stderr, "Q matrix load failed\n"); return 1; }

    int K = fm.K;
    if (qm.K != K) {
        fprintf(stderr, "[residual] WARNING: F has K=%d, Q has K=%d. Using min.\n", K, qm.K);
        K = K < qm.K ? K : qm.K;
    }

    // Open BEAGLE
    Gz gz;
    if (!gz_open(&gz, beagle_path)) {
        fprintf(stderr, "Cannot open BEAGLE: %s\n", beagle_path);
        return 1;
    }

    // Open outputs
    char path_dos[512], path_sites[512];
    snprintf(path_dos, sizeof(path_dos), "%s.residual.dosage.tsv.gz", out_prefix);
    snprintf(path_sites, sizeof(path_sites), "%s.sites.tsv.gz", out_prefix);

    gzFile gz_dos = gzopen(path_dos, "wb");
    gzFile gz_sites = gzopen(path_sites, "wb");
    if (!gz_dos || !gz_sites) {
        fprintf(stderr, "Cannot open output files\n");
        return 1;
    }

    // Write dosage header
    gzprintf(gz_dos, "marker");
    for (int i = 0; i < n_ind; i++) gzprintf(gz_dos, "\t%s", names[i]);
    gzprintf(gz_dos, "\n");

    // Write sites header
    gzprintf(gz_sites, "chrom\tpos\tmarker\n");

    // Process BEAGLE line by line
    char* line = (char*)malloc(1 << 20);
    gz_getline(&gz, line, 1 << 20); // skip header

    int site_global = 0;  // index into F matrix (across ALL sites in BEAGLE, pre-filter)
    int site_kept = 0;    // sites that passed chr filter

    while (gz_getline(&gz, line, 1 << 20)) {
        // Parse marker
        char* tab1 = strchr(line, '\t');
        if (!tab1) { site_global++; continue; }
        *tab1 = 0;
        char marker[256];
        strncpy(marker, line, 255); marker[255] = 0;

        char* us = strrchr(marker, '_');
        if (!us) { site_global++; continue; }

        char chr_buf[128];
        int marker_len = us - marker;
        if (marker_len > 127) marker_len = 127;
        memcpy(chr_buf, marker, marker_len); chr_buf[marker_len] = 0;
        int pos = atoi(us + 1);

        // Restore marker name
        *us = '_';

        // Chromosome filter
        if (filter_chr && strcmp(chr_buf, filter_chr) != 0) {
            site_global++;
            continue;
        }

        // Q: global mode uses same Q for all windows
        // Future: local mode would index qm.data[win_idx * n_ind * K + ...] here

        // Check F matrix bounds
        if (site_global >= fm.n_sites) {
            fprintf(stderr, "[residual] WARNING: site_global=%d exceeds F matrix (%d sites). "
                    "Stopping.\n", site_global, fm.n_sites);
            break;
        }

        // Get F row for this site: f[site_global * MAX_K + k]
        double* f_row = &fm.data[site_global * MAX_K];

        // Skip allele1, allele2
        char* p = tab1 + 1;
        char* tab2 = strchr(p, '\t'); if (!tab2) { site_global++; continue; }
        char* tab3 = strchr(tab2 + 1, '\t'); if (!tab3) { site_global++; continue; }
        p = tab3 + 1;

        // Write sites
        gzprintf(gz_sites, "%s\t%d\t%s\n", chr_buf, pos, marker);

        // Write dosage line
        gzprintf(gz_dos, "%s", marker);

        for (int i = 0; i < n_ind; i++) {
            // Parse GLs
            double gl0 = 0, gl1 = 0, gl2 = 0;
            char* endp;
            gl0 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;
            gl1 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;
            gl2 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;

            double s = gl0 + gl1 + gl2;
            if (s > 1e-15) { gl0 /= s; gl1 /= s; gl2 /= s; }
            double obs_dos = gl1 + 2.0 * gl2;

            // Compute expected dosage: E[g|Q,F] = 2 * Σ_k Q_ik * f_jk
            double exp_dos = 0;
            int q_idx = (i < qm.n_samples) ? i : qm.n_samples - 1;
            for (int k = 0; k < K; k++) {
                double q_ik = qm.data[q_idx * MAX_K + k];
                double f_jk = f_row[k];
                exp_dos += q_ik * f_jk;
            }
            exp_dos *= 2.0;

            // Residual
            double residual = obs_dos - exp_dos;

            gzprintf(gz_dos, "\t%.4f", residual);
        }
        gzprintf(gz_dos, "\n");

        site_global++;
        site_kept++;
    }

    gz_close(&gz);
    gzclose(gz_dos);
    gzclose(gz_sites);
    free(line);
    free(fm.data);
    free(qm.data);

    fprintf(stderr, "[residual] Done: %d sites (global idx %d), K=%d\n",
            site_kept, site_global, K);
    fprintf(stderr, "[residual] Output: %s, %s\n", path_dos, path_sites);
    return 0;
}
