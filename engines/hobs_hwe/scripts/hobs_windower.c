// =============================================================================
// hobs_windower.c — Fast multi-scale Hobs/F windowed summaries
//
// Reads ANGSD .hwe.gz output, computes site-level Hobs from F + hweFreq
// (Mérot logic: Hexp = 2*p*(1-p), Hobs = Hexp*(1-F)), then aggregates
// into sliding windows at 7 scales with mean, median, SD, outlier burden.
//
// This replaces Mérot's R windowscanr approach for speed when running
// 7 scales × 28 chromosomes × ~10 subsets.
//
// Compile: gcc -O3 -o hobs_windower hobs_windower.c -lz -lm
//
// Usage:
//   hobs_windower <input.hwe.gz> <output_prefix> <chrom_size> \
//     [--scales 5000:1000,10000:2000,...] [--mad_n 3.0]
//
// Input: ANGSD .hwe.gz with columns:
//   Chromo Position Major Minor hweFreq Freq F LRT p-value [hetFreq]
//
// Output per scale:
//   <prefix>.win<label>.tsv  — windowed summaries
//   <prefix>.sites.tsv       — site-level Hobs/F (once, shared)
//
// Citation:
//   Hobs sliding-window concept: Claire Mérot (angsd_pipeline)
//   ANGSD HWE: Korneliussen et al. (2014)
//   Patched ANGSD: github.com/Isoris/angsd_fixed_HWE
// =============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

#define MAX_SITES 5000000
#define MAX_SCALES 10

typedef struct {
    int pos;
    double hweFreq;    // allele freq under HWE
    double freq;       // allele freq (EM with F)
    double F;          // inbreeding coefficient
    double Hexp;       // 2*p*(1-p) where p = hweFreq
    double Hobs;       // Hexp * (1 - F)
    double pval;       // HWE p-value
} Site;

typedef struct {
    char label[16];
    int win_bp;
    int step_bp;
} Scale;

// ── Comparison for qsort (median) ──
int cmp_double(const void* a, const void* b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

double median_sorted(double* arr, int n) {
    if (n == 0) return 0.0/0.0;
    if (n % 2 == 1) return arr[n/2];
    return (arr[n/2 - 1] + arr[n/2]) / 2.0;
}

double mad_sorted(double* arr, int n) {
    // MAD = median(|x - median(x)|)
    double med = median_sorted(arr, n);
    double* dev = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) dev[i] = fabs(arr[i] - med);
    qsort(dev, n, sizeof(double), cmp_double);
    double m = median_sorted(dev, n);
    free(dev);
    return m;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: hobs_windower <input.hwe.gz> <output_prefix> <chrom_size> "
                "[--scales 5kb:5000:1000,...] [--mad_n 3.0]\n");
        return 1;
    }

    const char* infile = argv[1];
    const char* prefix = argv[2];
    int chrom_size = atoi(argv[3]);

    // Default scales
    Scale scales[MAX_SCALES];
    int n_scales = 7;
    const char* default_scales[] = {
        "5kb:5000:1000", "10kb:10000:2000", "50kb:50000:10000",
        "100kb:100000:20000", "250kb:250000:50000", "500kb:500000:100000",
        "1Mb:1000000:200000"
    };
    double mad_n = 3.0;

    // Parse optional args
    for (int i = 4; i < argc; i++) {
        if (!strcmp(argv[i], "--scales") && i+1 < argc) {
            // Parse comma-separated label:win:step
            n_scales = 0;
            char* tok = strtok(argv[++i], ",");
            while (tok && n_scales < MAX_SCALES) {
                sscanf(tok, "%[^:]:%d:%d", scales[n_scales].label,
                       &scales[n_scales].win_bp, &scales[n_scales].step_bp);
                n_scales++;
                tok = strtok(NULL, ",");
            }
        } else if (!strcmp(argv[i], "--mad_n") && i+1 < argc) {
            mad_n = atof(argv[++i]);
        }
    }

    // Apply defaults if not overridden
    if (scales[0].win_bp == 0) {
        for (int i = 0; i < 7; i++) {
            sscanf(default_scales[i], "%[^:]:%d:%d", scales[i].label,
                   &scales[i].win_bp, &scales[i].step_bp);
        }
    }

    // ── Load sites from .hwe.gz ──
    Site* sites = (Site*)malloc(MAX_SITES * sizeof(Site));
    int n_sites = 0;
    char chrom_name[256] = "";

    gzFile gz = gzopen(infile, "rb");
    if (!gz) { fprintf(stderr, "Cannot open %s\n", infile); return 1; }

    char line[4096];
    // Skip header
    gzgets(gz, line, sizeof(line));

    while (gzgets(gz, line, sizeof(line))) {
        if (n_sites >= MAX_SITES) break;

        char chr[128];
        int pos;
        char major, minor;
        double hweFreq, freq, F_val, lrt, pval;

        // Try parsing with hetFreq (10 cols) or without (9 cols)
        double hetFreq = -1;
        int nf = sscanf(line, "%s %d %c %c %lf %lf %lf %lf %lf %lf",
                         chr, &pos, &major, &minor,
                         &hweFreq, &freq, &F_val, &lrt, &pval, &hetFreq);
        if (nf < 9) continue;

        if (chrom_name[0] == '\0') strncpy(chrom_name, chr, sizeof(chrom_name)-1);

        Site* s = &sites[n_sites];
        s->pos = pos;
        s->hweFreq = hweFreq;
        s->freq = freq;
        s->F = F_val;
        s->pval = pval;

        // Mérot logic: Hexp = 2*p*(1-p), Hobs = Hexp - F*Hexp = Hexp*(1-F)
        s->Hexp = 2.0 * hweFreq * (1.0 - hweFreq);
        s->Hobs = s->Hexp * (1.0 - F_val);

        // Clamp to [0, 1]
        if (s->Hobs < 0) s->Hobs = 0;
        if (s->Hobs > 1) s->Hobs = 1;

        n_sites++;
    }
    gzclose(gz);

    fprintf(stderr, "[hobs_windower] %d sites from %s (chrom=%s)\n",
            n_sites, infile, chrom_name);

    if (n_sites == 0) {
        free(sites);
        return 0;
    }

    // ── Write site-level output ──
    char site_path[512];
    snprintf(site_path, sizeof(site_path), "%s.sites.tsv", prefix);
    FILE* fsite = fopen(site_path, "w");
    fprintf(fsite, "chrom\tposition\thweFreq\tfreq\tF\tHexp\tHobs\tpval\n");
    for (int i = 0; i < n_sites; i++) {
        Site* s = &sites[i];
        fprintf(fsite, "%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6e\n",
                chrom_name, s->pos, s->hweFreq, s->freq, s->F,
                s->Hexp, s->Hobs, s->pval);
    }
    fclose(fsite);

    // ── Compute genome-wide (subset-level) robust thresholds ──
    double* all_Hobs = (double*)malloc(n_sites * sizeof(double));
    double* all_F = (double*)malloc(n_sites * sizeof(double));
    for (int i = 0; i < n_sites; i++) {
        all_Hobs[i] = sites[i].Hobs;
        all_F[i] = sites[i].F;
    }
    qsort(all_Hobs, n_sites, sizeof(double), cmp_double);
    qsort(all_F, n_sites, sizeof(double), cmp_double);

    double med_Hobs = median_sorted(all_Hobs, n_sites);
    double mad_Hobs = mad_sorted(all_Hobs, n_sites) * 1.4826; // scale to SD
    double med_F = median_sorted(all_F, n_sites);
    double mad_F = mad_sorted(all_F, n_sites) * 1.4826;

    double lo_Hobs = med_Hobs - mad_n * mad_Hobs;
    double hi_Hobs = med_Hobs + mad_n * mad_Hobs;
    double lo_F = med_F - mad_n * mad_F;
    double hi_F = med_F + mad_n * mad_F;

    fprintf(stderr, "[hobs_windower] Thresholds: Hobs [%.4f, %.4f] (med=%.4f MAD=%.4f) "
            "F [%.4f, %.4f] (med=%.4f MAD=%.4f)\n",
            lo_Hobs, hi_Hobs, med_Hobs, mad_Hobs,
            lo_F, hi_F, med_F, mad_F);

    free(all_Hobs);
    free(all_F);

    // ── Windowed summaries per scale ──
    // Pre-allocate workspace
    double* buf_Hobs = (double*)malloc(n_sites * sizeof(double));
    double* buf_F = (double*)malloc(n_sites * sizeof(double));

    for (int sc = 0; sc < n_scales; sc++) {
        int win = scales[sc].win_bp;
        int step = scales[sc].step_bp;

        char out_path[512];
        snprintf(out_path, sizeof(out_path), "%s.win%s.tsv", prefix, scales[sc].label);
        FILE* fout = fopen(out_path, "w");

        fprintf(fout, "chrom\twindow_start\twindow_end\twindow_center\tn_sites\t"
                "mean_Hobs\tmedian_Hobs\tsd_Hobs\tmin_Hobs\tmax_Hobs\t"
                "mean_F\tmedian_F\tsd_F\tmin_F\tmax_F\t"
                "n_low_Hobs_outlier\tn_high_Hobs_outlier\t"
                "n_low_F_outlier\tn_high_F_outlier\t"
                "frac_low_Hobs_outlier\tfrac_high_Hobs_outlier\t"
                "frac_low_F_outlier\tfrac_high_F_outlier\n");

        int si = 0; // site index pointer
        for (int wstart = 1; wstart <= chrom_size; wstart += step) {
            int wend = wstart + win - 1;
            if (wend > chrom_size) wend = chrom_size;
            int wcenter = (wstart + wend) / 2;

            // Advance si to first site in window
            while (si > 0 && sites[si-1].pos >= wstart) si--;
            while (si < n_sites && sites[si].pos < wstart) si++;

            // Collect sites in window
            int n = 0;
            for (int j = si; j < n_sites && sites[j].pos <= wend; j++) {
                buf_Hobs[n] = sites[j].Hobs;
                buf_F[n] = sites[j].F;
                n++;
            }

            if (n == 0) continue;

            // Sort for median
            qsort(buf_Hobs, n, sizeof(double), cmp_double);
            qsort(buf_F, n, sizeof(double), cmp_double);

            // Compute stats
            double sum_h = 0, sum_f = 0, ssq_h = 0, ssq_f = 0;
            int n_lo_h = 0, n_hi_h = 0, n_lo_f = 0, n_hi_f = 0;

            for (int j = si; j < n_sites && sites[j].pos <= wend; j++) {
                double h = sites[j].Hobs, f = sites[j].F;
                sum_h += h; sum_f += f;
                ssq_h += h*h; ssq_f += f*f;
                if (h < lo_Hobs) n_lo_h++;
                if (h > hi_Hobs) n_hi_h++;
                if (f < lo_F) n_lo_f++;
                if (f > hi_F) n_hi_f++;
            }

            double mean_h = sum_h / n, mean_f = sum_f / n;
            double sd_h = sqrt(fmax(0, ssq_h/n - mean_h*mean_h));
            double sd_f = sqrt(fmax(0, ssq_f/n - mean_f*mean_f));

            fprintf(fout, "%s\t%d\t%d\t%d\t%d\t"
                    "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t"
                    "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t"
                    "%d\t%d\t%d\t%d\t"
                    "%.4f\t%.4f\t%.4f\t%.4f\n",
                    chrom_name, wstart, wend, wcenter, n,
                    mean_h, median_sorted(buf_Hobs, n), sd_h,
                    buf_Hobs[0], buf_Hobs[n-1],
                    mean_f, median_sorted(buf_F, n), sd_f,
                    buf_F[0], buf_F[n-1],
                    n_lo_h, n_hi_h, n_lo_f, n_hi_f,
                    (double)n_lo_h/n, (double)n_hi_h/n,
                    (double)n_lo_f/n, (double)n_hi_f/n);
        }
        fclose(fout);
        fprintf(stderr, "[hobs_windower] Wrote %s\n", out_path);
    }

    free(buf_Hobs);
    free(buf_F);
    free(sites);
    return 0;
}
