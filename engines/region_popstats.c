// =============================================================================
// region_popstats.c — Fast popstats from BEAGLE dosage for inversion confirmation
//
// Single binary computing Hudson Fst, dXY, dA, theta_pi, theta_W, Tajima's D
// from BEAGLE genotype likelihoods (expected dosage). Designed for rapid
// inversion confirmation: "does Fst/dXY spike in this region vs flanks?"
//
// NOT for publication-quality between-species Fst (use ANGSD realSFS for that).
// This is the correct estimator for within-population inversion genotype
// contrasts (Bhatia et al. 2013: Hudson's Fst unbiased for unequal N).
//
// Input:
//   --beagle <file.beagle.gz>     BEAGLE GL file
//   --sample_list <file>          Sample IDs (BAM list order)
//   --groups <g1:file,g2:file>    Group sample lists (e.g. HOM_REF:ref.txt,HET:het.txt)
//   --windows <file>              BED-like: chrom start end [window_id]
//   OR --fixed_win <win_bp:step_bp>
//   --chr <chr>                   Filter chromosome
//
// Output: one row per window with all stats.
//
// Compile: gcc -O3 -march=native -fopenmp -o region_popstats region_popstats.c -lz -lm
//
// Citation:
//   Hudson Fst: Hudson RR et al. (1992) Genetics 132:583
//   Bhatia G et al. (2013) Genome Research 23:1514
//   Watterson GA (1975) Theor Pop Biol 7:256
//   Tajima F (1989) Genetics 123:585
//   MI (nspope PR#208): H(p_pooled) - weighted_mean(H(p_group)), normalized by H_joint
// =============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <zlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_IND 500
#define MAX_GROUPS 10
#define MAX_SITES 2000000

#define SCHEMA_VERSION "region_popstats_v2"

// ── Data ──

typedef struct {
    int pos;
    char allele1, allele2;        // ANGSD numeric allele codes 0..3 (A,C,G,T)
    double dos[MAX_IND];          // expected dosage per individual
} SiteData;

typedef struct {
    char name[64];
    int indices[MAX_IND];         // indices into the full sample array
    int n;
} Group;

typedef struct {
    char id[128];
    int start, end;
} Window;

// Karyotype state per sample (-1 = unknown, 0 = HOM_A, 1 = HET, 2 = HOM_B).
typedef struct {
    int state[MAX_IND];           // per-sample arrangement state
    int n;                        // total samples set
    int n_hom_a, n_het, n_hom_b;  // hard counts across all loaded samples
} Karyotype;

// Per-variant emission context (passed through load_beagle_dosage).
typedef struct {
    FILE* fp;                     // output stream (or NULL)
    Group* groups;
    int n_groups;
    int pol_g1, pol_het, pol_g2;  // group indices for polarisation; -1 if N/A
    int pol_emit;                 // whether to emit polarisation cols
    double peak_threshold;        // default 0.85
    double peak_pass_frac;        // default 0.80
    const char* chrom_label;      // chrom label for output rows
} PerVariantCtx;

static const char* allele_code_to_str(int c) {
    switch (c) {
        case 0: return "A";
        case 1: return "C";
        case 2: return "G";
        case 3: return "T";
        default: return "N";
    }
}

// ── GZ reader ──

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
        if (g->lo_len + n > (int)sizeof(g->lo) - 1) {
            n = sizeof(g->lo) - 1 - g->lo_len;
        }
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

// ── Load group file ──

int load_group(const char* path, Group* g, char all_names[][64], int n_all) {
    FILE* f = fopen(path, "r");
    if (!f) return 0;
    g->n = 0;
    char line[256];
    while (fgets(line, sizeof(line), f) && g->n < MAX_IND) {
        int len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) len--;
        line[len] = 0;
        if (len == 0) continue;
        char* p = line;
        char* slash = strrchr(p, '/');
        if (slash) p = slash + 1;
        char tmp[256];
        int cplen = strlen(p);
        if (cplen > 255) cplen = 255;
        memcpy(tmp, p, cplen); tmp[cplen] = 0;
        char* ext;
        if ((ext = strstr(tmp, ".sorted.markdup.bam"))) *ext = 0;
        else if ((ext = strstr(tmp, ".bam"))) *ext = 0;

        for (int i = 0; i < n_all; i++) {
            if (!strcmp(all_names[i], tmp)) {
                g->indices[g->n++] = i;
                break;
            }
        }
    }
    fclose(f);
    return g->n;
}

// ── Load karyotype: TSV "sample_id<TAB>state" where state ∈ {HOM_A,HET,HOM_B,0,1,2} ──

static int parse_karyo_state(const char* s) {
    if (!strcasecmp(s, "HOM_A") || !strcasecmp(s, "HOM_STD") || !strcmp(s, "0")) return 0;
    if (!strcasecmp(s, "HET") || !strcmp(s, "1")) return 1;
    if (!strcasecmp(s, "HOM_B") || !strcasecmp(s, "HOM_INV") || !strcmp(s, "2")) return 2;
    return -1;
}

int load_karyotype(const char* path, Karyotype* k,
                   char all_names[][64], int n_all) {
    for (int i = 0; i < MAX_IND; i++) k->state[i] = -1;
    k->n = k->n_hom_a = k->n_het = k->n_hom_b = 0;
    FILE* f = fopen(path, "r");
    if (!f) return 0;
    char line[512];
    int matched = 0;
    while (fgets(line, sizeof(line), f)) {
        int len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = 0;
        if (len == 0 || line[0] == '#') continue;
        char sid[256] = "", st[64] = "";
        if (sscanf(line, "%255s %63s", sid, st) < 2) continue;
        char* slash = strrchr(sid, '/');
        char* nm = slash ? slash + 1 : sid;
        char tmp[256]; strncpy(tmp, nm, 255); tmp[255] = 0;
        char* ext;
        if ((ext = strstr(tmp, ".sorted.markdup.bam"))) *ext = 0;
        else if ((ext = strstr(tmp, ".bam"))) *ext = 0;
        int s = parse_karyo_state(st);
        if (s < 0) continue;
        for (int i = 0; i < n_all; i++) {
            if (!strcmp(all_names[i], tmp)) {
                k->state[i] = s;
                matched++;
                if (s == 0) k->n_hom_a++;
                else if (s == 1) k->n_het++;
                else if (s == 2) k->n_hom_b++;
                break;
            }
        }
    }
    fclose(f);
    k->n = matched;
    return matched;
}

// HWE statistics for one regime (= one window) under a fixed karyotype.
// All counts are hard counts of samples with known karyotype state.
typedef struct {
    int n;                        // samples with assigned state
    double Hobs;                  // n_het / n
    double Hexp;                  // 2 * p_A * p_B
    double FIS;                   // 1 - Hobs / Hexp
    double pvalue;                // chi-square HWE
    const char* label;            // HET_excess / HWE_like / HET_deficit
    int n_hom_a, n_het, n_hom_b;
} HWEStats;

static double chi_square_p1df(double chi2) {
    // Survival function of chi^2 with 1 df = erfc(sqrt(chi2/2))
    if (chi2 <= 0) return 1.0;
    return erfc(sqrt(chi2 / 2.0));
}

static const char* fis_label(double FIS, double thr) {
    if (FIS <= -thr) return "HET_excess";
    if (FIS >=  thr) return "HET_deficit";
    return "HWE_like";
}

static void compute_hwe(const Karyotype* k, double thr, HWEStats* out) {
    out->n_hom_a = k->n_hom_a;
    out->n_het = k->n_het;
    out->n_hom_b = k->n_hom_b;
    int N = k->n_hom_a + k->n_het + k->n_hom_b;
    out->n = N;
    if (N < 2) {
        out->Hobs = out->Hexp = out->FIS = 0;
        out->pvalue = 1.0;
        out->label = "HWE_like";
        return;
    }
    double p_A = (2.0 * k->n_hom_a + k->n_het) / (2.0 * N);
    double p_B = 1.0 - p_A;
    out->Hobs = (double)k->n_het / N;
    out->Hexp = 2.0 * p_A * p_B;
    if (out->Hexp > 1e-15) out->FIS = 1.0 - out->Hobs / out->Hexp;
    else out->FIS = 0;
    double e_aa = N * p_A * p_A;
    double e_ab = 2.0 * N * p_A * p_B;
    double e_bb = N * p_B * p_B;
    double chi2 = 0;
    if (e_aa > 1e-9) chi2 += (k->n_hom_a - e_aa) * (k->n_hom_a - e_aa) / e_aa;
    if (e_ab > 1e-9) chi2 += (k->n_het   - e_ab) * (k->n_het   - e_ab) / e_ab;
    if (e_bb > 1e-9) chi2 += (k->n_hom_b - e_bb) * (k->n_hom_b - e_bb) / e_bb;
    out->pvalue = chi_square_p1df(chi2);
    out->label = fis_label(out->FIS, thr);
}

// ── Per-variant row emission ──

static void emit_per_variant_row(PerVariantCtx* pv,
                                 const char* chrom, int pos,
                                 int allele1, int allele2,
                                 const double* per_sample_dos,
                                 const double* per_sample_peak,
                                 int n_ind) {
    if (!pv || !pv->fp) return;
    double g_E[MAX_GROUPS] = {0}, g_peak_sum[MAX_GROUPS] = {0};
    int g_peak_pass[MAX_GROUPS] = {0};
    for (int g = 0; g < pv->n_groups; g++) {
        int ng = pv->groups[g].n;
        for (int k = 0; k < ng; k++) {
            int s = pv->groups[g].indices[k];
            if (s < 0 || s >= n_ind) continue;
            g_E[g] += per_sample_dos[s];
            g_peak_sum[g] += per_sample_peak[s];
            if (per_sample_peak[s] >= pv->peak_threshold) g_peak_pass[g]++;
        }
    }
    for (int g = 0; g < pv->n_groups; g++) {
        int ng = pv->groups[g].n;
        if (ng <= 0) continue;
        double E_freq = g_E[g] / (2.0 * ng);
        double mean_peak = g_peak_sum[g] / ng;
        fprintf(pv->fp, "%s\t%d\t%s\t%s\t%s\t%d\t%.6f\t%.6f\t%.6f",
                chrom, pos,
                allele_code_to_str(allele1), allele_code_to_str(allele2),
                pv->groups[g].name, ng,
                g_E[g], E_freq, mean_peak);
        if (pv->pol_emit) {
            double f1 = (pv->groups[pv->pol_g1].n > 0)
                ? g_E[pv->pol_g1] / (2.0 * pv->groups[pv->pol_g1].n) : 0;
            double f2 = (pv->groups[pv->pol_g2].n > 0)
                ? g_E[pv->pol_g2] / (2.0 * pv->groups[pv->pol_g2].n) : 0;
            double pol = fabs(f1 - f2);
            int het_intermediate = -1;
            if (pv->pol_het >= 0 && pv->groups[pv->pol_het].n > 0) {
                double fh = g_E[pv->pol_het] / (2.0 * pv->groups[pv->pol_het].n);
                double lo = (f1 < f2) ? f1 : f2;
                double hi = (f1 < f2) ? f2 : f1;
                het_intermediate = (fh >= lo && fh <= hi) ? 1 : 0;
            }
            int all_pass = 1;
            for (int gp = 0; gp < pv->n_groups; gp++) {
                if (gp != pv->pol_g1 && gp != pv->pol_g2 && gp != pv->pol_het) continue;
                int ngp = pv->groups[gp].n;
                if (ngp <= 0) { all_pass = 0; break; }
                if ((double)g_peak_pass[gp] / ngp < pv->peak_pass_frac) { all_pass = 0; break; }
            }
            fprintf(pv->fp, "\t%.6f\t%s\t%d", pol,
                    (het_intermediate < 0) ? "NA" : (het_intermediate ? "1" : "0"),
                    all_pass);
        }
        fputc('\n', pv->fp);
    }
}

// ── Load BEAGLE dosage for one chromosome ──

int load_beagle_dosage(const char* path, const char* filter_chr,
                        SiteData* sites, int n_ind_expected,
                        PerVariantCtx* pv) {
    Gz gz;
    if (!gz_open(&gz, path)) { fprintf(stderr, "Cannot open %s\n", path); return 0; }

    char line[1<<20];
    gz_getline(&gz, line, sizeof(line)); // skip header

    // Scratch buffer for per-variant emission (peak per sample).
    double* sample_peak = (pv && pv->fp) ? (double*)malloc(n_ind_expected * sizeof(double)) : NULL;

    int n = 0;
    while (gz_getline(&gz, line, sizeof(line)) && n < MAX_SITES) {
        char* tab1 = strchr(line, '\t');
        if (!tab1) continue;
        *tab1 = 0;
        char* us = strrchr(line, '_');
        if (!us) continue;
        *us = 0;
        char* chr = line;
        int pos = atoi(us + 1);

        if (filter_chr && strcmp(chr, filter_chr) != 0) { *tab1 = '\t'; continue; }

        char* p = tab1 + 1;
        char* tab2 = strchr(p, '\t'); if (!tab2) continue;
        int a1_code = atoi(p);
        char* p2 = tab2 + 1;
        char* tab3 = strchr(p2, '\t'); if (!tab3) continue;
        int a2_code = atoi(p2);
        p = tab3 + 1;

        sites[n].pos = pos;
        sites[n].allele1 = (char)a1_code;
        sites[n].allele2 = (char)a2_code;

        for (int i = 0; i < n_ind_expected; i++) {
            double gl0 = 0, gl1 = 0, gl2 = 0;
            char* endp;

            gl0 = strtod(p, &endp); p = endp;
            while (*p == '\t' || *p == ' ') p++;
            gl1 = strtod(p, &endp); p = endp;
            while (*p == '\t' || *p == ' ') p++;
            gl2 = strtod(p, &endp); p = endp;
            while (*p == '\t' || *p == ' ') p++;

            double s = gl0 + gl1 + gl2;
            if (s > 1e-15) {
                gl0 /= s; gl1 /= s; gl2 /= s;
            }
            sites[n].dos[i] = gl1 + 2.0 * gl2;
            if (sample_peak) {
                double pk = gl0;
                if (gl1 > pk) pk = gl1;
                if (gl2 > pk) pk = gl2;
                sample_peak[i] = pk;
            }
        }

        if (pv && pv->fp) {
            emit_per_variant_row(pv, chr, pos, a1_code, a2_code,
                                 sites[n].dos, sample_peak, n_ind_expected);
        }

        n++;
    }
    if (sample_peak) free(sample_peak);
    gz_close(&gz);
    fprintf(stderr, "[popstats] Loaded %d sites for %s (%d samples)\n",
            n, filter_chr ? filter_chr : "all", n_ind_expected);
    return n;
}

// ── Compute all stats for a window ──

typedef struct {
    char id[128];
    int start, end, n_sites;
    double theta_pi[MAX_GROUPS];
    double hudson_fst[MAX_GROUPS * MAX_GROUPS];
    double dxy[MAX_GROUPS * MAX_GROUPS];
    double dA[MAX_GROUPS * MAX_GROUPS];
    double MI[MAX_GROUPS * MAX_GROUPS];
    double MI_norm[MAX_GROUPS * MAX_GROUPS];
    double theta_pi_all, theta_W_all, tajima_D;
    double Hp;
    int S;
} WinStats;

void compute_window(SiteData* sites, int n_sites, int n_ind,
                    Group* groups, int n_groups,
                    WinStats* out) {
    out->n_sites = n_sites;
    if (n_sites < 5) return;

    // ── Whole-sample theta_pi and theta_W ──
    double sum_pi = 0;
    int S = 0;
    double a1 = 0;
    for (int k = 1; k < n_ind; k++) a1 += 1.0 / k;

    double Hp_num = 0, Hp_den = 0;

    for (int j = 0; j < n_sites; j++) {
        double sum_dos = 0;
        int n_ok = 0;
        for (int i = 0; i < n_ind; i++) {
            sum_dos += sites[j].dos[i];
            n_ok++;
        }
        if (n_ok < 5) continue;
        double p = sum_dos / (2.0 * n_ok);
        double pi_site = 2.0 * p * (1.0 - p) * n_ok / (n_ok - 1.0);
        sum_pi += pi_site;
        if (p > 0.01 && p < 0.99) S++;

        double n_alleles = 2.0 * n_ok;
        Hp_num += 2.0 * p * (1.0 - p) * n_alleles;
        Hp_den += n_alleles;
    }
    out->theta_pi_all = sum_pi / n_sites;
    out->theta_W_all = (a1 > 0) ? (double)S / (a1 * n_sites) : 0;
    out->S = S;
    out->Hp = (Hp_den > 0) ? Hp_num / Hp_den : 0;

    // Tajima's D
    if (S > 1 && n_ind > 2) {
        double a2 = 0;
        for (int k = 1; k < n_ind; k++) a2 += 1.0 / ((double)k * k);
        double b1 = (n_ind + 1.0) / (3.0 * (n_ind - 1.0));
        double b2 = 2.0 * (n_ind * n_ind + n_ind + 3.0) / (9.0 * n_ind * (n_ind - 1.0));
        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (n_ind + 2.0) / (a1 * n_ind) + a2 / (a1 * a1);
        double e1 = c1 / a1;
        double e2 = c2 / (a1 * a1 + a2);
        double D_num = sum_pi - (double)S / a1;
        double D_den = sqrt(e1 * S + e2 * S * (S - 1.0));
        out->tajima_D = (D_den > 1e-15) ? D_num / D_den : 0;
    } else {
        out->tajima_D = 0;
    }

    // ── Per-group and pairwise stats ──
    // Dynamic alloc for per-group allele frequencies per site
    double** pg = (double**)malloc(n_groups * sizeof(double*));
    for (int g = 0; g < n_groups; g++) pg[g] = (double*)malloc(n_sites * sizeof(double));

    for (int j = 0; j < n_sites; j++) {
        for (int g = 0; g < n_groups; g++) {
            double s = 0;
            int ng = groups[g].n;
            for (int k = 0; k < ng; k++) {
                s += sites[j].dos[groups[g].indices[k]];
            }
            pg[g][j] = s / (2.0 * ng);
        }
    }

    // Per-group theta_pi
    for (int g = 0; g < n_groups; g++) {
        int ng = groups[g].n;
        if (ng < 2) { out->theta_pi[g] = 0; continue; }
        double s = 0;
        for (int j = 0; j < n_sites; j++) {
            s += 2.0 * pg[g][j] * (1.0 - pg[g][j]) * ng / (ng - 1.0);
        }
        out->theta_pi[g] = s / n_sites;
    }

    // Pairwise Hudson Fst + dXY + dA + MI (nspope PR#208)
    for (int a = 0; a < n_groups; a++) {
        for (int b = a + 1; b < n_groups; b++) {
            int idx = a * MAX_GROUPS + b;
            int na = groups[a].n, nb = groups[b].n;
            if (na < 2 || nb < 2) {
                out->hudson_fst[idx] = 0; out->dxy[idx] = 0; out->dA[idx] = 0;
                out->MI[idx] = 0; out->MI_norm[idx] = 0;
                continue;
            }

            double sum_num = 0, sum_den = 0, sum_dxy = 0;
            double sum_mi = 0, sum_joint_H = 0;

            for (int j = 0; j < n_sites; j++) {
                double p1 = pg[a][j], p2 = pg[b][j];

                // Hudson Fst numerator/denominator (Bhatia 2013 unbiased)
                double num = (p1 - p2) * (p1 - p2) - p1*(1-p1)/(na-1.0) - p2*(1-p2)/(nb-1.0);
                double den = p1*(1-p2) + p2*(1-p1);
                sum_num += num;
                sum_den += den;

                // dXY
                sum_dxy += p1*(1-p2) + p2*(1-p1);

                // MI (nspope PR#208, Jaquemet & Aguilar 2017)
                double w1 = (double)na / (na + nb);
                double w2 = (double)nb / (na + nb);
                double p_pool = w1 * p1 + w2 * p2;

                double H_pool = 0, H1 = 0, H2 = 0;
                if (p_pool > 1e-10 && p_pool < 1-1e-10)
                    H_pool = -p_pool * log(p_pool) - (1-p_pool) * log(1-p_pool);
                if (p1 > 1e-10 && p1 < 1-1e-10)
                    H1 = -p1 * log(p1) - (1-p1) * log(1-p1);
                if (p2 > 1e-10 && p2 < 1-1e-10)
                    H2 = -p2 * log(p2) - (1-p2) * log(1-p2);

                double mi_site = H_pool - (w1 * H1 + w2 * H2);
                if (mi_site < 0) mi_site = 0;
                sum_mi += mi_site;
                sum_joint_H += H_pool;
            }

            out->hudson_fst[idx] = (sum_den > 1e-15) ? fmax(0, sum_num / sum_den) : 0;
            out->dxy[idx] = sum_dxy / n_sites;
            double within_pi = (out->theta_pi[a] + out->theta_pi[b]) / 2.0;
            out->dA[idx] = out->dxy[idx] - within_pi;
            out->MI[idx] = sum_mi / n_sites;
            out->MI_norm[idx] = (sum_joint_H > 1e-15) ? sum_mi / sum_joint_H : 0;
        }
    }

    for (int g = 0; g < n_groups; g++) free(pg[g]);
    free(pg);
}

// ── Main ──

void print_usage() {
    fprintf(stderr,
        "region_popstats — Fast Hudson Fst / dXY / theta from BEAGLE dosage\n\n"
        "  --beagle <f>          BEAGLE GL file (.beagle.gz)\n"
        "  --sample_list <f>     Sample IDs (BAM list order)\n"
        "  --groups <spec>       Group files: g1:file1,g2:file2[,g3:file3]\n"
        "  --windows <f>         Window coordinates: chrom start end [id]\n"
        "  --fixed_win W:S       Generate fixed bp windows\n"
        "  --chr <n>             Filter chromosome\n"
        "  --out <f>             Output file (default: stdout)\n"
        "  --ncores <n>          Threads (default: 1)\n"
        "\n"
        "Resolution control:\n"
        "  --downsample N        Use every Nth site (1=all, 5=5x, 10=10x, etc)\n"
        "\n"
        "Window anchoring (ANGSD-compatible -type codes):\n"
        "  --type 0              First window where data fills complete window\n"
        "                        (same centers across datasets — ANGSD default)\n"
        "  --type 1              Start at first data position\n"
        "  --type 2              Fixed grid from chromosome position 0\n"
        "\n"
        "Range restriction:\n"
        "  --range START:END     Only process sites in this genomic range (bp)\n"
        "  --site_idx FROM:TO    Only use site indices FROM..TO (0-based)\n"
        "\n"
        "HWE / inversion-context (per-window, treats each window as one regime):\n"
        "  --karyotype <f>             TSV: sample_id<TAB>{HOM_A|HET|HOM_B|0|1|2}\n"
        "                              Adds HWE_Hobs/HWE_Hexp/HWE_FIS/HWE_pvalue/\n"
        "                              HWE_HET_excess_deficit columns to per-window output.\n"
        "  --hwe_fis_threshold F       Threshold for HET_excess/HWE_like/HET_deficit (default 0.10)\n"
        "  --output_naming disambiguated\n"
        "                              Add arrangement_FST_like column = Fst between\n"
        "                              HOM_A and HOM_B (or HOM_STD / HOM_INV) groups.\n"
        "                              Adds naming_convention_version=disambiguated_v1 to header.\n"
        "\n"
        "Per-variant per-group dosage (POD burden, polarisation):\n"
        "  --emit_per_variant_group_dosage <f.tsv>\n"
        "                              Emit one row per (variant × group): chrom, pos_bp,\n"
        "                              ref, alt, group, n_samples_in_group,\n"
        "                              E_alt_dosage_sum, E_alt_freq, mean_GL_peak.\n"
        "  --polarisation_groups A[,H],B\n"
        "                              When 2 names: emit polarisation = |E_alt_freq(A) - E_alt_freq(B)|\n"
        "                              When 3 names (A,H,B): also emit het_intermediate (0/1) and\n"
        "                              gl_peak_pass (samples with peak >= 0.85 in each group >= 80%%).\n"
        "\n"
        "Input format (currently only BEAGLE GL is supported):\n"
        "  --input_format beagle_gl    (accepted for forward compatibility; default)\n"
        "\n"
        "Examples:\n"
        "  # Precise, exact triangle windows\n"
        "  region_popstats --beagle ... --windows triangles.bed --downsample 1\n"
        "\n"
        "  # Coarse 10x, quick scan 100kb, ANGSD-style type 2 grid\n"
        "  region_popstats --beagle ... --fixed_win 100000:20000 --downsample 10 --type 2\n"
        "\n"
        "  # Per-regime HWE with karyotype, disambiguated naming, per-variant burden\n"
        "  region_popstats --beagle ... --regions regimes.bed --karyotype karyo.tsv \\\n"
        "    --groups HOM_A:hom_a.txt,HET:het.txt,HOM_B:hom_b.txt \\\n"
        "    --output_naming disambiguated \\\n"
        "    --emit_per_variant_group_dosage variants.tsv \\\n"
        "    --polarisation_groups HOM_A,HET,HOM_B\n");
}

int main(int argc, char** argv) {
    const char *beagle_path = NULL, *sample_path = NULL, *groups_spec = NULL;
    const char *win_path = NULL, *out_path = NULL, *filter_chr = NULL;
    const char *karyo_path = NULL;
    const char *per_variant_path = NULL;
    const char *polarisation_groups_spec = NULL;
    const char *output_naming = NULL;
    int fixed_win = 0, fixed_step = 0, ncores = 1;
    int downsample = 1;
    int win_type = 0;
    int range_start = -1, range_end = -1;
    int site_idx_from = -1, site_idx_to = -1;
    double hwe_fis_threshold = 0.10;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--beagle") && i+1<argc)      beagle_path = argv[++i];
        else if (!strcmp(argv[i], "--sample_list") && i+1<argc) sample_path = argv[++i];
        else if (!strcmp(argv[i], "--groups") && i+1<argc)      groups_spec = argv[++i];
        else if (!strcmp(argv[i], "--windows") && i+1<argc)     win_path = argv[++i];
        else if (!strcmp(argv[i], "--regions") && i+1<argc)     win_path = argv[++i];
        else if (!strcmp(argv[i], "--fixed_win") && i+1<argc)   sscanf(argv[++i], "%d:%d", &fixed_win, &fixed_step);
        else if (!strcmp(argv[i], "--chr") && i+1<argc)         filter_chr = argv[++i];
        else if (!strcmp(argv[i], "--out") && i+1<argc)         out_path = argv[++i];
        else if (!strcmp(argv[i], "--ncores") && i+1<argc)      ncores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--downsample") && i+1<argc)  downsample = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--type") && i+1<argc)        win_type = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--range") && i+1<argc)       sscanf(argv[++i], "%d:%d", &range_start, &range_end);
        else if (!strcmp(argv[i], "--site_idx") && i+1<argc)    sscanf(argv[++i], "%d:%d", &site_idx_from, &site_idx_to);
        else if (!strcmp(argv[i], "--karyotype") && i+1<argc)   karyo_path = argv[++i];
        else if (!strcmp(argv[i], "--emit_per_variant_group_dosage") && i+1<argc)
                                                                per_variant_path = argv[++i];
        else if (!strcmp(argv[i], "--polarisation_groups") && i+1<argc)
                                                                polarisation_groups_spec = argv[++i];
        else if (!strcmp(argv[i], "--output_naming") && i+1<argc) output_naming = argv[++i];
        else if (!strcmp(argv[i], "--hwe_fis_threshold") && i+1<argc) hwe_fis_threshold = atof(argv[++i]);
        else if (!strcmp(argv[i], "--input_format") && i+1<argc) ++i; // accepted for forward compat (beagle_gl is the only supported format)
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(); return 0;
        }
    }
    if (downsample < 1) downsample = 1;
    int disambiguated = (output_naming && !strcmp(output_naming, "disambiguated"));

    if (!beagle_path || !sample_path) {
        fprintf(stderr, "Need --beagle and --sample_list\n");
        print_usage(); return 1;
    }

    #ifdef _OPENMP
    omp_set_num_threads(ncores);
    #endif

    static char names[MAX_IND][64];
    int n_ind = load_samples(sample_path, names);
    fprintf(stderr, "[popstats] %d samples\n", n_ind);

    Group groups[MAX_GROUPS];
    int n_groups = 0;

    if (groups_spec) {
        char spec_copy[4096];
        strncpy(spec_copy, groups_spec, sizeof(spec_copy)-1);
        spec_copy[sizeof(spec_copy)-1] = 0;
        char* tok = strtok(spec_copy, ",");
        while (tok && n_groups < MAX_GROUPS) {
            char* colon = strchr(tok, ':');
            if (!colon) { tok = strtok(NULL, ","); continue; }
            *colon = 0;
            strncpy(groups[n_groups].name, tok, 63);
            groups[n_groups].name[63] = 0;
            int ng = load_group(colon + 1, &groups[n_groups], names, n_ind);
            fprintf(stderr, "[popstats] Group %s: %d samples\n", tok, ng);
            n_groups++;
            tok = strtok(NULL, ",");
        }
    }

    // Karyotype (optional, for HWE per-window).
    Karyotype karyo;
    int have_karyo = 0;
    if (karyo_path) {
        int m = load_karyotype(karyo_path, &karyo, names, n_ind);
        if (m > 0) {
            have_karyo = 1;
            fprintf(stderr, "[popstats] Karyotype: %d samples (HOM_A=%d HET=%d HOM_B=%d)\n",
                    karyo.n, karyo.n_hom_a, karyo.n_het, karyo.n_hom_b);
        } else {
            fprintf(stderr, "[popstats] WARN: --karyotype %s loaded 0 samples (ignored)\n", karyo_path);
        }
    }

    // Locate HOM_A/HOM_B group indices for disambiguated FST output.
    int idx_hom_a = -1, idx_hom_b = -1;
    if (disambiguated) {
        for (int g = 0; g < n_groups; g++) {
            if (!strcasecmp(groups[g].name, "HOM_A") || !strcasecmp(groups[g].name, "HOM_STD"))
                idx_hom_a = g;
            else if (!strcasecmp(groups[g].name, "HOM_B") || !strcasecmp(groups[g].name, "HOM_INV"))
                idx_hom_b = g;
        }
        if (idx_hom_a < 0 || idx_hom_b < 0) {
            fprintf(stderr, "[popstats] WARN: --output_naming disambiguated requested but groups "
                    "HOM_A/HOM_STD and HOM_B/HOM_INV not both found in --groups\n");
        }
    }

    // Per-variant emission context.
    PerVariantCtx pv = {0};
    pv.pol_g1 = pv.pol_het = pv.pol_g2 = -1;
    pv.peak_threshold = 0.85;
    pv.peak_pass_frac = 0.80;
    pv.groups = groups;
    pv.n_groups = n_groups;
    pv.chrom_label = filter_chr ? filter_chr : ".";
    if (per_variant_path) {
        pv.fp = fopen(per_variant_path, "w");
        if (!pv.fp) {
            fprintf(stderr, "[popstats] Cannot open --emit_per_variant_group_dosage %s\n",
                    per_variant_path);
            return 1;
        }
        // Resolve polarisation groups (2 names: A,B; 3 names: A,H,B).
        if (polarisation_groups_spec) {
            char buf[256]; strncpy(buf, polarisation_groups_spec, sizeof(buf)-1); buf[sizeof(buf)-1]=0;
            char* names_pol[3] = {NULL,NULL,NULL}; int nn = 0;
            char* tok = strtok(buf, ",");
            while (tok && nn < 3) { names_pol[nn++] = tok; tok = strtok(NULL, ","); }
            int idx[3] = {-1,-1,-1};
            for (int j = 0; j < nn; j++) {
                for (int g = 0; g < n_groups; g++) {
                    if (!strcasecmp(groups[g].name, names_pol[j])) { idx[j] = g; break; }
                }
            }
            if (nn == 2) { pv.pol_g1 = idx[0]; pv.pol_g2 = idx[1]; pv.pol_het = -1; }
            else if (nn == 3) { pv.pol_g1 = idx[0]; pv.pol_het = idx[1]; pv.pol_g2 = idx[2]; }
            if (pv.pol_g1 < 0 || pv.pol_g2 < 0) {
                fprintf(stderr, "[popstats] WARN: --polarisation_groups names not found in --groups; "
                        "polarisation columns disabled\n");
                pv.pol_emit = 0;
            } else {
                pv.pol_emit = 1;
                fprintf(stderr, "[popstats] Polarisation: %s%s%s vs %s\n",
                        groups[pv.pol_g1].name,
                        (pv.pol_het >= 0) ? "," : "",
                        (pv.pol_het >= 0) ? groups[pv.pol_het].name : "",
                        groups[pv.pol_g2].name);
            }
        }
        // Header for per-variant output.
        fprintf(pv.fp, "# schema_version=%s_per_variant_v1\n", SCHEMA_VERSION);
        fprintf(pv.fp, "chrom\tpos_bp\tref\talt\tgroup\tn_samples_in_group\t"
                       "E_alt_dosage_sum\tE_alt_freq\tmean_GL_peak");
        if (pv.pol_emit) fprintf(pv.fp, "\tpolarisation\thet_intermediate\tgl_peak_pass");
        fputc('\n', pv.fp);
    }

    SiteData* sites = (SiteData*)malloc(MAX_SITES * sizeof(SiteData));
    int n_sites_raw = load_beagle_dosage(beagle_path, filter_chr, sites, n_ind,
                                         pv.fp ? &pv : NULL);
    if (pv.fp) { fclose(pv.fp); pv.fp = NULL; }
    if (n_sites_raw == 0) { free(sites); return 0; }

    // Apply downsampling
    int n_sites = 0;
    if (downsample > 1) {
        for (int i = 0; i < n_sites_raw; i++) {
            if (i % downsample == 0) {
                if (n_sites != i) sites[n_sites] = sites[i];
                n_sites++;
            }
        }
        fprintf(stderr, "[popstats] Downsampled %dx: %d → %d sites\n",
                downsample, n_sites_raw, n_sites);
    } else {
        n_sites = n_sites_raw;
    }

    // Apply genomic range filter
    if (range_start >= 0 && range_end >= 0) {
        int j = 0;
        for (int i = 0; i < n_sites; i++) {
            if (sites[i].pos >= range_start && sites[i].pos <= range_end) {
                if (j != i) sites[j] = sites[i];
                j++;
            }
        }
        fprintf(stderr, "[popstats] Range %d-%d: %d → %d sites\n",
                range_start, range_end, n_sites, j);
        n_sites = j;
    }

    // Apply site index filter
    if (site_idx_from >= 0 && site_idx_to >= 0) {
        if (site_idx_to >= n_sites) site_idx_to = n_sites - 1;
        if (site_idx_from < n_sites) {
            int new_n = site_idx_to - site_idx_from + 1;
            if (site_idx_from > 0)
                memmove(sites, sites + site_idx_from, new_n * sizeof(SiteData));
            fprintf(stderr, "[popstats] Site idx %d-%d: %d sites\n",
                    site_idx_from, site_idx_to, new_n);
            n_sites = new_n;
        }
    }

    if (n_sites == 0) { fprintf(stderr, "[popstats] No sites after filtering\n"); free(sites); return 0; }

    // Load or generate windows
    Window* wins = NULL;
    int n_wins = 0;

    if (win_path) {
        FILE* f = fopen(win_path, "r");
        if (!f) { fprintf(stderr, "Cannot open %s\n", win_path); return 1; }
        int cap = 100000;
        wins = (Window*)malloc(cap * sizeof(Window));
        char line[512];
        while (fgets(line, sizeof(line), f) && n_wins < cap) {
            char chr[128], id[128] = "";
            int s, e;
            int nf = sscanf(line, "%s %d %d %127s", chr, &s, &e, id);
            if (nf < 3) continue;
            if (filter_chr && strcmp(chr, filter_chr) != 0) continue;
            if (chr[0] >= '0' && chr[0] <= '9') {
                wins[n_wins].start = atoi(chr);
                wins[n_wins].end = s;
                snprintf(wins[n_wins].id, sizeof(wins[n_wins].id), "w%d", n_wins+1);
            } else {
                wins[n_wins].start = s;
                wins[n_wins].end = e;
                if (id[0]) { strncpy(wins[n_wins].id, id, 127); wins[n_wins].id[127] = 0; }
                else snprintf(wins[n_wins].id, sizeof(wins[n_wins].id), "w%d", n_wins+1);
            }
            n_wins++;
        }
        fclose(f);
    } else if (fixed_win > 0) {
        int max_pos = sites[n_sites-1].pos;
        int first_pos = sites[0].pos;
        int win_start;

        switch (win_type) {
            case 0: win_start = first_pos; break;
            case 1: win_start = first_pos; break;
            case 2: win_start = 0; break;
            default: win_start = 0;
        }

        int cap = ((max_pos - win_start) / fixed_step) + 2;
        wins = (Window*)malloc(cap * sizeof(Window));
        for (int s = win_start; s <= max_pos && n_wins < cap; s += fixed_step) {
            wins[n_wins].start = s;
            wins[n_wins].end = s + fixed_win - 1;
            snprintf(wins[n_wins].id, sizeof(wins[n_wins].id),
                     "%s_w%d", filter_chr ? filter_chr : "chr", n_wins+1);
            n_wins++;
        }
    } else {
        n_wins = 1;
        wins = (Window*)malloc(sizeof(Window));
        wins[0].start = sites[0].pos;
        wins[0].end = sites[n_sites-1].pos;
        snprintf(wins[0].id, sizeof(wins[0].id), "%s_all", filter_chr ? filter_chr : "all");
    }

    fprintf(stderr, "[popstats] %d windows, %d groups, downsample=%dx, type=%d\n",
            n_wins, n_groups, downsample, win_type);

    FILE* fout = out_path ? fopen(out_path, "w") : stdout;

    // Schema version (first line).
    fprintf(fout, "# schema_version=%s\n", SCHEMA_VERSION);
    if (disambiguated) fprintf(fout, "# naming_convention_version=disambiguated_v1\n");

    // Metadata comment
    fprintf(fout, "# region_popstats downsample=%d type=%d n_sites_loaded=%d n_sites_used=%d",
            downsample, win_type, n_sites_raw, n_sites);
    if (range_start >= 0) fprintf(fout, " range=%d:%d", range_start, range_end);
    if (site_idx_from >= 0) fprintf(fout, " site_idx=%d:%d", site_idx_from, site_idx_to);
    if (have_karyo) fprintf(fout, " karyotype_N=%d (HOM_A=%d HET=%d HOM_B=%d) hwe_fis_threshold=%.3f",
                            karyo.n, karyo.n_hom_a, karyo.n_het, karyo.n_hom_b, hwe_fis_threshold);
    fprintf(fout, "\n");

    // Header
    fprintf(fout, "window_id\tchrom\tstart\tend\tn_sites\tn_sites_used\tS\t"
                  "theta_pi_all\ttheta_W_all\tTajima_D\tHp");
    for (int g = 0; g < n_groups; g++)
        fprintf(fout, "\ttheta_pi_%s", groups[g].name);
    for (int a = 0; a < n_groups; a++)
        for (int b = a+1; b < n_groups; b++)
            fprintf(fout, "\tFst_%s_%s\tdXY_%s_%s\tdA_%s_%s\tMI_%s_%s\tMInorm_%s_%s",
                    groups[a].name, groups[b].name,
                    groups[a].name, groups[b].name,
                    groups[a].name, groups[b].name,
                    groups[a].name, groups[b].name,
                    groups[a].name, groups[b].name);
    if (disambiguated) fprintf(fout, "\tarrangement_FST_like");
    if (have_karyo)
        fprintf(fout, "\tHWE_n\tHWE_n_HOM_A\tHWE_n_HET\tHWE_n_HOM_B\t"
                      "HWE_Hobs\tHWE_Hexp\tHWE_FIS\tHWE_pvalue\tHWE_HET_excess_deficit");
    fprintf(fout, "\n");

    // Process windows
    for (int wi = 0; wi < n_wins; wi++) {
        int ws = wins[wi].start, we = wins[wi].end;

        int lo = 0, hi = n_sites;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (sites[mid].pos < ws) lo = mid + 1;
            else hi = mid;
        }
        int first = lo;

        int last = first;
        while (last < n_sites && sites[last].pos <= we) last++;
        int n_win = last - first;

        WinStats ws_out;
        memset(&ws_out, 0, sizeof(ws_out));

        if (n_win >= 5) {
            compute_window(sites + first, n_win, n_ind, groups, n_groups, &ws_out);
        }

        fprintf(fout, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.4f\t%.6f",
                wins[wi].id, filter_chr ? filter_chr : ".",
                ws, we, n_win,
                (downsample > 1) ? n_win : n_win,
                ws_out.S,
                ws_out.theta_pi_all, ws_out.theta_W_all, ws_out.tajima_D,
                ws_out.Hp);

        for (int g = 0; g < n_groups; g++)
            fprintf(fout, "\t%.6f", ws_out.theta_pi[g]);

        for (int a = 0; a < n_groups; a++)
            for (int b = a+1; b < n_groups; b++) {
                int idx = a * MAX_GROUPS + b;
                fprintf(fout, "\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f",
                        ws_out.hudson_fst[idx], ws_out.dxy[idx], ws_out.dA[idx],
                        ws_out.MI[idx], ws_out.MI_norm[idx]);
            }

        if (disambiguated) {
            double arrangement_fst = 0;
            if (idx_hom_a >= 0 && idx_hom_b >= 0 && idx_hom_a != idx_hom_b) {
                int a = (idx_hom_a < idx_hom_b) ? idx_hom_a : idx_hom_b;
                int b = (idx_hom_a < idx_hom_b) ? idx_hom_b : idx_hom_a;
                arrangement_fst = ws_out.hudson_fst[a * MAX_GROUPS + b];
            }
            fprintf(fout, "\t%.6f", arrangement_fst);
        }

        if (have_karyo) {
            HWEStats h;
            compute_hwe(&karyo, hwe_fis_threshold, &h);
            fprintf(fout, "\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6e\t%s",
                    h.n, h.n_hom_a, h.n_het, h.n_hom_b,
                    h.Hobs, h.Hexp, h.FIS, h.pvalue, h.label);
        }

        fprintf(fout, "\n");
    }

    if (out_path) fclose(fout);
    free(sites); free(wins);
    fprintf(stderr, "[popstats] Done: %d windows\n", n_wins);
    return 0;
}
