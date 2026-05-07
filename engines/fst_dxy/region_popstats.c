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
#include <math.h>
#include <zlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_IND 500
#define MAX_GROUPS 10
#define MAX_SITES 2000000

// ── Data ──

typedef struct {
    int pos;
    double dos[MAX_IND];  // expected dosage per individual
} SiteData;

typedef struct {
    char name[64];
    int indices[MAX_IND]; // indices into the full sample array
    int n;
} Group;

typedef struct {
    char id[128];
    int start, end;
} Window;

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

// ── Load BEAGLE dosage for one chromosome ──

int load_beagle_dosage(const char* path, const char* filter_chr,
                        SiteData* sites, int n_ind_expected) {
    Gz gz;
    if (!gz_open(&gz, path)) { fprintf(stderr, "Cannot open %s\n", path); return 0; }

    char line[1<<20];
    gz_getline(&gz, line, sizeof(line)); // skip header

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
        char* tab3 = strchr(tab2 + 1, '\t'); if (!tab3) continue;
        p = tab3 + 1;

        sites[n].pos = pos;

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
        }
        n++;
    }
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
        "Examples:\n"
        "  # Precise, exact triangle windows\n"
        "  region_popstats --beagle ... --windows triangles.bed --downsample 1\n"
        "\n"
        "  # Coarse 10x, quick scan 100kb, ANGSD-style type 2 grid\n"
        "  region_popstats --beagle ... --fixed_win 100000:20000 --downsample 10 --type 2\n"
        "\n"
        "  # Only a 2Mb region around a candidate\n"
        "  region_popstats --beagle ... --fixed_win 50000:10000 --range 35000000:37000000\n");
}

int main(int argc, char** argv) {
    const char *beagle_path = NULL, *sample_path = NULL, *groups_spec = NULL;
    const char *win_path = NULL, *out_path = NULL, *filter_chr = NULL;
    int fixed_win = 0, fixed_step = 0, ncores = 1;
    int downsample = 1;
    int win_type = 0;
    int range_start = -1, range_end = -1;
    int site_idx_from = -1, site_idx_to = -1;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--beagle") && i+1<argc)      beagle_path = argv[++i];
        else if (!strcmp(argv[i], "--sample_list") && i+1<argc) sample_path = argv[++i];
        else if (!strcmp(argv[i], "--groups") && i+1<argc)      groups_spec = argv[++i];
        else if (!strcmp(argv[i], "--windows") && i+1<argc)     win_path = argv[++i];
        else if (!strcmp(argv[i], "--fixed_win") && i+1<argc)   sscanf(argv[++i], "%d:%d", &fixed_win, &fixed_step);
        else if (!strcmp(argv[i], "--chr") && i+1<argc)         filter_chr = argv[++i];
        else if (!strcmp(argv[i], "--out") && i+1<argc)         out_path = argv[++i];
        else if (!strcmp(argv[i], "--ncores") && i+1<argc)      ncores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--downsample") && i+1<argc)  downsample = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--type") && i+1<argc)        win_type = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--range") && i+1<argc)       sscanf(argv[++i], "%d:%d", &range_start, &range_end);
        else if (!strcmp(argv[i], "--site_idx") && i+1<argc)    sscanf(argv[++i], "%d:%d", &site_idx_from, &site_idx_to);
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(); return 0;
        }
    }
    if (downsample < 1) downsample = 1;

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

    SiteData* sites = (SiteData*)malloc(MAX_SITES * sizeof(SiteData));
    int n_sites_raw = load_beagle_dosage(beagle_path, filter_chr, sites, n_ind);
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

    // Metadata comment
    fprintf(fout, "# region_popstats downsample=%d type=%d n_sites_loaded=%d n_sites_used=%d",
            downsample, win_type, n_sites_raw, n_sites);
    if (range_start >= 0) fprintf(fout, " range=%d:%d", range_start, range_end);
    if (site_idx_from >= 0) fprintf(fout, " site_idx=%d:%d", site_idx_from, site_idx_to);
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

        fprintf(fout, "\n");
    }

    if (out_path) fclose(fout);
    free(sites); free(wins);
    fprintf(stderr, "[popstats] Done: %d windows\n", n_wins);
    return 0;
}
