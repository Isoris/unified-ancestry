// =============================================================================
// fdM.c — Patterson's D + Malinsky 2015 fdM introgression statistic, per
//         genomic window, from BEAGLE genotype likelihoods.
//
// ── What Patterson's D measures (ABBA-BABA test, four-taxon) ────────────────
//
// Species tree assumed: ((P1, P2), P3), Outgroup.
//   P1, P2 are sister groups; P3 is more distant; the outgroup defines the
//   ancestral allele (A). Derived = B (anything different from O).
//
// At each biallelic SNP, classify the 4-pop allele state:
//   ABBA  : P1=A  P2=B  P3=B  O=A   (P2 + P3 share derived allele)
//   BABA  : P1=B  P2=A  P3=B  O=A   (P1 + P3 share derived allele)
//
//   D = (Σ ABBA − Σ BABA) / (Σ ABBA + Σ BABA)
//
//   D ≈ 0  : symmetric — incomplete-lineage-sorting only, tree-like.
//   D > 0  : excess ABBA → P3 shares more derived alleles with P2
//                          (consistent with P2↔P3 gene flow).
//   D < 0  : excess BABA → P3 shares more derived alleles with P1
//                          (consistent with P1↔P3 gene flow).
//
// Why it works: under the null (no gene flow), ILS produces ABBA and BABA at
// equal expected rates. Introgression P3↔P2 pushes more sites into ABBA;
// P3↔P1 pushes more into BABA. The imbalance is the signal.
//
// fdM (Malinsky 2015) is the modified-fd that's symmetric around 0 and bounded
// in [-1, 1]: same numerator as D but divided by the maximum-introgression
// expectation, with sign matching D so |fdM| is interpretable as the fraction
// of admixture relative to a "complete introgression" reference.
//
// Catfish-atlas examples:
//   - P1, P2 = two hatchery families; P3 = wild lineage; O = distant species
//     → tests for differential introgression of wild ancestry into families.
//   - Run inside vs outside an inversion candidate to ask whether the
//     inversion carries introgressed ancestry.
//   - Genome-wide jackknife (--jackknife_blocks N) gives the SE / Z / p for
//     "is the genome-wide signal significantly non-zero?"
//
// Mental shortcut: D asks whether P3 shares more derived alleles with P1 or P2.
//                   fdM is the same question, but signed and on [-1, 1].
//
// ── Genealogy summary ───────────────────────────────────────────────────────
//
// Genealogy: ((P1, (P2, P3)), O), where O is the outgroup.
//   fdM > 0  → gene flow P3 → P2
//   fdM < 0  → gene flow P3 → P1
//   fdM ∈ [-1, 1]
//
// Per-site contributions (p_k = expected allele1 frequency in pop k from GLs;
// p_O is treated as the ancestral allele freq, no enforcement of fixation):
//
//   num_i  = p3 · (1 - p_O) · (p2 - p1)             (= ABBA − BABA)
//   ABBA   = (1 - p1) · p2 · p3 · (1 - p_O)
//   BABA   = p1 · (1 - p2) · p3 · (1 - p_O)
//   if (p2 ≥ p1):  denom_i = max(p2, p3) · (1 - p_O) · (max(p2, p3) − p1)
//   else:          denom_i = max(p1, p3) · (1 - p_O) · (p2 − max(p1, p3))
//
//   D_window   = Σ(ABBA - BABA) / Σ(ABBA + BABA)
//   fdM_window = Σ num_i / Σ denom_i
//
// Works for cross-species AND within-species introgression (e.g. P1 = family A,
// P2 = family B, P3 = candidate donor family, O = distant family C).
//
// Biallelic assumption: BEAGLE GLs are inherently biallelic (GL triplet per
// individual). We do NOT impose extra biallelic-only filtering on top — sites
// whose major/minor split is uncertain pass through with whatever frequency
// the GLs produce.
//
// Reference: Malinsky M et al. (2015) Science 350:1493
//            (built on Green 2010 D, Martin 2015 f_d).
//
// Compile: gcc -O3 -march=native -fopenmp -o fdM fdM.c -lz -lm
// =============================================================================

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <zlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_IND 1024
#define MAX_SITES 5000000
#define SCHEMA_VERSION "fdM_v1"

// ── Gz reader (line-buffered) ───────────────────────────────────────────────

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

// ── Sample-list loader (BAM-path-aware, like other engines) ─────────────────

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
        int cplen = len < 255 ? len : 255;
        memcpy(tmp, p, cplen); tmp[cplen] = 0;
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

// ── Pop loader: same format as region_popstats --groups, but exactly 4 pops ──

typedef struct {
    char name[16];
    int  indices[MAX_IND];
    int  n;
} Pop;

static int load_pop(const char* path, Pop* p, char all_names[][64], int n_all) {
    FILE* f = fopen(path, "r");
    if (!f) return 0;
    p->n = 0;
    char line[256];
    while (fgets(line, sizeof(line), f) && p->n < MAX_IND) {
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
        if      ((ext = strstr(tmp, ".sorted.markdup.bam"))) *ext = 0;
        else if ((ext = strstr(tmp, ".bam")))                 *ext = 0;
        for (int i = 0; i < n_all; i++) {
            if (!strcmp(all_names[i], tmp)) { p->indices[p->n++] = i; break; }
        }
    }
    fclose(f);
    return p->n;
}

// Spec like "P1:file1.txt,P2:file2.txt,P3:file3.txt,O:fileO.txt"
static int parse_pops_spec(const char* spec, Pop pops[4], char names[][64], int n_ind) {
    char buf[4096];
    strncpy(buf, spec, sizeof(buf) - 1); buf[sizeof(buf) - 1] = 0;
    int p1 = -1, p2 = -1, p3 = -1, po = -1;
    char* tok = strtok(buf, ",");
    while (tok) {
        char* colon = strchr(tok, ':');
        if (!colon) { tok = strtok(NULL, ","); continue; }
        *colon = 0;
        const char* name = tok;
        const char* fpath = colon + 1;
        int slot = -1;
        if      (!strcasecmp(name, "P1")) slot = 0;
        else if (!strcasecmp(name, "P2")) slot = 1;
        else if (!strcasecmp(name, "P3")) slot = 2;
        else if (!strcasecmp(name, "O") || !strcasecmp(name, "OUT") || !strcasecmp(name, "OUTGROUP")) slot = 3;
        else {
            fprintf(stderr, "[fdM] WARN: unknown pop key '%s' (expected P1/P2/P3/O)\n", name);
            tok = strtok(NULL, ","); continue;
        }
        strncpy(pops[slot].name, name, sizeof(pops[slot].name) - 1);
        pops[slot].name[sizeof(pops[slot].name) - 1] = 0;
        int nl = load_pop(fpath, &pops[slot], names, n_ind);
        fprintf(stderr, "[fdM] Pop %s: %d samples (from %s)\n", name, nl, fpath);
        if      (slot == 0) p1 = nl;
        else if (slot == 1) p2 = nl;
        else if (slot == 2) p3 = nl;
        else                po = nl;
        tok = strtok(NULL, ",");
    }
    if (p1 <= 0 || p2 <= 0 || p3 <= 0 || po <= 0) {
        fprintf(stderr, "[fdM] All four pops (P1, P2, P3, O) must be non-empty.\n");
        return 0;
    }
    return 1;
}

// ── Per-site population frequencies ─────────────────────────────────────────

typedef struct {
    char chrom[32];
    int  pos;
    double p1, p2, p3, po;
} SiteFreqs;

// Compute p_k = sum(dosage) / (2 * N) over samples in pop k.
static inline double pop_freq(const double* dos, const Pop* pop) {
    if (pop->n <= 0) return NAN;
    double s = 0;
    for (int k = 0; k < pop->n; k++) s += dos[pop->indices[k]];
    return s / (2.0 * pop->n);
}

static int load_beagle_freqs(const char* path, const char* filter_chr,
                              SiteFreqs* sites, int n_ind, const Pop pops[4]) {
    Gz gz;
    if (!gz_open(&gz, path)) { fprintf(stderr, "[fdM] Cannot open %s\n", path); return 0; }
    char line[1<<20];
    gz_getline(&gz, line, sizeof(line)); // skip header

    double* dos = (double*)malloc((size_t)n_ind * sizeof(double));
    int n = 0;
    while (gz_getline(&gz, line, sizeof(line)) && n < MAX_SITES) {
        char* tab1 = strchr(line, '\t');
        if (!tab1) continue;
        *tab1 = 0;
        char* us = strrchr(line, '_');
        if (!us) continue;
        *us = 0;
        const char* chr = line;
        int pos = atoi(us + 1);
        if (filter_chr && strcmp(chr, filter_chr) != 0) { *tab1 = '\t'; continue; }

        char* p = tab1 + 1;
        char* tab2 = strchr(p, '\t');       if (!tab2) continue;
        char* tab3 = strchr(tab2 + 1, '\t'); if (!tab3) continue;
        p = tab3 + 1;

        for (int i = 0; i < n_ind; i++) {
            double gl0 = 0, gl1 = 0, gl2 = 0;
            char* endp;
            gl0 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;
            gl1 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;
            gl2 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;
            double s = gl0 + gl1 + gl2;
            if (s > 1e-15) { gl0 /= s; gl1 /= s; gl2 /= s; }
            dos[i] = gl1 + 2.0 * gl2;
        }

        strncpy(sites[n].chrom, chr, sizeof(sites[n].chrom) - 1);
        sites[n].chrom[sizeof(sites[n].chrom) - 1] = 0;
        sites[n].pos = pos;
        sites[n].p1 = pop_freq(dos, &pops[0]);
        sites[n].p2 = pop_freq(dos, &pops[1]);
        sites[n].p3 = pop_freq(dos, &pops[2]);
        sites[n].po = pop_freq(dos, &pops[3]);
        n++;
    }
    free(dos);
    gz_close(&gz);
    fprintf(stderr, "[fdM] Loaded %d sites for %s\n", n, filter_chr ? filter_chr : "all");
    return n;
}

// ── Per-window fdM computation ──────────────────────────────────────────────

typedef struct {
    int n_sites;
    int n_informative;
    double sum_ABBA, sum_BABA;
    double sum_num, sum_denom;
    double D;
    double fdM;
} WindowFdM;

static void compute_fdM_window(const SiteFreqs* sites, int first, int last,
                                WindowFdM* out) {
    memset(out, 0, sizeof(*out));
    out->D = out->fdM = NAN;
    for (int j = first; j < last; j++) {
        double p1 = sites[j].p1, p2 = sites[j].p2, p3 = sites[j].p3, po = sites[j].po;
        if (!isfinite(p1) || !isfinite(p2) || !isfinite(p3) || !isfinite(po)) continue;
        out->n_sites++;

        double q_o = 1.0 - po;
        double ABBA = (1.0 - p1) * p2 * p3 * q_o;
        double BABA = p1 * (1.0 - p2) * p3 * q_o;
        double num  = p3 * q_o * (p2 - p1);

        // Malinsky 2015 fdM: denom is the "maximum-introgression" expectation,
        // kept non-negative so fdM sign comes from num (P3→P2 positive, P3→P1 negative).
        double denom;
        if (p2 >= p1) {
            double pd = (p2 > p3) ? p2 : p3;
            denom = pd * q_o * (pd - p1);
        } else {
            double pd = (p1 > p3) ? p1 : p3;
            denom = pd * q_o * (pd - p2);
        }

        out->sum_ABBA  += ABBA;
        out->sum_BABA  += BABA;
        out->sum_num   += num;
        out->sum_denom += denom;
        if (fabs(num) > 1e-15) out->n_informative++;
    }

    double total = out->sum_ABBA + out->sum_BABA;
    if (total > 0)         out->D   = (out->sum_ABBA - out->sum_BABA) / total;
    if (out->sum_denom != 0) out->fdM = out->sum_num / out->sum_denom;
}

// ── Windows ────────────────────────────────────────────────────────────────

typedef struct { char id[128]; int start, end; } Window;

// ── Main ───────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "fdM — Patterson's D + Malinsky 2015 fdM per window, from BEAGLE GLs.\n"
        "\n"
        "  --beagle <f.beagle.gz>     Required.\n"
        "  --sample_list <f>          Required. Sample IDs (one per line, same order as BEAGLE).\n"
        "  --pops <spec>              Required. P1:file,P2:file,P3:file,O:file\n"
        "                              Genealogy: ((P1, (P2, P3)), O).\n"
        "                              fdM > 0 → P3 → P2; fdM < 0 → P3 → P1.\n"
        "                              Same files work for cross-species or within-species\n"
        "                              (e.g. P1=familyA, P2=familyB, P3=candidate donor, O=distant family).\n"
        "  --chr <name>               Filter to a single chromosome.\n"
        "  --windows <bed>            BED of windows: chrom start end [id].\n"
        "  --fixed_win W:S            Fixed-size windows of W bp, step S.\n"
        "  --jackknife_blocks N       Block jackknife for the GENOME-wide fdM.\n"
        "                              Sites split into N equally-sized blocks.\n"
        "                              Emits an extra 'GENOME' row with jackknife_SE,\n"
        "                              jackknife_Z, jackknife_p. Per-window rows have NA\n"
        "                              in those columns. Typical N=50–100.\n"
        "  --out <f>                  Output TSV (default stdout).\n"
        "\n"
        "Note: BEAGLE GL triplets are inherently biallelic; we don't apply extra\n"
        "      biallelic-only filtering. Multi-allelic loci are not in scope.\n"
        "\n"
        "Output columns (schema_version=" SCHEMA_VERSION "):\n"
        "  window_id chrom start end n_sites n_informative\n"
        "  D fdM sum_ABBA sum_BABA sum_num sum_denom\n"
        "  jackknife_SE jackknife_Z jackknife_p   (NA except in GENOME row)\n");
}

int main(int argc, char** argv) {
    const char *beagle = NULL, *sample_path = NULL, *pops_spec = NULL,
               *filter_chr = NULL, *win_path = NULL, *out_path = NULL;
    int fixed_win = 0, fixed_step = 0;
    int jackknife_blocks = 0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--beagle")           && i+1<argc) beagle = argv[++i];
        else if (!strcmp(argv[i], "--sample_list")      && i+1<argc) sample_path = argv[++i];
        else if (!strcmp(argv[i], "--pops")             && i+1<argc) pops_spec = argv[++i];
        else if (!strcmp(argv[i], "--chr")              && i+1<argc) filter_chr = argv[++i];
        else if (!strcmp(argv[i], "--windows")          && i+1<argc) win_path = argv[++i];
        else if (!strcmp(argv[i], "--fixed_win")        && i+1<argc) sscanf(argv[++i], "%d:%d", &fixed_win, &fixed_step);
        else if (!strcmp(argv[i], "--jackknife_blocks") && i+1<argc) jackknife_blocks = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out")              && i+1<argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[fdM] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }

    if (!beagle || !sample_path || !pops_spec) { print_usage(); return 1; }

    static char names[MAX_IND][64];
    int n_ind = load_samples(sample_path, names);
    fprintf(stderr, "[fdM] %d samples\n", n_ind);
    if (n_ind == 0) return 2;

    Pop pops[4]; memset(pops, 0, sizeof(pops));
    if (!parse_pops_spec(pops_spec, pops, names, n_ind)) return 3;

    SiteFreqs* sites = (SiteFreqs*)malloc((size_t)MAX_SITES * sizeof(SiteFreqs));
    int n_sites = load_beagle_freqs(beagle, filter_chr, sites, n_ind, pops);
    if (n_sites == 0) { free(sites); return 0; }

    // Build window list
    Window* wins = NULL;
    int n_wins = 0;
    if (win_path) {
        FILE* f = fopen(win_path, "r");
        if (!f) { fprintf(stderr, "[fdM] Cannot open --windows %s\n", win_path); free(sites); return 4; }
        int cap = 100000;
        wins = (Window*)malloc(cap * sizeof(Window));
        char line[512];
        while (fgets(line, sizeof(line), f) && n_wins < cap) {
            char chr[64], id[128] = "";
            int s, e;
            int nf = sscanf(line, "%63s %d %d %127s", chr, &s, &e, id);
            if (nf < 3) continue;
            if (filter_chr && strcmp(chr, filter_chr) != 0) continue;
            wins[n_wins].start = s; wins[n_wins].end = e;
            if (id[0]) { strncpy(wins[n_wins].id, id, 127); wins[n_wins].id[127] = 0; }
            else        snprintf(wins[n_wins].id, sizeof(wins[n_wins].id), "w%d", n_wins+1);
            n_wins++;
        }
        fclose(f);
    } else if (fixed_win > 0) {
        int max_pos = sites[n_sites-1].pos;
        int cap = (max_pos / fixed_step) + 2;
        wins = (Window*)malloc(cap * sizeof(Window));
        for (int s = 0; s <= max_pos && n_wins < cap; s += fixed_step) {
            wins[n_wins].start = s;
            wins[n_wins].end = s + fixed_win - 1;
            snprintf(wins[n_wins].id, sizeof(wins[n_wins].id), "w%d", n_wins + 1);
            n_wins++;
        }
    } else {
        n_wins = 1;
        wins = (Window*)malloc(sizeof(Window));
        wins[0].start = sites[0].pos;
        wins[0].end = sites[n_sites-1].pos;
        snprintf(wins[0].id, sizeof(wins[0].id), "%s_all", filter_chr ? filter_chr : "all");
    }

    fprintf(stderr, "[fdM] %d windows over %d sites\n", n_wins, n_sites);

    FILE* fout = out_path ? fopen(out_path, "w") : stdout;
    if (!fout) { fprintf(stderr, "[fdM] Cannot open --out %s\n", out_path); free(sites); free(wins); return 5; }

    fprintf(fout, "# schema_version=%s P1=%s P2=%s P3=%s O=%s",
            SCHEMA_VERSION, pops[0].name, pops[1].name, pops[2].name, pops[3].name);
    if (jackknife_blocks > 0) fprintf(fout, " jackknife_blocks=%d", jackknife_blocks);
    fputc('\n', fout);
    fprintf(fout, "window_id\tchrom\tstart\tend\tn_sites\tn_informative\t"
                  "D\tfdM\tsum_ABBA\tsum_BABA\tsum_num\tsum_denom\t"
                  "jackknife_SE\tjackknife_Z\tjackknife_p\n");

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

        WindowFdM w;
        compute_fdM_window(sites, first, last, &w);

        fprintf(fout, "%s\t%s\t%d\t%d\t%d\t%d",
                wins[wi].id, filter_chr ? filter_chr : ".",
                ws, we, w.n_sites, w.n_informative);
        #define EM(v) do { \
            if (isnan(v))      fprintf(fout, "\tNA"); \
            else if (isinf(v)) fprintf(fout, "\tInf"); \
            else               fprintf(fout, "\t%.6g", (double)(v)); \
        } while (0)
        EM(w.D);
        EM(w.fdM);
        EM(w.sum_ABBA);
        EM(w.sum_BABA);
        EM(w.sum_num);
        EM(w.sum_denom);
        fputs("\tNA\tNA\tNA", fout);  // jackknife cols filled only in GENOME row
        #undef EM
        fputc('\n', fout);
    }

    // Genome-wide block jackknife (delete-1 over equal-sized site blocks).
    if (jackknife_blocks > 0) {
        WindowFdM gw;
        compute_fdM_window(sites, 0, n_sites, &gw);

        int B = jackknife_blocks;
        if (B > n_sites) B = n_sites;
        if (B < 2) {
            fprintf(stderr, "[fdM] WARN: --jackknife_blocks needs ≥ 2 blocks; got %d sites → %d blocks. Skipping.\n", n_sites, B);
        } else {
            int sites_per_block = n_sites / B;
            double* fdM_dropped = (double*)malloc((size_t)B * sizeof(double));
            int n_effective = 0;
            double sum_dropped = 0;

            for (int b = 0; b < B; b++) {
                int first = b * sites_per_block;
                int last  = (b == B - 1) ? n_sites : (b + 1) * sites_per_block;
                WindowFdM blk;
                compute_fdM_window(sites, first, last, &blk);
                double dropped_denom = gw.sum_denom - blk.sum_denom;
                double dropped_fdM = (dropped_denom != 0)
                                     ? (gw.sum_num - blk.sum_num) / dropped_denom : NAN;
                fdM_dropped[b] = dropped_fdM;
                if (isfinite(dropped_fdM)) { sum_dropped += dropped_fdM; n_effective++; }
            }

            double SE = NAN, Z = NAN, p = NAN;
            if (n_effective >= 2 && isfinite(gw.fdM)) {
                double mean_dropped = sum_dropped / n_effective;
                double var_sum = 0;
                for (int b = 0; b < B; b++) if (isfinite(fdM_dropped[b])) {
                    double d = fdM_dropped[b] - mean_dropped;
                    var_sum += d * d;
                }
                // Standard delete-1 jackknife SE: SE = sqrt((n-1)/n · Σ (θᵢ - mean(θ))²)
                SE = sqrt((double)(n_effective - 1) / n_effective * var_sum);
                if (SE > 0) {
                    Z = gw.fdM / SE;
                    p = erfc(fabs(Z) * M_SQRT1_2);
                }
            }
            free(fdM_dropped);

            fprintf(fout, "GENOME\t%s\t%d\t%d\t%d\t%d",
                    filter_chr ? filter_chr : ".",
                    sites[0].pos, sites[n_sites-1].pos,
                    gw.n_sites, gw.n_informative);
            #define EM2(v) do { \
                if (isnan(v))      fprintf(fout, "\tNA"); \
                else if (isinf(v)) fprintf(fout, "\tInf"); \
                else               fprintf(fout, "\t%.6g", (double)(v)); \
            } while (0)
            EM2(gw.D); EM2(gw.fdM);
            EM2(gw.sum_ABBA); EM2(gw.sum_BABA);
            EM2(gw.sum_num); EM2(gw.sum_denom);
            EM2(SE); EM2(Z); EM2(p);
            #undef EM2
            fputc('\n', fout);
            fprintf(stderr, "[fdM] Genome jackknife (%d blocks): fdM=%.4g SE=%.4g Z=%.3f p=%.4g\n",
                    n_effective, gw.fdM, SE, Z, p);
        }
    }

    if (out_path) fclose(fout);
    free(sites); free(wins);
    fprintf(stderr, "[fdM] Done: %d windows\n", n_wins);
    return 0;
}
