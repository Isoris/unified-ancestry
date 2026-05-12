// =============================================================================
// codon_stats.c — 4-fold degenerate sites + Nei-Gojobori (NG86) and
//                  Yang-Nielsen (YN00) dN/dS / Ka/Ks for aligned CDS pairs.
//
// Pure C, no deps beyond libc + libm. Built for "click a gene → see evolution"
// speed: precomputed 64×64 substitution tables, single linear pass per pair,
// OpenMP across pairs.
//
// Methods:
//   NG86 — Nei & Gojobori 1986. Uniform path averaging, Jukes-Cantor distance.
//   YN00 — Yang & Nielsen 2000 (simplified). κ-weighted path averaging,
//          K80 (Kimura two-parameter) distance. κ estimated from the data
//          (K80 moment estimator on aligned codons), or user-supplied
//          via --kappa.
//
// Both methods always emit per-pair transition / transversion breakdowns
// (sd_ts, sd_tv, nd_ts, nd_tv, n_4d_ts, n_4d_tv) so downstream code can
// recompute distances any way it likes. dS / dN / d4d columns hold the
// method-appropriate corrected distance.
//
// Input (one of):
//   --pairs <pairs.tsv>   TSV: pair_id<TAB>seqA<TAB>seqB (sequences inline)
//   --fasta <align.fa>    Multi-FASTA, consecutive records form pairs
//                          (>gene_A then >gene_B; trailing A/B/1/2 stripped)
//
// Output: TSV to stdout (or --out <file>); JSON schema in
//   engines/schemas/codon_stats.output.schema.json (schema_version=codon_stats_v2)
//
// Compile: gcc -O3 -march=native -fopenmp -o codon_stats codon_stats.c -lm
//
// Citations:
//   Nei M, Gojobori T (1986) Mol Biol Evol 3:418-426
//   Yang Z, Nielsen R (2000) Mol Biol Evol 17:32-43
//   Kimura M (1980) J Mol Evol 16:111-120
//   4D-site neutral rate proxy: Bulmer M (1991) Genetics 129:897
// =============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define SCHEMA_VERSION "codon_stats_v2"
#define N_CODONS 64

// ── Genetic code (NCBI table 1, standard). codon = (nt1<<4)|(nt2<<2)|nt3,
//    A=0 C=1 G=2 T=3.
//                AAx       ACx       AGx       ATx       CAx       CCx       CGx       CTx
//                GAx       GCx       GGx       GTx       TAx       TCx       TGx       TTx
static const char STD_CODE[N_CODONS + 1] =
    "KNKN" "TTTT" "RSRS" "IIMI"
    "QHQH" "PPPP" "RRRR" "LLLL"
    "EDED" "AAAA" "GGGG" "VVVV"
    "*Y*Y" "SSSS" "*CWC" "LFLF";

static char  aa_of[N_CODONS];
static int   is_stop[N_CODONS];
static int   is_4d_third[N_CODONS];   // 1 if 3rd-position is 4-fold degenerate

typedef struct {
    double sd_ts, sd_tv;  // synonymous transitions / transversions
    double nd_ts, nd_tv;  // nonsynonymous transitions / transversions
} CodonSubs;

// Tables for each method (computed once at startup; YN00 tables rebuilt after κ is known).
static double    S_codon_ng[N_CODONS];        // NG86 synonymous-site count per codon
static CodonSubs sub_ng[N_CODONS][N_CODONS];  // NG86 path-averaged substitution counts
static double    S_codon_yn[N_CODONS];        // YN00 κ-weighted synonymous-site count
static CodonSubs sub_yn[N_CODONS][N_CODONS];  // YN00 κ-weighted substitution counts

typedef enum { METHOD_NG86 = 0, METHOD_YN00 = 1 } Method;

static inline int nt_enc(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1;
    }
}
static inline int is_transition(int a, int b) { return ((a ^ b) == 2); }

// ── Path enumeration shared between NG86 (uniform weight) and YN00 (κ weight) ──

static CodonSubs compute_pair_subs(int c1, int c2, double kappa) {
    CodonSubs r = {0, 0, 0, 0};
    if (c1 == c2 || is_stop[c1] || is_stop[c2]) return r;
    int n1[3] = { (c1 >> 4) & 3, (c1 >> 2) & 3, c1 & 3 };
    int n2[3] = { (c2 >> 4) & 3, (c2 >> 2) & 3, c2 & 3 };
    int diff[3], nd = 0;
    for (int p = 0; p < 3; p++) if (n1[p] != n2[p]) diff[nd++] = p;

    if (nd == 1) {
        int p = diff[0];
        int ts = is_transition(n1[p], n2[p]);
        int syn = (aa_of[c1] == aa_of[c2]);
        if (syn && ts)       r.sd_ts = 1;
        else if (syn)        r.sd_tv = 1;
        else if (ts)         r.nd_ts = 1;
        else                 r.nd_tv = 1;
        return r;
    }

    static const int PERM2[2][2] = {{0,1},{1,0}};
    static const int PERM3[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    int n_perms = (nd == 2) ? 2 : 6;

    double tot_w = 0;
    double sum_sts = 0, sum_stv = 0, sum_nts = 0, sum_ntv = 0;

    for (int pi = 0; pi < n_perms; pi++) {
        int cur[3] = { n1[0], n1[1], n1[2] };
        int cur_c = c1;
        double w = 1.0;
        double sts = 0, stv = 0, nts = 0, ntv = 0;
        int valid = 1;
        for (int j = 0; j < nd; j++) {
            int pos = (nd == 2) ? diff[PERM2[pi][j]] : diff[PERM3[pi][j]];
            int from = cur[pos], to = n2[pos];
            cur[pos] = to;
            int next_c = (cur[0] << 4) | (cur[1] << 2) | cur[2];
            if (is_stop[next_c]) { valid = 0; break; }
            int ts = is_transition(from, to);
            w *= (ts ? kappa : 1.0);
            int syn = (aa_of[cur_c] == aa_of[next_c]);
            if (syn && ts) sts++;
            else if (syn)  stv++;
            else if (ts)   nts++;
            else           ntv++;
            cur_c = next_c;
        }
        if (valid) {
            tot_w += w;
            sum_sts += w * sts; sum_stv += w * stv;
            sum_nts += w * nts; sum_ntv += w * ntv;
        }
    }
    if (tot_w > 0) {
        r.sd_ts = sum_sts / tot_w;
        r.sd_tv = sum_stv / tot_w;
        r.nd_ts = sum_nts / tot_w;
        r.nd_tv = sum_ntv / tot_w;
    }
    return r;
}

// ── S per codon: weighted average of (syn 1-step neighbours) over the 3 alt nts at each position.
//                NG86 sets weights=1 (denominator 3). YN00 weights transitions by κ, transversions by 1
//                (denominator κ+2). Stop neighbours are excluded from the numerator only — matches the
//                standard convention used by both NG86 and YN00.
static double compute_S_codon(int c, double kappa) {
    if (is_stop[c]) return 0;
    int nt[3] = { (c >> 4) & 3, (c >> 2) & 3, c & 3 };
    double S = 0;
    double denom = kappa + 2.0;  // sum of weights over the 3 alternative nts at one position
    for (int pos = 0; pos < 3; pos++) {
        double w_syn = 0;
        for (int alt = 0; alt < 4; alt++) {
            if (alt == nt[pos]) continue;
            int nn[3] = { nt[0], nt[1], nt[2] };
            nn[pos] = alt;
            int cc = (nn[0] << 4) | (nn[1] << 2) | nn[2];
            if (is_stop[cc]) continue;
            double w = is_transition(nt[pos], alt) ? kappa : 1.0;
            if (aa_of[cc] == aa_of[c]) w_syn += w;
        }
        S += w_syn / denom;
    }
    return S;
}

static void build_tables(double kappa, double* S_dst, CodonSubs (*sub_dst)[N_CODONS]) {
    for (int c = 0; c < N_CODONS; c++) S_dst[c] = compute_S_codon(c, kappa);
    for (int c1 = 0; c1 < N_CODONS; c1++)
        for (int c2 = 0; c2 < N_CODONS; c2++)
            sub_dst[c1][c2] = compute_pair_subs(c1, c2, kappa);
}

static void init_static_tables(void) {
    for (int c = 0; c < N_CODONS; c++) {
        aa_of[c] = STD_CODE[c];
        is_stop[c] = (STD_CODE[c] == '*');
    }
    for (int c = 0; c < N_CODONS; c++) {
        int base = c & 0x3C; // mask out position 3
        if (is_stop[base]) { is_4d_third[c] = 0; continue; }
        char aa = aa_of[base];
        int all_same = 1;
        for (int x = 0; x < 4; x++) {
            int cc = base | x;
            if (is_stop[cc] || aa_of[cc] != aa) { all_same = 0; break; }
        }
        is_4d_third[c] = all_same;
    }
    build_tables(1.0, S_codon_ng, sub_ng);  // NG86 = κ=1
    // YN00 tables built later, after κ is known.
}

// ── Distance corrections ────────────────────────────────────────────────────

static double jc_correct(double p) {
    if (isnan(p)) return NAN;
    if (p < 0) return NAN;
    if (p >= 0.75) return INFINITY;
    return -0.75 * log(1.0 - (4.0 / 3.0) * p);
}

static double k80_correct(double P, double Q) {
    if (isnan(P) || isnan(Q)) return NAN;
    double a = 1.0 - 2.0 * P - Q;
    double b = 1.0 - 2.0 * Q;
    if (a <= 0 || b <= 0) return INFINITY;
    return -0.5 * log(a) - 0.25 * log(b);
}

// ── κ estimation from aligned codons (K80 moment estimator) ─────────────────

typedef struct {
    char* id;
    char* a;
    char* b;
    int   len;
} Pair;
typedef struct { Pair* items; int n, cap; } PairList;

static double estimate_kappa(const PairList* pl) {
    long ts = 0, tv = 0, sites = 0;
    for (int k = 0; k < pl->n; k++) {
        const char* sa = pl->items[k].a;
        const char* sb = pl->items[k].b;
        int len = pl->items[k].len;
        for (int i = 0; i + 2 < len; i += 3) {
            int n1[3], n2[3], gap = 0;
            for (int p = 0; p < 3; p++) {
                n1[p] = nt_enc(sa[i+p]); n2[p] = nt_enc(sb[i+p]);
                if (n1[p] < 0 || n2[p] < 0) { gap = 1; break; }
            }
            if (gap) continue;
            int c1 = (n1[0]<<4)|(n1[1]<<2)|n1[2];
            int c2 = (n2[0]<<4)|(n2[1]<<2)|n2[2];
            if (is_stop[c1] || is_stop[c2]) continue;
            for (int p = 0; p < 3; p++) {
                sites++;
                if (n1[p] != n2[p]) {
                    if (is_transition(n1[p], n2[p])) ts++;
                    else                              tv++;
                }
            }
        }
    }
    if (sites < 10) return 2.0;
    double P = (double)ts / sites;
    double Q = (double)tv / sites;
    if (1.0 - 2.0*P - Q <= 0 || 1.0 - 2.0*Q <= 0) return 2.0;
    double at = -0.5 * log(1.0 - 2.0*P - Q) + 0.25 * log(1.0 - 2.0*Q);
    double bt = -0.25 * log(1.0 - 2.0*Q);
    if (bt <= 0) return 2.0;
    double k = at / bt;
    if (!isfinite(k) || k <= 0) return 2.0;
    return k;
}

// ── Per-pair compute ────────────────────────────────────────────────────────

typedef struct {
    int  n_aligned_bp;
    int  n_codons_total;
    int  n_codons_used;
    int  n_codons_gap;
    int  n_codons_stop;
    int  n_4d_sites;
    int  n_4d_diffs;
    int  n_4d_ts;
    int  n_4d_tv;
    double S, N;
    double sd_ts, sd_tv;
    double nd_ts, nd_tv;
    double sd, nd;
    double pS, pN;
    double dS, dN;
    double omega;
    double p4d, d4d;
} PairStats;

static void compute_pair(const char* sa, const char* sb, int len,
                         Method method,
                         const double* S_tab, const CodonSubs (*sub_tab)[N_CODONS],
                         PairStats* out) {
    memset(out, 0, sizeof(*out));
    out->n_aligned_bp = len;
    out->n_codons_total = len / 3;
    if (len < 3) return;
    for (int i = 0; i + 2 < len; i += 3) {
        int n1[3], n2[3], gap = 0;
        for (int p = 0; p < 3; p++) {
            n1[p] = nt_enc(sa[i+p]); n2[p] = nt_enc(sb[i+p]);
            if (n1[p] < 0 || n2[p] < 0) { gap = 1; break; }
        }
        if (gap) { out->n_codons_gap++; continue; }
        int c1 = (n1[0]<<4)|(n1[1]<<2)|n1[2];
        int c2 = (n2[0]<<4)|(n2[1]<<2)|n2[2];
        if (is_stop[c1] || is_stop[c2]) { out->n_codons_stop++; continue; }
        out->n_codons_used++;
        out->S    += (S_tab[c1] + S_tab[c2]) / 2.0;
        out->sd_ts += sub_tab[c1][c2].sd_ts;
        out->sd_tv += sub_tab[c1][c2].sd_tv;
        out->nd_ts += sub_tab[c1][c2].nd_ts;
        out->nd_tv += sub_tab[c1][c2].nd_tv;
        if (is_4d_third[c1] && is_4d_third[c2] && n1[0] == n2[0] && n1[1] == n2[1]) {
            out->n_4d_sites++;
            if (n1[2] != n2[2]) {
                out->n_4d_diffs++;
                if (is_transition(n1[2], n2[2])) out->n_4d_ts++;
                else                              out->n_4d_tv++;
            }
        }
    }
    out->N  = 3.0 * out->n_codons_used - out->S;
    out->sd = out->sd_ts + out->sd_tv;
    out->nd = out->nd_ts + out->nd_tv;
    out->pS = (out->S > 0) ? out->sd / out->S : NAN;
    out->pN = (out->N > 0) ? out->nd / out->N : NAN;
    if (method == METHOD_NG86) {
        out->dS = jc_correct(out->pS);
        out->dN = jc_correct(out->pN);
    } else {
        double Ps = (out->S > 0) ? out->sd_ts / out->S : NAN;
        double Qs = (out->S > 0) ? out->sd_tv / out->S : NAN;
        double Pn = (out->N > 0) ? out->nd_ts / out->N : NAN;
        double Qn = (out->N > 0) ? out->nd_tv / out->N : NAN;
        out->dS = k80_correct(Ps, Qs);
        out->dN = k80_correct(Pn, Qn);
    }
    out->omega = (isfinite(out->dS) && out->dS > 0) ? out->dN / out->dS : NAN;
    out->p4d = (out->n_4d_sites > 0) ? (double)out->n_4d_diffs / out->n_4d_sites : NAN;
    if (method == METHOD_NG86 || out->n_4d_sites == 0) {
        out->d4d = jc_correct(out->p4d);
    } else {
        double P4 = (double)out->n_4d_ts / out->n_4d_sites;
        double Q4 = (double)out->n_4d_tv / out->n_4d_sites;
        out->d4d = k80_correct(P4, Q4);
    }
}

// ── I/O helpers ─────────────────────────────────────────────────────────────

static void plist_push(PairList* L, const char* id, char* a, char* b, int len) {
    if (L->n == L->cap) {
        L->cap = L->cap ? L->cap * 2 : 64;
        L->items = (Pair*)realloc(L->items, L->cap * sizeof(Pair));
    }
    L->items[L->n].id = strdup(id);
    L->items[L->n].a = a;
    L->items[L->n].b = b;
    L->items[L->n].len = len;
    L->n++;
}

static char* read_all(FILE* f, size_t* out_len) {
    size_t cap = 1 << 16, n = 0;
    char* buf = (char*)malloc(cap);
    for (;;) {
        if (n == cap) { cap *= 2; buf = (char*)realloc(buf, cap); }
        size_t r = fread(buf + n, 1, cap - n, f);
        if (r == 0) break;
        n += r;
    }
    *out_len = n;
    return buf;
}

static int read_pairs_tsv(const char* path, PairList* out) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[codon_stats] Cannot open %s\n", path); return 0; }
    char* line = NULL;
    size_t cap = 0;
    ssize_t got;
    int header_skipped = 0;
    while ((got = getline(&line, &cap, f)) != -1) {
        if (got > 0 && line[got - 1] == '\n') line[--got] = 0;
        if (got > 0 && line[got - 1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        if (!header_skipped && strncasecmp(line, "pair_id", 7) == 0) {
            header_skipped = 1; continue;
        }
        header_skipped = 1;
        char* t1 = strchr(line, '\t'); if (!t1) continue;
        char* t2 = strchr(t1 + 1, '\t'); if (!t2) continue;
        *t1 = 0; *t2 = 0;
        const char* id = line;
        char* a = strdup(t1 + 1);
        char* b = strdup(t2 + 1);
        int la = (int)strlen(a), lb = (int)strlen(b);
        if (la != lb) {
            fprintf(stderr, "[codon_stats] WARN: pair '%s' lengths differ (%d vs %d) — skipped\n",
                    id, la, lb);
            free(a); free(b); continue;
        }
        plist_push(out, id, a, b, la);
    }
    free(line);
    fclose(f);
    return out->n;
}

static int read_pairs_fasta(const char* path, PairList* out) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[codon_stats] Cannot open %s\n", path); return 0; }
    size_t flen;
    char* buf = read_all(f, &flen);
    fclose(f);

    typedef struct { char* name; char* seq; int slen; } Rec;
    Rec* recs = NULL;
    int nr = 0, cap = 0;

    size_t i = 0;
    while (i < flen) {
        while (i < flen && (buf[i] == '\n' || buf[i] == '\r' || buf[i] == ' ' || buf[i] == '\t')) i++;
        if (i >= flen) break;
        if (buf[i] != '>') { i++; continue; }
        i++;
        size_t hs = i;
        while (i < flen && buf[i] != '\n' && buf[i] != '\r') i++;
        size_t he = i;
        size_t name_end = hs;
        while (name_end < he && buf[name_end] != ' ' && buf[name_end] != '\t') name_end++;
        char* name = (char*)malloc(name_end - hs + 1);
        memcpy(name, buf + hs, name_end - hs); name[name_end - hs] = 0;
        while (i < flen && (buf[i] == '\n' || buf[i] == '\r')) i++;
        size_t scap = 256, sn = 0;
        char* seq = (char*)malloc(scap);
        while (i < flen && buf[i] != '>') {
            char c = buf[i++];
            if (c == '\n' || c == '\r' || c == ' ' || c == '\t') continue;
            if (sn == scap) { scap *= 2; seq = (char*)realloc(seq, scap); }
            seq[sn++] = c;
        }
        if (sn + 1 > scap) seq = (char*)realloc(seq, sn + 1);
        seq[sn] = 0;
        if (nr == cap) { cap = cap ? cap * 2 : 64; recs = (Rec*)realloc(recs, cap * sizeof(Rec)); }
        recs[nr].name = name;
        recs[nr].seq = seq;
        recs[nr].slen = (int)sn;
        nr++;
    }
    free(buf);

    for (int k = 0; k + 1 < nr; k += 2) {
        char id[1024];
        strncpy(id, recs[k].name, sizeof(id) - 1);
        id[sizeof(id) - 1] = 0;
        int L = (int)strlen(id);
        if (L > 2 && id[L - 2] == '_' &&
            (id[L - 1] == 'A' || id[L - 1] == 'a' || id[L - 1] == 'B' || id[L - 1] == 'b' ||
             id[L - 1] == '1' || id[L - 1] == '2')) {
            id[L - 2] = 0;
        }
        if (recs[k].slen != recs[k + 1].slen) {
            fprintf(stderr, "[codon_stats] WARN: pair %s lengths differ (%d vs %d) — skipped\n",
                    id, recs[k].slen, recs[k + 1].slen);
            free(recs[k].name); free(recs[k].seq);
            free(recs[k + 1].name); free(recs[k + 1].seq);
            continue;
        }
        plist_push(out, id, recs[k].seq, recs[k + 1].seq, recs[k].slen);
        free(recs[k].name); free(recs[k + 1].name);
    }
    if (nr & 1) {
        fprintf(stderr, "[codon_stats] WARN: FASTA had odd number of records (%d); last one ignored\n", nr);
        free(recs[nr - 1].name); free(recs[nr - 1].seq);
    }
    free(recs);
    return out->n;
}

static void print_header(FILE* fout, Method method, double kappa) {
    fprintf(fout, "# schema_version=%s method=%s",
            SCHEMA_VERSION, method == METHOD_NG86 ? "NG86" : "YN00");
    if (method == METHOD_YN00) fprintf(fout, " kappa=%.6f", kappa);
    fputc('\n', fout);
    fprintf(fout,
        "pair_id\tn_aligned_bp\tn_codons_total\tn_codons_used\tn_codons_gap\tn_codons_stop\t"
        "n_4d_sites\tn_4d_diffs\tn_4d_ts\tn_4d_tv\tS\tN\tsd\tnd\tsd_ts\tsd_tv\tnd_ts\tnd_tv\t"
        "pS\tpN\tdS\tdN\tomega\tp4d\td4d\n");
}

static void emit_field_double(FILE* fout, double v) {
    if (isnan(v))      fprintf(fout, "\tNA");
    else if (isinf(v)) fprintf(fout, "\tInf");
    else               fprintf(fout, "\t%.6f", v);
}

static void emit_row(FILE* fout, const char* pair_id, const PairStats* s) {
    fprintf(fout, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
            pair_id, s->n_aligned_bp, s->n_codons_total, s->n_codons_used,
            s->n_codons_gap, s->n_codons_stop,
            s->n_4d_sites, s->n_4d_diffs, s->n_4d_ts, s->n_4d_tv,
            s->S, s->N, s->sd, s->nd,
            s->sd_ts, s->sd_tv, s->nd_ts, s->nd_tv);
    emit_field_double(fout, s->pS);
    emit_field_double(fout, s->pN);
    emit_field_double(fout, s->dS);
    emit_field_double(fout, s->dN);
    emit_field_double(fout, s->omega);
    emit_field_double(fout, s->p4d);
    emit_field_double(fout, s->d4d);
    fputc('\n', fout);
}

// ── Main ────────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "codon_stats — 4D sites + NG86/YN00 dN/dS/Ka/Ks (fast, pure C)\n"
        "\n"
        "Input (one of):\n"
        "  --pairs <pairs.tsv>      TSV: pair_id<TAB>seqA<TAB>seqB (sequences inline)\n"
        "  --fasta <align.fa>       Multi-FASTA, consecutive records form pairs\n"
        "\n"
        "Options:\n"
        "  --method ng86|yn00       NG86 = Nei-Gojobori (default; JC distance, κ=1)\n"
        "                            YN00 = Yang-Nielsen (κ-weighted paths, K80 distance)\n"
        "  --kappa <value>          YN00 only: user-supplied ts/tv ratio. If omitted,\n"
        "                            estimated from the data (K80 moments on aligned codons).\n"
        "  --out <file>             Output TSV (default: stdout)\n"
        "  --ncores N               OpenMP threads across pairs (default: 1)\n"
        "  --no_header              Skip header line\n"
        "  -h, --help               This help\n"
        "\n"
        "Output columns (TSV; schema in engines/schemas/codon_stats.output.schema.json):\n"
        "  pair_id n_aligned_bp n_codons_total n_codons_used n_codons_gap n_codons_stop\n"
        "  n_4d_sites n_4d_diffs n_4d_ts n_4d_tv S N sd nd sd_ts sd_tv nd_ts nd_tv\n"
        "  pS pN dS dN omega p4d d4d\n");
}

int main(int argc, char** argv) {
    const char *pairs_path = NULL, *fasta_path = NULL, *out_path = NULL;
    const char *method_str = "ng86";
    double kappa_user = -1.0;
    int ncores = 1, no_header = 0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--pairs") && i + 1 < argc)   pairs_path = argv[++i];
        else if (!strcmp(argv[i], "--fasta") && i + 1 < argc)   fasta_path = argv[++i];
        else if (!strcmp(argv[i], "--method") && i + 1 < argc)  method_str = argv[++i];
        else if (!strcmp(argv[i], "--kappa") && i + 1 < argc)   kappa_user = atof(argv[++i]);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc)     out_path = argv[++i];
        else if (!strcmp(argv[i], "--ncores") && i + 1 < argc)  ncores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--no_header"))               no_header = 1;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
    }
    if (!pairs_path && !fasta_path) { print_usage(); return 1; }

    Method method;
    if (!strcasecmp(method_str, "ng86"))      method = METHOD_NG86;
    else if (!strcasecmp(method_str, "yn00")) method = METHOD_YN00;
    else { fprintf(stderr, "[codon_stats] Unknown method: %s (use ng86 or yn00)\n", method_str); return 2; }

    #ifdef _OPENMP
    omp_set_num_threads(ncores > 0 ? ncores : 1);
    #else
    (void)ncores;
    #endif

    init_static_tables();

    PairList pairs = {0};
    if (pairs_path) read_pairs_tsv(pairs_path, &pairs);
    if (fasta_path) read_pairs_fasta(fasta_path, &pairs);
    fprintf(stderr, "[codon_stats] %d pair(s) loaded\n", pairs.n);
    if (pairs.n == 0) return 0;

    double kappa_used = 1.0;
    const double* S_tab; const CodonSubs (*sub_tab)[N_CODONS];
    if (method == METHOD_NG86) {
        S_tab = S_codon_ng;
        sub_tab = sub_ng;
    } else {
        kappa_used = (kappa_user > 0) ? kappa_user : estimate_kappa(&pairs);
        fprintf(stderr, "[codon_stats] YN00: κ = %.4f (%s)\n", kappa_used,
                kappa_user > 0 ? "user-supplied" : "estimated from data");
        build_tables(kappa_used, S_codon_yn, sub_yn);
        S_tab = S_codon_yn;
        sub_tab = sub_yn;
    }

    PairStats* results = (PairStats*)calloc(pairs.n, sizeof(PairStats));

    #pragma omp parallel for schedule(static)
    for (int k = 0; k < pairs.n; k++) {
        compute_pair(pairs.items[k].a, pairs.items[k].b, pairs.items[k].len,
                     method, S_tab, sub_tab, &results[k]);
    }

    FILE* fout = out_path ? fopen(out_path, "w") : stdout;
    if (!fout) { fprintf(stderr, "[codon_stats] Cannot open %s for writing\n", out_path); return 3; }
    if (!no_header) print_header(fout, method, kappa_used);
    for (int k = 0; k < pairs.n; k++) emit_row(fout, pairs.items[k].id, &results[k]);
    if (out_path) fclose(fout);

    for (int k = 0; k < pairs.n; k++) {
        free(pairs.items[k].id);
        free(pairs.items[k].a);
        free(pairs.items[k].b);
    }
    free(pairs.items);
    free(results);
    fprintf(stderr, "[codon_stats] Done: %d pair(s)\n", pairs.n);
    return 0;
}
