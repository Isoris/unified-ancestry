// =============================================================================
// pi_NS.c — π_S, π_N, π_N/π_S per CDS locus, optionally per group.
//
// Population-level Nei-Gojobori-style nucleotide diversity broken down by
// site class. Conceptual analog of dNdSpiNpiS (PopPhyl) but in pure C.
//
// For each locus (= one multi-FASTA file, sequences are haplotypes/CDS samples):
//   For each group (default: ALL; optional --groups for HOM_STD / HET / HOM_INV
//                   or any partition):
//     For each codon position:
//       - Count valid haplotypes at this codon (no gaps, not a stop).
//       - mean_S per codon = average of NG86 S over valid haplotypes.
//       - Pairwise sd_ij, nd_ij over all valid pairs at this codon.
//       - Normalize pairwise sums by n_pairs_at_this_codon (handles missing
//         data per codon — same logic as Nei-Gojobori 1986 applied per pair).
//     π_S = Σ_codons (mean_sd_per_pair) / Σ_codons mean_S
//     π_N = Σ_codons (mean_nd_per_pair) / Σ_codons mean_N
//     π_N/π_S = ratio (NA if π_S == 0)
//
// Input:
//   --fasta <multi.fa>             One locus. Sequences are haplotypes.
//   --fasta_list <list.tsv>        Many loci. TSV: locus_id<TAB>fasta_path
//   --locus_id <id>                Locus id when using --fasta (default: filename).
//   --groups <file>                Optional. TSV: seq_name<TAB>group_label.
//                                   Sequences not in the file are skipped from
//                                   group rows but still counted in ALL.
//   --no_all                       Skip the "ALL" row (only emit per-group).
//   --out <f>                      Default stdout.
//   --ncores N                     OpenMP across loci.
//
// Output (TSV, one row per locus × group):
//   schema_version=pi_NS_v1
//   locus_id group n_seqs n_codons_total n_codons_used
//   S N total_sd total_nd pi_S pi_N pi_NS_ratio
//
// Citation: Nei M, Gojobori T (1986) Mol Biol Evol 3:418-426 (per-pair sd/nd).
//           PopPhyl dNdSpiNpiS (Galtier lab) for the population-level
//           formulation.
//
// Compile: gcc -O3 -march=native -fopenmp -o pi_NS pi_NS.c -lm
// =============================================================================

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define SCHEMA_VERSION "pi_NS_v2"
#define N_CODONS 64

// ── Genetic code (NCBI table 1). codon = (nt1<<4)|(nt2<<2)|nt3, A=0 C=1 G=2 T=3
static const char STD_CODE[N_CODONS + 1] =
    "KNKN" "TTTT" "RSRS" "IIMI"
    "QHQH" "PPPP" "RRRR" "LLLL"
    "EDED" "AAAA" "GGGG" "VVVV"
    "*Y*Y" "SSSS" "*CWC" "LFLF";

static char  aa_of[N_CODONS];
static int   is_stop[N_CODONS];
static double S_codon[N_CODONS];

// Per-codon-position degeneracy = number of synonymous 1-step alt nts.
//   0 → 0-fold (any nt change is nonsyn) — π0_fold proxy for πN
//   3 → 4-fold (all nt changes are syn)   — π4_fold proxy for πS (most-neutral)
//   1 or 2 → intermediate, ignored by π0/π4
// degeneracy[c][p] = -1 if codon c is a stop.
static int   degeneracy[N_CODONS][3];

typedef struct { double sd; double nd; } CodonSubs;
static CodonSubs sub_table[N_CODONS][N_CODONS];

static inline int nt_enc(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1;
    }
}

// NG86 path-averaged sd/nd (κ=1 uniform paths, stop intermediates dropped).
static CodonSubs compute_pair_subs(int c1, int c2) {
    CodonSubs r = {0, 0};
    if (c1 == c2 || is_stop[c1] || is_stop[c2]) return r;
    int n1[3] = { (c1 >> 4) & 3, (c1 >> 2) & 3, c1 & 3 };
    int n2[3] = { (c2 >> 4) & 3, (c2 >> 2) & 3, c2 & 3 };
    int diff[3], nd = 0;
    for (int p = 0; p < 3; p++) if (n1[p] != n2[p]) diff[nd++] = p;

    if (nd == 1) {
        if (aa_of[c1] == aa_of[c2]) r.sd = 1;
        else                        r.nd = 1;
        return r;
    }

    static const int PERM2[2][2] = {{0,1},{1,0}};
    static const int PERM3[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    int n_perms = (nd == 2) ? 2 : 6;
    double sum_s = 0, sum_n = 0;
    int n_paths = 0;
    for (int pi = 0; pi < n_perms; pi++) {
        int cur[3] = { n1[0], n1[1], n1[2] };
        int cur_c = c1;
        double ps = 0, pn = 0;
        int valid = 1;
        for (int j = 0; j < nd; j++) {
            int pos = (nd == 2) ? diff[PERM2[pi][j]] : diff[PERM3[pi][j]];
            cur[pos] = n2[pos];
            int next_c = (cur[0] << 4) | (cur[1] << 2) | cur[2];
            if (is_stop[next_c]) { valid = 0; break; }
            if (aa_of[cur_c] == aa_of[next_c]) ps++;
            else                                pn++;
            cur_c = next_c;
        }
        if (valid) { sum_s += ps; sum_n += pn; n_paths++; }
    }
    if (n_paths > 0) { r.sd = sum_s / n_paths; r.nd = sum_n / n_paths; }
    return r;
}

static double compute_S_codon(int c) {
    if (is_stop[c]) return 0;
    int nt[3] = { (c >> 4) & 3, (c >> 2) & 3, c & 3 };
    double S = 0;
    for (int pos = 0; pos < 3; pos++) {
        int n_syn = 0;
        for (int alt = 0; alt < 4; alt++) {
            if (alt == nt[pos]) continue;
            int nn[3] = { nt[0], nt[1], nt[2] };
            nn[pos] = alt;
            int cc = (nn[0] << 4) | (nn[1] << 2) | nn[2];
            if (is_stop[cc]) continue;
            if (aa_of[cc] == aa_of[c]) n_syn++;
        }
        S += (double)n_syn / 3.0;
    }
    return S;
}

static void init_codon_tables(void) {
    for (int c = 0; c < N_CODONS; c++) {
        aa_of[c]  = STD_CODE[c];
        is_stop[c] = (STD_CODE[c] == '*');
    }
    for (int c = 0; c < N_CODONS; c++) S_codon[c] = compute_S_codon(c);
    for (int c1 = 0; c1 < N_CODONS; c1++)
        for (int c2 = 0; c2 < N_CODONS; c2++)
            sub_table[c1][c2] = compute_pair_subs(c1, c2);
    // Per-position degeneracy: count syn 1-step alt nts at each codon position.
    for (int c = 0; c < N_CODONS; c++) {
        if (is_stop[c]) {
            for (int p = 0; p < 3; p++) degeneracy[c][p] = -1;
            continue;
        }
        int nt[3] = { (c >> 4) & 3, (c >> 2) & 3, c & 3 };
        for (int p = 0; p < 3; p++) {
            int n_syn = 0;
            for (int alt = 0; alt < 4; alt++) {
                if (alt == nt[p]) continue;
                int nn[3] = { nt[0], nt[1], nt[2] };
                nn[p] = alt;
                int cc = (nn[0] << 4) | (nn[1] << 2) | nn[2];
                if (is_stop[cc]) continue;
                if (aa_of[cc] == aa_of[c]) n_syn++;
            }
            degeneracy[c][p] = n_syn;
        }
    }
}

// ── Haplotype storage ─────────────────────────────────────────────────────

typedef struct {
    char* name;
    char* seq;
    int   len;
} Hap;

typedef struct {
    Hap* items;
    int  n;
    int  cap;
} HapList;

static void hl_push(HapList* L, char* name, char* seq, int len) {
    if (L->n == L->cap) {
        L->cap = L->cap ? L->cap * 2 : 64;
        L->items = (Hap*)realloc(L->items, L->cap * sizeof(Hap));
    }
    L->items[L->n].name = name;
    L->items[L->n].seq = seq;
    L->items[L->n].len = len;
    L->n++;
}

static void hl_free(HapList* L) {
    for (int i = 0; i < L->n; i++) {
        free(L->items[i].name);
        free(L->items[i].seq);
    }
    free(L->items);
    L->items = NULL; L->n = L->cap = 0;
}

// Read a multi-FASTA from `path` into haplotype list. Returns count.
static int read_multi_fasta(const char* path, HapList* out) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[pi_NS] Cannot open %s\n", path); return 0; }
    fseek(f, 0, SEEK_END);
    long flen = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buf = (char*)malloc((size_t)flen + 1);
    if (fread(buf, 1, (size_t)flen, f) != (size_t)flen) {
        fprintf(stderr, "[pi_NS] short read on %s\n", path);
        free(buf); fclose(f); return 0;
    }
    buf[flen] = 0;
    fclose(f);

    size_t i = 0, n_added = 0;
    while ((long)i < flen) {
        while ((long)i < flen && (buf[i] == '\n' || buf[i] == '\r' || buf[i] == ' ' || buf[i] == '\t')) i++;
        if ((long)i >= flen) break;
        if (buf[i] != '>') { i++; continue; }
        i++;
        size_t hs = i;
        while ((long)i < flen && buf[i] != '\n' && buf[i] != '\r') i++;
        size_t he = i;
        size_t name_end = hs;
        while (name_end < he && buf[name_end] != ' ' && buf[name_end] != '\t') name_end++;
        char* name = (char*)malloc(name_end - hs + 1);
        memcpy(name, buf + hs, name_end - hs); name[name_end - hs] = 0;
        while ((long)i < flen && (buf[i] == '\n' || buf[i] == '\r')) i++;

        size_t scap = 256, sn = 0;
        char* seq = (char*)malloc(scap);
        while ((long)i < flen && buf[i] != '>') {
            char c = buf[i++];
            if (c == '\n' || c == '\r' || c == ' ' || c == '\t') continue;
            if (sn == scap) { scap *= 2; seq = (char*)realloc(seq, scap); }
            seq[sn++] = c;
        }
        if (sn + 1 > scap) seq = (char*)realloc(seq, sn + 1);
        seq[sn] = 0;

        hl_push(out, name, seq, (int)sn);
        n_added++;
    }
    free(buf);
    return (int)n_added;
}

// ── Group mapping ─────────────────────────────────────────────────────────

typedef struct {
    char name[64];
    int* members;   // indices into HapList
    int  n;
    int  cap;
} HGroup;

typedef struct {
    HGroup* items;
    int     n;
    int     cap;
} GroupList;

static int group_find_or_add(GroupList* gl, const char* name) {
    for (int i = 0; i < gl->n; i++) if (!strcmp(gl->items[i].name, name)) return i;
    if (gl->n == gl->cap) {
        gl->cap = gl->cap ? gl->cap * 2 : 8;
        gl->items = (HGroup*)realloc(gl->items, gl->cap * sizeof(HGroup));
    }
    HGroup* g = &gl->items[gl->n];
    strncpy(g->name, name, sizeof(g->name) - 1); g->name[sizeof(g->name) - 1] = 0;
    g->members = NULL; g->n = 0; g->cap = 0;
    return gl->n++;
}

static void group_push_member(HGroup* g, int idx) {
    if (g->n == g->cap) {
        g->cap = g->cap ? g->cap * 2 : 16;
        g->members = (int*)realloc(g->members, g->cap * sizeof(int));
    }
    g->members[g->n++] = idx;
}

static void gl_free(GroupList* gl) {
    for (int i = 0; i < gl->n; i++) free(gl->items[i].members);
    free(gl->items);
    gl->items = NULL; gl->n = gl->cap = 0;
}

// Loads seq_name → group label. Returns map count.
typedef struct { char name[256]; char label[64]; } NameLabel;
static int load_group_file(const char* path, NameLabel** out, int* n_out) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[pi_NS] Cannot open --groups %s\n", path); return 0; }
    int cap = 64, n = 0;
    NameLabel* arr = (NameLabel*)malloc(cap * sizeof(NameLabel));
    char* line = NULL; size_t lcap = 0; ssize_t got;
    int header_seen = 0;
    while ((got = getline(&line, &lcap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        char nm[256] = "", lab[64] = "";
        if (sscanf(line, "%255s %63s", nm, lab) < 2) continue;
        if (!header_seen && (!strcasecmp(nm, "seq_name") || !strcasecmp(nm, "name") ||
                              !strcasecmp(nm, "sample") || !strcasecmp(nm, "haplotype"))) {
            header_seen = 1; continue;
        }
        header_seen = 1;
        if (n == cap) { cap *= 2; arr = (NameLabel*)realloc(arr, cap * sizeof(NameLabel)); }
        strncpy(arr[n].name, nm, sizeof(arr[n].name) - 1);  arr[n].name[sizeof(arr[n].name) - 1] = 0;
        strncpy(arr[n].label, lab, sizeof(arr[n].label) - 1); arr[n].label[sizeof(arr[n].label) - 1] = 0;
        n++;
    }
    free(line);
    fclose(f);
    *out = arr; *n_out = n;
    return n;
}

// ── π_S / π_N per group ────────────────────────────────────────────────────

typedef struct {
    int n_seqs;
    int n_codons_total;
    int n_codons_used;
    double S_total;
    double N_total;
    double total_sd;
    double total_nd;
    double pi_S;
    double pi_N;
    double pi_NS_ratio;
    // π0-fold / π4-fold (per-codon-position pairwise nt diffs at sites where the
    // reference codon's degeneracy is 0 or 3). π4 ≈ neutral, π0 ≈ all nonsyn.
    int    n_4fold_sites;
    int    n_0fold_sites;
    double total_diffs_4fold;          // sum over 4-fold sites of mean pairwise nt diff
    double total_diffs_0fold;
    double pi_4fold;
    double pi_0fold;
    double pi_04_ratio;                // π0 / π4, ≈ πN/πS proxy
    // Coverage correction.
    double coverage_factor;
    double pi_S_corrected;
    double pi_N_corrected;
    double pi_4fold_corrected;
    double pi_0fold_corrected;
    // Bootstrap CIs (set by bootstrap_pi; NaN if --bootstrap 0).
    int    bootstrap_reps;
    double pi_S_lo, pi_S_hi;
    double pi_N_lo, pi_N_hi;
    double pi_NS_lo, pi_NS_hi;
    double pi_4fold_lo, pi_4fold_hi;
    double pi_0fold_lo, pi_0fold_hi;
    double pi_04_lo, pi_04_hi;
    // Per-codon contributions, kept for bootstrap (free with free_codon_contribs).
    double* codon_S;
    double* codon_sd;
    double* codon_nd;
    double* codon_n4;                  // # 4-fold positions contributed by this codon (0..3)
    double* codon_d4;                  // sum of mean pairwise nt diffs at those positions
    double* codon_n0;
    double* codon_d0;
    int     n_codon_contribs;
} PiStats;

static void free_codon_contribs(PiStats* s) {
    free(s->codon_S);  s->codon_S = NULL;
    free(s->codon_sd); s->codon_sd = NULL;
    free(s->codon_nd); s->codon_nd = NULL;
    free(s->codon_n4); s->codon_n4 = NULL;
    free(s->codon_d4); s->codon_d4 = NULL;
    free(s->codon_n0); s->codon_n0 = NULL;
    free(s->codon_d0); s->codon_d0 = NULL;
    s->n_codon_contribs = 0;
}

static void compute_pi_for_group(const Hap* haps, const int* members, int n_members,
                                 int seq_len, int keep_codon_contribs, PiStats* out) {
    memset(out, 0, sizeof(*out));
    out->n_seqs = n_members;
    out->pi_S = NAN; out->pi_N = NAN; out->pi_NS_ratio = NAN;
    out->pi_4fold = NAN; out->pi_0fold = NAN; out->pi_04_ratio = NAN;
    out->coverage_factor = 1.0;
    out->pi_S_corrected = NAN; out->pi_N_corrected = NAN;
    out->pi_4fold_corrected = NAN; out->pi_0fold_corrected = NAN;
    out->pi_S_lo = out->pi_S_hi = NAN;
    out->pi_N_lo = out->pi_N_hi = NAN;
    out->pi_NS_lo = out->pi_NS_hi = NAN;
    out->pi_4fold_lo = out->pi_4fold_hi = NAN;
    out->pi_0fold_lo = out->pi_0fold_hi = NAN;
    out->pi_04_lo = out->pi_04_hi = NAN;
    if (n_members < 2 || seq_len < 3) return;
    out->n_codons_total = seq_len / 3;

    int max_codons = out->n_codons_total;
    if (keep_codon_contribs && max_codons > 0) {
        out->codon_S  = (double*)malloc((size_t)max_codons * sizeof(double));
        out->codon_sd = (double*)malloc((size_t)max_codons * sizeof(double));
        out->codon_nd = (double*)malloc((size_t)max_codons * sizeof(double));
        out->codon_n4 = (double*)malloc((size_t)max_codons * sizeof(double));
        out->codon_d4 = (double*)malloc((size_t)max_codons * sizeof(double));
        out->codon_n0 = (double*)malloc((size_t)max_codons * sizeof(double));
        out->codon_d0 = (double*)malloc((size_t)max_codons * sizeof(double));
    }

    int* codons = (int*)malloc((size_t)n_members * sizeof(int));

    for (int pos = 0; pos + 2 < seq_len; pos += 3) {
        int n_valid = 0;
        int ref_codon = -1;  // first valid codon — used for degeneracy classification
        for (int k = 0; k < n_members; k++) {
            int idx = members[k];
            if (pos + 2 >= haps[idx].len) { codons[k] = -1; continue; }
            int n1 = nt_enc(haps[idx].seq[pos]);
            int n2 = nt_enc(haps[idx].seq[pos+1]);
            int n3 = nt_enc(haps[idx].seq[pos+2]);
            if (n1 < 0 || n2 < 0 || n3 < 0) { codons[k] = -1; continue; }
            int cc = (n1 << 4) | (n2 << 2) | n3;
            if (is_stop[cc]) { codons[k] = -1; continue; }
            codons[k] = cc;
            if (ref_codon < 0) ref_codon = cc;
            n_valid++;
        }
        if (n_valid < 2) continue;

        // πS / πN (NG86 path-averaged, codon-pair level).
        double mean_S = 0;
        for (int k = 0; k < n_members; k++) if (codons[k] >= 0) mean_S += S_codon[codons[k]];
        mean_S /= n_valid;

        double sd_pos = 0, nd_pos = 0;
        for (int i = 0; i < n_members; i++) {
            if (codons[i] < 0) continue;
            for (int j = i + 1; j < n_members; j++) {
                if (codons[j] < 0) continue;
                sd_pos += sub_table[codons[i]][codons[j]].sd;
                nd_pos += sub_table[codons[i]][codons[j]].nd;
            }
        }
        double n_pairs_here = 0.5 * n_valid * (n_valid - 1);
        double mean_sd = (n_pairs_here > 0) ? sd_pos / n_pairs_here : 0;
        double mean_nd = (n_pairs_here > 0) ? nd_pos / n_pairs_here : 0;

        out->S_total += mean_S;
        out->N_total += 3.0 - mean_S;
        out->total_sd += mean_sd;
        out->total_nd += mean_nd;

        // π0-fold / π4-fold (per-codon-position pairwise nt-diff at sites where
        // the reference codon's position is 0-fold or 4-fold degenerate).
        double n4_contrib = 0, d4_contrib = 0;
        double n0_contrib = 0, d0_contrib = 0;
        for (int p = 0; p < 3; p++) {
            int deg = degeneracy[ref_codon][p];
            if (deg != 0 && deg != 3) continue;
            int nt_pos = pos + p;
            int pair_diffs = 0, pair_count = 0;
            for (int i = 0; i < n_members; i++) {
                if (codons[i] < 0) continue;
                int ni = nt_enc(haps[members[i]].seq[nt_pos]);
                if (ni < 0) continue;
                for (int j = i + 1; j < n_members; j++) {
                    if (codons[j] < 0) continue;
                    int nj = nt_enc(haps[members[j]].seq[nt_pos]);
                    if (nj < 0) continue;
                    pair_count++;
                    if (ni != nj) pair_diffs++;
                }
            }
            if (pair_count == 0) continue;
            double mean_diff = (double)pair_diffs / pair_count;
            if (deg == 3) {
                out->n_4fold_sites++;
                out->total_diffs_4fold += mean_diff;
                n4_contrib += 1.0;
                d4_contrib += mean_diff;
            } else {
                out->n_0fold_sites++;
                out->total_diffs_0fold += mean_diff;
                n0_contrib += 1.0;
                d0_contrib += mean_diff;
            }
        }

        if (keep_codon_contribs && out->codon_S) {
            int k = out->n_codon_contribs;
            out->codon_S[k]  = mean_S;
            out->codon_sd[k] = mean_sd;
            out->codon_nd[k] = mean_nd;
            out->codon_n4[k] = n4_contrib;
            out->codon_d4[k] = d4_contrib;
            out->codon_n0[k] = n0_contrib;
            out->codon_d0[k] = d0_contrib;
            out->n_codon_contribs++;
        }
        out->n_codons_used++;
    }

    free(codons);

    if (out->n_codons_used == 0) return;
    out->pi_S = (out->S_total > 0) ? out->total_sd / out->S_total : NAN;
    out->pi_N = (out->N_total > 0) ? out->total_nd / out->N_total : NAN;
    if (isfinite(out->pi_S) && out->pi_S > 0) out->pi_NS_ratio = out->pi_N / out->pi_S;

    out->pi_4fold = (out->n_4fold_sites > 0) ? out->total_diffs_4fold / out->n_4fold_sites : NAN;
    out->pi_0fold = (out->n_0fold_sites > 0) ? out->total_diffs_0fold / out->n_0fold_sites : NAN;
    if (isfinite(out->pi_4fold) && out->pi_4fold > 0)
        out->pi_04_ratio = out->pi_0fold / out->pi_4fold;
}

// Apply coverage correction in place. Called after compute_pi_for_group.
static void apply_coverage_correction(PiStats* s, double coverage_factor) {
    s->coverage_factor = (coverage_factor > 0) ? coverage_factor : 1.0;
    double cf = s->coverage_factor;
    s->pi_S_corrected      = (isfinite(s->pi_S)      && cf > 0) ? s->pi_S      / cf : NAN;
    s->pi_N_corrected      = (isfinite(s->pi_N)      && cf > 0) ? s->pi_N      / cf : NAN;
    s->pi_4fold_corrected  = (isfinite(s->pi_4fold)  && cf > 0) ? s->pi_4fold  / cf : NAN;
    s->pi_0fold_corrected  = (isfinite(s->pi_0fold)  && cf > 0) ? s->pi_0fold  / cf : NAN;
}

static int dbl_cmp(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    if (x < y) return -1;
    if (x > y) return  1;
    return 0;
}

static double pct_after_dropna(double* arr, int n, double p) {
    int m = 0;
    for (int i = 0; i < n; i++) if (isfinite(arr[i])) arr[m++] = arr[i];
    if (m == 0) return NAN;
    qsort(arr, (size_t)m, sizeof(double), dbl_cmp);
    int idx = (int)floor(p * (m - 1) + 0.5);
    if (idx < 0) idx = 0;
    if (idx >= m) idx = m - 1;
    return arr[idx];
}

// Bootstrap 95% CIs by resampling codons (with replacement) from the stored
// per-codon contributions. Applies the coverage correction to each replicate.
static void bootstrap_pi(PiStats* s, int n_boot, double coverage_factor, unsigned int* rng_state) {
    s->bootstrap_reps = 0;
    if (n_boot <= 0 || s->n_codon_contribs < 2) return;
    s->bootstrap_reps = n_boot;
    int nc = s->n_codon_contribs;
    double cov = (coverage_factor > 0) ? coverage_factor : 1.0;

    double* repS  = (double*)malloc((size_t)n_boot * sizeof(double));
    double* repN  = (double*)malloc((size_t)n_boot * sizeof(double));
    double* repR  = (double*)malloc((size_t)n_boot * sizeof(double));
    double* rep4  = (double*)malloc((size_t)n_boot * sizeof(double));
    double* rep0  = (double*)malloc((size_t)n_boot * sizeof(double));
    double* rep04 = (double*)malloc((size_t)n_boot * sizeof(double));

    for (int b = 0; b < n_boot; b++) {
        double S_total = 0, N_total = 0, sd_total = 0, nd_total = 0;
        double n4_total = 0, d4_total = 0, n0_total = 0, d0_total = 0;
        for (int c = 0; c < nc; c++) {
            int k = (int)(rand_r(rng_state) % (unsigned int)nc);
            S_total += s->codon_S[k];
            N_total += 3.0 - s->codon_S[k];
            sd_total += s->codon_sd[k];
            nd_total += s->codon_nd[k];
            n4_total += s->codon_n4[k];
            d4_total += s->codon_d4[k];
            n0_total += s->codon_n0[k];
            d0_total += s->codon_d0[k];
        }
        double pS  = (S_total > 0)  ? (sd_total / S_total) / cov : NAN;
        double pN  = (N_total > 0)  ? (nd_total / N_total) / cov : NAN;
        double p4  = (n4_total > 0) ? (d4_total / n4_total) / cov : NAN;
        double p0  = (n0_total > 0) ? (d0_total / n0_total) / cov : NAN;
        repS[b] = pS;  repN[b] = pN;
        repR[b]  = (isfinite(pS) && pS > 0) ? pN / pS : NAN;
        rep4[b] = p4;  rep0[b] = p0;
        rep04[b] = (isfinite(p4) && p4 > 0) ? p0 / p4 : NAN;
    }

    s->pi_S_lo = pct_after_dropna(repS, n_boot, 0.025);
    s->pi_S_hi = pct_after_dropna(repS, n_boot, 0.975);
    s->pi_N_lo = pct_after_dropna(repN, n_boot, 0.025);
    s->pi_N_hi = pct_after_dropna(repN, n_boot, 0.975);
    s->pi_NS_lo = pct_after_dropna(repR, n_boot, 0.025);
    s->pi_NS_hi = pct_after_dropna(repR, n_boot, 0.975);
    s->pi_4fold_lo = pct_after_dropna(rep4, n_boot, 0.025);
    s->pi_4fold_hi = pct_after_dropna(rep4, n_boot, 0.975);
    s->pi_0fold_lo = pct_after_dropna(rep0, n_boot, 0.025);
    s->pi_0fold_hi = pct_after_dropna(rep0, n_boot, 0.975);
    s->pi_04_lo = pct_after_dropna(rep04, n_boot, 0.025);
    s->pi_04_hi = pct_after_dropna(rep04, n_boot, 0.975);

    free(repS); free(repN); free(repR);
    free(rep4); free(rep0); free(rep04);
}

// ── Output ────────────────────────────────────────────────────────────────

static void emit_num(FILE* f, double v) {
    if (isnan(v))      fprintf(f, "\tNA");
    else if (isinf(v)) fprintf(f, "\tInf");
    else               fprintf(f, "\t%.6g", v);
}

static void emit_row(FILE* fout, const char* locus_id, const char* group,
                     const PiStats* s, int with_bootstrap) {
    fprintf(fout, "%s\t%s\t%d\t%d\t%d",
            locus_id, group, s->n_seqs, s->n_codons_total, s->n_codons_used);
    emit_num(fout, s->S_total);
    emit_num(fout, s->N_total);
    emit_num(fout, s->total_sd);
    emit_num(fout, s->total_nd);
    emit_num(fout, s->pi_S);
    emit_num(fout, s->pi_N);
    emit_num(fout, s->pi_NS_ratio);
    fprintf(fout, "\t%d\t%d", s->n_4fold_sites, s->n_0fold_sites);
    emit_num(fout, s->total_diffs_4fold);
    emit_num(fout, s->total_diffs_0fold);
    emit_num(fout, s->pi_4fold);
    emit_num(fout, s->pi_0fold);
    emit_num(fout, s->pi_04_ratio);
    emit_num(fout, s->coverage_factor);
    emit_num(fout, s->pi_S_corrected);
    emit_num(fout, s->pi_N_corrected);
    emit_num(fout, s->pi_4fold_corrected);
    emit_num(fout, s->pi_0fold_corrected);
    if (with_bootstrap) {
        fprintf(fout, "\t%d", s->bootstrap_reps);
        emit_num(fout, s->pi_S_lo);
        emit_num(fout, s->pi_S_hi);
        emit_num(fout, s->pi_N_lo);
        emit_num(fout, s->pi_N_hi);
        emit_num(fout, s->pi_NS_lo);
        emit_num(fout, s->pi_NS_hi);
        emit_num(fout, s->pi_4fold_lo);
        emit_num(fout, s->pi_4fold_hi);
        emit_num(fout, s->pi_0fold_lo);
        emit_num(fout, s->pi_0fold_hi);
        emit_num(fout, s->pi_04_lo);
        emit_num(fout, s->pi_04_hi);
    }
    fputc('\n', fout);
}

// ── Processing one locus ──────────────────────────────────────────────────

typedef struct { char locus_id[256]; char fasta_path[512]; } LocusEntry;

// Per-locus coverage factor (paper: invariant-sites-derived expected fraction
// of bases called per CDS). Look up by locus_id; default 1.0 if not found.
typedef struct { char locus_id[256]; double factor; } LocusCoverage;

static double lookup_coverage(const LocusCoverage* cov, int n_cov, const char* locus_id) {
    for (int i = 0; i < n_cov; i++) if (!strcmp(cov[i].locus_id, locus_id)) return cov[i].factor;
    return 1.0;
}

static int load_coverage_tsv(const char* path, LocusCoverage** out, int* n_out) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[pi_NS] Cannot open --coverage_tsv %s\n", path); return 0; }
    int cap = 64, n = 0;
    LocusCoverage* arr = (LocusCoverage*)malloc((size_t)cap * sizeof(LocusCoverage));
    char* line = NULL; size_t lcap = 0; ssize_t got;
    int header_seen = 0;
    while ((got = getline(&line, &lcap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        char id[256]; double v;
        if (sscanf(line, "%255s %lf", id, &v) < 2) continue;
        if (!header_seen && !strcasecmp(id, "locus_id")) { header_seen = 1; continue; }
        header_seen = 1;
        if (n == cap) { cap *= 2; arr = (LocusCoverage*)realloc(arr, (size_t)cap * sizeof(LocusCoverage)); }
        strncpy(arr[n].locus_id, id, sizeof(arr[n].locus_id) - 1);
        arr[n].locus_id[sizeof(arr[n].locus_id) - 1] = 0;
        arr[n].factor = v;
        n++;
    }
    free(line);
    fclose(f);
    *out = arr; *n_out = n;
    return n;
}

static void process_locus(const char* locus_id, const char* fasta_path,
                          const NameLabel* group_map, int group_map_n,
                          int no_all,
                          double cov_factor, int n_bootstrap, unsigned int* rng_state,
                          FILE* fout) {
    HapList haps = {0};
    int n = read_multi_fasta(fasta_path, &haps);
    if (n == 0) { fprintf(stderr, "[pi_NS] %s: no sequences in %s\n", locus_id, fasta_path); return; }

    int min_len = haps.items[0].len, max_len = min_len;
    for (int i = 1; i < n; i++) {
        if (haps.items[i].len < min_len) min_len = haps.items[i].len;
        if (haps.items[i].len > max_len) max_len = haps.items[i].len;
    }
    if (min_len != max_len) {
        fprintf(stderr, "[pi_NS] %s: ragged alignment (lengths %d..%d). Using min=%d.\n",
                locus_id, min_len, max_len, min_len);
    }
    int seq_len = min_len - (min_len % 3);

    GroupList gl = {0};
    int all_idx = group_find_or_add(&gl, "ALL");
    for (int i = 0; i < n; i++) group_push_member(&gl.items[all_idx], i);

    if (group_map && group_map_n > 0) {
        for (int i = 0; i < n; i++) {
            const char* hn = haps.items[i].name;
            for (int j = 0; j < group_map_n; j++) {
                if (!strcmp(group_map[j].name, hn)) {
                    int gi = group_find_or_add(&gl, group_map[j].label);
                    group_push_member(&gl.items[gi], i);
                    break;
                }
            }
        }
    }

    int keep_codons = (n_bootstrap > 0) ? 1 : 0;
    for (int gi = 0; gi < gl.n; gi++) {
        if (no_all && gi == all_idx) continue;
        PiStats s;
        compute_pi_for_group(haps.items, gl.items[gi].members,
                              gl.items[gi].n, seq_len, keep_codons, &s);
        apply_coverage_correction(&s, cov_factor);
        if (n_bootstrap > 0) {
            bootstrap_pi(&s, n_bootstrap, cov_factor, rng_state);
        }
        emit_row(fout, locus_id, gl.items[gi].name, &s, n_bootstrap > 0);
        free_codon_contribs(&s);
    }

    gl_free(&gl);
    hl_free(&haps);
}

// ── Main ──────────────────────────────────────────────────────────────────

static void print_usage(void) {
    fprintf(stderr,
        "pi_NS — π_S, π_N, π_N/π_S per CDS locus, optionally per group\n"
        "\n"
        "Input (one of):\n"
        "  --fasta <multi.fa>       One locus. Sequences are haplotypes.\n"
        "  --fasta_list <list.tsv>  Many loci. TSV: locus_id<TAB>fasta_path\n"
        "  --locus_id <id>          Locus id (only with --fasta; default: basename)\n"
        "\n"
        "Grouping:\n"
        "  --groups <file>          Optional. TSV: seq_name<TAB>group_label\n"
        "  --no_all                 Skip the ALL row (only per-group rows).\n"
        "\n"
        "Coverage correction (paper: invariant-sites-derived fraction of bases called):\n"
        "  --coverage_factor F      Global per-locus factor (0..1). Default 1.0.\n"
        "  --coverage_tsv <file>    Per-locus TSV: locus_id<TAB>factor. Looked up\n"
        "                            by locus_id; falls back to --coverage_factor.\n"
        "                            π_X_corrected = π_X / factor.\n"
        "\n"
        "Bootstrap 95%% CI (codon resampling with replacement):\n"
        "  --bootstrap N            N bootstrap replicates (default 0 = off, paper used 1000).\n"
        "                            Emits pi_S_lo / pi_S_hi / pi_N_lo / pi_N_hi /\n"
        "                            pi_NS_lo / pi_NS_hi.\n"
        "  --seed K                 RNG seed for reproducibility (default time(NULL)).\n"
        "\n"
        "Other:\n"
        "  --out <file>             Output TSV (default: stdout)\n"
        "  --ncores N               OpenMP across loci (default: 1)\n"
        "\n"
        "Output columns (TSV; schema_version=" SCHEMA_VERSION "):\n"
        "  locus_id group n_seqs n_codons_total n_codons_used\n"
        "  S N total_sd total_nd pi_S pi_N pi_NS_ratio\n"
        "  n_4fold_sites n_0fold_sites total_diffs_4fold total_diffs_0fold\n"
        "  pi_4fold pi_0fold pi_04_ratio\n"
        "  coverage_factor pi_S_corrected pi_N_corrected pi_4fold_corrected pi_0fold_corrected\n"
        "  [if --bootstrap > 0:] bootstrap_reps\n"
        "                        pi_S_lo pi_S_hi pi_N_lo pi_N_hi pi_NS_lo pi_NS_hi\n"
        "                        pi_4fold_lo pi_4fold_hi pi_0fold_lo pi_0fold_hi pi_04_lo pi_04_hi\n"
        "\n"
        "π4-fold = pairwise nt diversity at codon positions where the reference codon's\n"
        "          degeneracy is 3 (all 3 alt nts synonymous). ≈ neutral rate proxy.\n"
        "π0-fold = pairwise nt diversity at codon positions where degeneracy = 0\n"
        "          (any change is nonsynonymous). ≈ conservative πN proxy.\n"
        "π0/π4   = direct functional-burden ratio (parallel to πN/πS).\n");
}

static int load_fasta_list(const char* path, LocusEntry** out, int* n_out) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[pi_NS] Cannot open --fasta_list %s\n", path); return 0; }
    int cap = 16, n = 0;
    LocusEntry* arr = (LocusEntry*)malloc(cap * sizeof(LocusEntry));
    char* line = NULL; size_t lcap = 0; ssize_t got;
    int header_seen = 0;
    while ((got = getline(&line, &lcap, f)) != -1) {
        if (got > 0 && line[got-1] == '\n') line[--got] = 0;
        if (got > 0 && line[got-1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        char id[256] = "", p[512] = "";
        if (sscanf(line, "%255s %511s", id, p) < 2) continue;
        if (!header_seen && !strcasecmp(id, "locus_id")) { header_seen = 1; continue; }
        header_seen = 1;
        if (n == cap) { cap *= 2; arr = (LocusEntry*)realloc(arr, cap * sizeof(LocusEntry)); }
        strncpy(arr[n].locus_id, id, sizeof(arr[n].locus_id) - 1);
        arr[n].locus_id[sizeof(arr[n].locus_id) - 1] = 0;
        strncpy(arr[n].fasta_path, p, sizeof(arr[n].fasta_path) - 1);
        arr[n].fasta_path[sizeof(arr[n].fasta_path) - 1] = 0;
        n++;
    }
    free(line);
    fclose(f);
    *out = arr; *n_out = n;
    return n;
}

static const char* basename_of(const char* p) {
    const char* s = strrchr(p, '/');
    return s ? s + 1 : p;
}

#include <time.h>

int main(int argc, char** argv) {
    const char *fasta = NULL, *fasta_list = NULL, *locus_id = NULL,
               *groups_path = NULL, *out_path = NULL, *coverage_tsv = NULL;
    int ncores = 1, no_all = 0, n_bootstrap = 0;
    double coverage_factor_global = 1.0;
    unsigned int seed = (unsigned int)time(NULL);
    int seed_set = 0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--fasta")            && i+1<argc) fasta = argv[++i];
        else if (!strcmp(argv[i], "--fasta_list")       && i+1<argc) fasta_list = argv[++i];
        else if (!strcmp(argv[i], "--locus_id")         && i+1<argc) locus_id = argv[++i];
        else if (!strcmp(argv[i], "--groups")           && i+1<argc) groups_path = argv[++i];
        else if (!strcmp(argv[i], "--no_all"))                       no_all = 1;
        else if (!strcmp(argv[i], "--coverage_factor")  && i+1<argc) coverage_factor_global = atof(argv[++i]);
        else if (!strcmp(argv[i], "--coverage_tsv")     && i+1<argc) coverage_tsv = argv[++i];
        else if (!strcmp(argv[i], "--bootstrap")        && i+1<argc) n_bootstrap = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed")             && i+1<argc) { seed = (unsigned int)atoi(argv[++i]); seed_set = 1; }
        else if (!strcmp(argv[i], "--out")              && i+1<argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "--ncores")           && i+1<argc) ncores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
        else { fprintf(stderr, "[pi_NS] Unknown arg: %s\n", argv[i]); print_usage(); return 1; }
    }

    if (!fasta && !fasta_list) { print_usage(); return 1; }

    #ifdef _OPENMP
    omp_set_num_threads(ncores > 0 ? ncores : 1);
    #else
    (void)ncores;
    #endif

    if (n_bootstrap > 0) {
        fprintf(stderr, "[pi_NS] bootstrap=%d, seed=%u%s\n",
                n_bootstrap, seed, seed_set ? "" : " (time-based)");
    }

    init_codon_tables();

    NameLabel* group_map = NULL;
    int group_map_n = 0;
    if (groups_path) {
        load_group_file(groups_path, &group_map, &group_map_n);
        fprintf(stderr, "[pi_NS] --groups: %d sample→group entries\n", group_map_n);
    }

    LocusCoverage* cov_arr = NULL;
    int cov_n = 0;
    if (coverage_tsv) {
        load_coverage_tsv(coverage_tsv, &cov_arr, &cov_n);
        fprintf(stderr, "[pi_NS] --coverage_tsv: %d per-locus factors\n", cov_n);
    }

    FILE* fout = out_path ? fopen(out_path, "w") : stdout;
    if (!fout) { fprintf(stderr, "[pi_NS] Cannot open --out %s\n", out_path); return 2; }

    fprintf(fout, "# schema_version=%s method=NG86_pi", SCHEMA_VERSION);
    if (coverage_tsv || coverage_factor_global != 1.0)
        fprintf(fout, " coverage_global=%.6g%s", coverage_factor_global, coverage_tsv ? " coverage_tsv=yes" : "");
    if (n_bootstrap > 0) fprintf(fout, " bootstrap_reps=%d seed=%u", n_bootstrap, seed);
    fputc('\n', fout);
    fprintf(fout, "locus_id\tgroup\tn_seqs\tn_codons_total\tn_codons_used\t"
                  "S\tN\ttotal_sd\ttotal_nd\tpi_S\tpi_N\tpi_NS_ratio\t"
                  "n_4fold_sites\tn_0fold_sites\ttotal_diffs_4fold\ttotal_diffs_0fold\t"
                  "pi_4fold\tpi_0fold\tpi_04_ratio\t"
                  "coverage_factor\tpi_S_corrected\tpi_N_corrected\t"
                  "pi_4fold_corrected\tpi_0fold_corrected");
    if (n_bootstrap > 0)
        fprintf(fout, "\tbootstrap_reps\tpi_S_lo\tpi_S_hi\tpi_N_lo\tpi_N_hi\tpi_NS_lo\tpi_NS_hi"
                       "\tpi_4fold_lo\tpi_4fold_hi\tpi_0fold_lo\tpi_0fold_hi\tpi_04_lo\tpi_04_hi");
    fputc('\n', fout);

    unsigned int rng_state = seed;

    if (fasta) {
        char lid[256];
        if (locus_id) { strncpy(lid, locus_id, sizeof(lid) - 1); lid[sizeof(lid) - 1] = 0; }
        else {
            strncpy(lid, basename_of(fasta), sizeof(lid) - 1); lid[sizeof(lid) - 1] = 0;
            char* dot = strrchr(lid, '.'); if (dot) *dot = 0;
        }
        double cov = (cov_n > 0) ? lookup_coverage(cov_arr, cov_n, lid) : coverage_factor_global;
        process_locus(lid, fasta, group_map, group_map_n, no_all,
                      cov, n_bootstrap, &rng_state, fout);
    }
    if (fasta_list) {
        LocusEntry* loci = NULL;
        int n_loci = 0;
        load_fasta_list(fasta_list, &loci, &n_loci);
        fprintf(stderr, "[pi_NS] --fasta_list: %d loci\n", n_loci);
        for (int i = 0; i < n_loci; i++) {
            double cov = (cov_n > 0)
                         ? lookup_coverage(cov_arr, cov_n, loci[i].locus_id)
                         : coverage_factor_global;
            process_locus(loci[i].locus_id, loci[i].fasta_path,
                          group_map, group_map_n, no_all,
                          cov, n_bootstrap, &rng_state, fout);
        }
        free(loci);
    }

    if (out_path) fclose(fout);
    free(group_map);
    free(cov_arr);
    return 0;
}
