// =============================================================================
// codon_stats.c — 4-fold degenerate sites + Nei-Gojobori (1986) dN/dS, Ka/Ks.
//
// Pure C, no deps beyond libc. Built for "click a gene → see evolution" speed:
// precomputed 64×64 substitution tables, single linear pass per CDS pair,
// optional OpenMP across pairs.
//
// Input (one of):
//   --pairs <pairs.tsv>   TSV with columns: pair_id<TAB>seqA<TAB>seqB
//                         (sequences inline, same length, aligned, uppercase
//                         ACGT and '-' for gaps; 'N' or other → that codon skipped)
//   --fasta <align.fa>    Multi-FASTA where consecutive sequences are paired
//                         (>id_A then >id_B). Pair id = first header w/o suffix.
//
// Output: TSV to stdout (or --out <file>) with columns described in
//   engines/schemas/codon_stats.output.schema.json.
//
// Currently implements method=NG86 (Nei & Gojobori 1986). YN00 (Yang & Nielsen
// 2000, κ-corrected) is a planned follow-up; --method yn00 errors for now.
//
// Compile: gcc -O3 -march=native -fopenmp -o codon_stats codon_stats.c -lm
//
// Citation:
//   Nei M, Gojobori T (1986) Mol Biol Evol 3:418-426
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

#define SCHEMA_VERSION "codon_stats_v1"
#define N_CODONS 64

// ── Genetic code (NCBI table 1, standard). codon = (nt1<<4)|(nt2<<2)|nt3, A=0 C=1 G=2 T=3
//                AAx       ACx       AGx       ATx       CAx       CCx       CGx       CTx
//                GAx       GCx       GGx       GTx       TAx       TCx       TGx       TTx
static const char STD_CODE[N_CODONS + 1] =
    "KNKN" "TTTT" "RSRS" "IIMI"
    "QHQH" "PPPP" "RRRR" "LLLL"
    "EDED" "AAAA" "GGGG" "VVVV"
    "*Y*Y" "SSSS" "*CWC" "LFLF";

static char  aa_of[N_CODONS];
static int   is_stop[N_CODONS];
static double S_per_codon[N_CODONS];  // synonymous-site count (0..3)
static int   is_4d_third[N_CODONS];   // 1 if 3rd position is 4-fold degenerate
static double sd_table[N_CODONS][N_CODONS];  // NG86 path-averaged syn substitutions
static double nd_table[N_CODONS][N_CODONS];  // NG86 path-averaged nonsyn substitutions

static inline int nt_enc(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1;
    }
}

static void init_tables(void) {
    for (int c = 0; c < N_CODONS; c++) {
        aa_of[c] = STD_CODE[c];
        is_stop[c] = (STD_CODE[c] == '*');
    }

    for (int c = 0; c < N_CODONS; c++) {
        if (is_stop[c]) { S_per_codon[c] = 0; continue; }
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
        S_per_codon[c] = S;
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

    static const int PERM2[2][2] = {{0,1},{1,0}};
    static const int PERM3[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};

    for (int c1 = 0; c1 < N_CODONS; c1++) {
        for (int c2 = 0; c2 < N_CODONS; c2++) {
            sd_table[c1][c2] = 0;
            nd_table[c1][c2] = 0;
            if (c1 == c2 || is_stop[c1] || is_stop[c2]) continue;
            int n1[3] = { (c1 >> 4) & 3, (c1 >> 2) & 3, c1 & 3 };
            int n2[3] = { (c2 >> 4) & 3, (c2 >> 2) & 3, c2 & 3 };
            int diff[3], nd = 0;
            for (int p = 0; p < 3; p++) if (n1[p] != n2[p]) diff[nd++] = p;
            if (nd == 1) {
                if (aa_of[c1] == aa_of[c2]) sd_table[c1][c2] = 1;
                else                        nd_table[c1][c2] = 1;
                continue;
            }
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
            if (n_paths > 0) {
                sd_table[c1][c2] = sum_s / n_paths;
                nd_table[c1][c2] = sum_n / n_paths;
            }
        }
    }
}

// ── Per-pair stats ──────────────────────────────────────────────────────────

typedef struct {
    int  n_aligned_bp;
    int  n_codons_total;
    int  n_codons_used;
    int  n_codons_gap;
    int  n_codons_stop;
    int  n_4d_sites;
    int  n_4d_diffs;
    double S, N;            // synonymous, nonsynonymous sites (NG86)
    double sd, nd;          // synonymous, nonsynonymous diffs
    double pS, pN;          // raw proportions
    double dS, dN;          // Jukes-Cantor corrected
    double omega;           // dN/dS
    double p4d;             // raw 4D divergence
    double d4d;             // JC-corrected 4D divergence
} PairStats;

static double jc_correct(double p) {
    if (p < 0) return NAN;
    if (p >= 0.75) return INFINITY;
    return -0.75 * log(1.0 - (4.0 / 3.0) * p);
}

static void compute_pair(const char* sa, const char* sb, int len, PairStats* out) {
    memset(out, 0, sizeof(*out));
    if (len < 3) return;
    out->n_aligned_bp = len;
    out->n_codons_total = len / 3;
    for (int i = 0; i + 2 < len; i += 3) {
        int n1[3], n2[3];
        int gap = 0;
        for (int p = 0; p < 3; p++) {
            n1[p] = nt_enc(sa[i+p]);
            n2[p] = nt_enc(sb[i+p]);
            if (n1[p] < 0 || n2[p] < 0) { gap = 1; break; }
        }
        if (gap) { out->n_codons_gap++; continue; }
        int c1 = (n1[0] << 4) | (n1[1] << 2) | n1[2];
        int c2 = (n2[0] << 4) | (n2[1] << 2) | n2[2];
        if (is_stop[c1] || is_stop[c2]) { out->n_codons_stop++; continue; }
        out->n_codons_used++;
        out->S  += (S_per_codon[c1] + S_per_codon[c2]) / 2.0;
        out->sd += sd_table[c1][c2];
        out->nd += nd_table[c1][c2];
        if (is_4d_third[c1] && is_4d_third[c2] && n1[0] == n2[0] && n1[1] == n2[1]) {
            out->n_4d_sites++;
            if (n1[2] != n2[2]) out->n_4d_diffs++;
        }
    }
    out->N = 3.0 * out->n_codons_used - out->S;
    out->pS = (out->S > 0) ? out->sd / out->S : NAN;
    out->pN = (out->N > 0) ? out->nd / out->N : NAN;
    out->dS = jc_correct(out->pS);
    out->dN = jc_correct(out->pN);
    out->omega = (isfinite(out->dS) && out->dS > 0) ? out->dN / out->dS : NAN;
    out->p4d = (out->n_4d_sites > 0) ? (double)out->n_4d_diffs / out->n_4d_sites : NAN;
    out->d4d = jc_correct(out->p4d);
}

// ── I/O helpers ─────────────────────────────────────────────────────────────

typedef struct {
    char* id;
    char* a;
    char* b;
    int   len;   // length of a (must equal length of b)
} Pair;

typedef struct {
    Pair* items;
    int   n;
    int   cap;
} PairList;

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
    int nl = 0, header_skipped = 0;
    while ((got = getline(&line, &cap, f)) != -1) {
        nl++;
        if (got > 0 && line[got - 1] == '\n') line[--got] = 0;
        if (got > 0 && line[got - 1] == '\r') line[--got] = 0;
        if (got == 0 || line[0] == '#') continue;
        // Skip header if it starts with "pair_id"
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

    // Parse records.
    typedef struct { char* name; char* seq; int slen; } Rec;
    Rec* recs = NULL;
    int nr = 0, cap = 0;

    size_t i = 0;
    while (i < flen) {
        while (i < flen && (buf[i] == '\n' || buf[i] == '\r' || buf[i] == ' ' || buf[i] == '\t')) i++;
        if (i >= flen) break;
        if (buf[i] != '>') { i++; continue; }
        i++;
        // Read header until newline (use first whitespace-delimited token as name)
        size_t hs = i;
        while (i < flen && buf[i] != '\n' && buf[i] != '\r') i++;
        size_t he = i;
        size_t name_end = hs;
        while (name_end < he && buf[name_end] != ' ' && buf[name_end] != '\t') name_end++;
        char* name = (char*)malloc(name_end - hs + 1);
        memcpy(name, buf + hs, name_end - hs); name[name_end - hs] = 0;
        // Skip newline
        while (i < flen && (buf[i] == '\n' || buf[i] == '\r')) i++;
        // Read sequence until next '>'
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

    // Pair consecutive entries; pair_id = first name, strip trailing _A/_B/_a/_b/_1/_2 suffix
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

static void print_header(FILE* fout) {
    fprintf(fout, "# schema_version=%s method=NG86\n", SCHEMA_VERSION);
    fprintf(fout,
        "pair_id\tn_aligned_bp\tn_codons_total\tn_codons_used\tn_codons_gap\tn_codons_stop\t"
        "n_4d_sites\tn_4d_diffs\tS\tN\tsd\tnd\tpS\tpN\tdS\tdN\tomega\tp4d\td4d\n");
}

static void emit_field_double(FILE* fout, double v) {
    if (isnan(v))      fprintf(fout, "\tNA");
    else if (isinf(v)) fprintf(fout, "\tInf");
    else               fprintf(fout, "\t%.6f", v);
}

static void emit_row(FILE* fout, const char* pair_id, const PairStats* s) {
    fprintf(fout, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f",
            pair_id, s->n_aligned_bp, s->n_codons_total, s->n_codons_used,
            s->n_codons_gap, s->n_codons_stop, s->n_4d_sites, s->n_4d_diffs,
            s->S, s->N, s->sd, s->nd);
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
        "codon_stats — 4D sites + Nei-Gojobori dN/dS/Ka/Ks (fast, pure C)\n"
        "\n"
        "Input (one of):\n"
        "  --pairs <pairs.tsv>      TSV: pair_id<TAB>seqA<TAB>seqB (sequences inline)\n"
        "  --fasta <align.fa>       Multi-FASTA, consecutive records form pairs\n"
        "                            (e.g. >gene_A then >gene_B; suffix stripped)\n"
        "\n"
        "Options:\n"
        "  --method ng86            Nei-Gojobori 1986 (default; only one currently)\n"
        "  --method yn00            Yang-Nielsen 2000 κ-corrected (planned, errors for now)\n"
        "  --out <file>             Output TSV (default: stdout)\n"
        "  --ncores N               OpenMP threads across pairs (default: 1)\n"
        "  --no_header              Skip header line\n"
        "  -h, --help               This help\n"
        "\n"
        "Output columns (TSV, JSON schema in engines/schemas/codon_stats.output.schema.json):\n"
        "  pair_id n_aligned_bp n_codons_total n_codons_used n_codons_gap n_codons_stop\n"
        "  n_4d_sites n_4d_diffs S N sd nd pS pN dS dN omega p4d d4d\n");
}

int main(int argc, char** argv) {
    const char *pairs_path = NULL, *fasta_path = NULL, *out_path = NULL;
    const char *method = "ng86";
    int ncores = 1, no_header = 0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--pairs") && i + 1 < argc)   pairs_path = argv[++i];
        else if (!strcmp(argv[i], "--fasta") && i + 1 < argc)   fasta_path = argv[++i];
        else if (!strcmp(argv[i], "--method") && i + 1 < argc)  method = argv[++i];
        else if (!strcmp(argv[i], "--out") && i + 1 < argc)     out_path = argv[++i];
        else if (!strcmp(argv[i], "--ncores") && i + 1 < argc)  ncores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--no_header"))               no_header = 1;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { print_usage(); return 0; }
    }
    if (!pairs_path && !fasta_path) { print_usage(); return 1; }
    if (strcasecmp(method, "ng86") != 0) {
        fprintf(stderr, "[codon_stats] method '%s' not implemented yet (only ng86). Aborting.\n", method);
        return 2;
    }

    #ifdef _OPENMP
    omp_set_num_threads(ncores > 0 ? ncores : 1);
    #else
    (void)ncores;
    #endif

    init_tables();

    PairList pairs = {0};
    if (pairs_path) read_pairs_tsv(pairs_path, &pairs);
    if (fasta_path) read_pairs_fasta(fasta_path, &pairs);
    fprintf(stderr, "[codon_stats] %d pair(s) loaded\n", pairs.n);
    if (pairs.n == 0) return 0;

    PairStats* results = (PairStats*)calloc(pairs.n, sizeof(PairStats));

    #pragma omp parallel for schedule(static)
    for (int k = 0; k < pairs.n; k++) {
        compute_pair(pairs.items[k].a, pairs.items[k].b, pairs.items[k].len, &results[k]);
    }

    FILE* fout = out_path ? fopen(out_path, "w") : stdout;
    if (!fout) { fprintf(stderr, "[codon_stats] Cannot open %s for writing\n", out_path); return 3; }
    if (!no_header) print_header(fout);
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
