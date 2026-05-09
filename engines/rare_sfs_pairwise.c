// =============================================================================
// rare_sfs_pairwise.c — Pairwise rare allele sharing from BEAGLE dosage
//
// For each site, determines the minor allele count (MAC) across all samples.
// Then for sites in each MAC bin (doubleton=2, tripleton=3, quadrupleton=4, ...),
// counts how many times each pair of samples both carry ≥1 copy of the minor
// allele. Outputs a symmetric N×N matrix per MAC bin.
//
// This is the Bergström et al. (2020) Science Fig. 2 computation:
// "pairwise counts of doubleton alleles (alleles observed exactly twice
// across the dataset) among all 929 individuals, grouped by population."
//
// Input:
//   --beagle <file.beagle.gz>     BEAGLE GL file
//   --sample_list <file>          Sample IDs (BAM list order)
//   --chr <chr>                   Filter chromosome (optional)
//   --groups <file>               Group file: sample_id<tab>group_id
//   --bins 2,3,4                  MAC bins (default: 2,3,4,5)
//   --range START:END             Genomic range filter
//   --het_threshold <f>           Dosage threshold for "carries allele" (default: 0.5)
//   --out_prefix <prefix>         Output prefix (writes <prefix>.bin<N>.tsv per bin)
//
// Dosage logic:
//   Expected dosage d_i = GL1 + 2*GL2 (after normalization).
//   A sample "carries the minor allele" if d_i >= het_threshold.
//   MAC = count of carriers across all samples.
//   For MAC == bin, increment sharing matrix M[i][j] for all carrier pairs.
//
// Output per bin:
//   TSV with header row = sample IDs, then N rows of counts.
//   Also writes <prefix>.group_summary.bin<N>.tsv with group-level aggregates.
//
// Memory: O(N^2) for the sharing matrix + O(sites) for site data.
//   For N=226, matrix = 226*226*4 bytes = ~200KB per bin. Trivial.
//
// Compile: gcc -O3 -march=native -o rare_sfs_pairwise rare_sfs_pairwise.c -lz -lm
//
// Citation: Bergström A et al. (2020) Science 367:eaay5012
// =============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

#define MAX_IND 500
#define MAX_BINS 20
#define MAX_SITES 5000000

// ── GZ reader (same as region_popstats) ──

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

// ── Sample loading ──

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

// ── Group loading ──

typedef struct {
    char group[64];
    int group_idx;
} SampleGroup;

int load_groups(const char* path, char sample_names[][64], int n_samples,
                SampleGroup* sg, char group_names[][64]) {
    FILE* f = fopen(path, "r");
    if (!f) return 0;

    // Initialize all to "ungrouped"
    for (int i = 0; i < n_samples; i++) {
        sg[i].group_idx = -1;
        sg[i].group[0] = 0;
    }

    int n_groups = 0;
    char line[512];
    while (fgets(line, sizeof(line), f)) {
        int len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) len--;
        line[len] = 0;
        if (len == 0) continue;

        // Parse: sample_id\tgroup_id
        char sid[128] = "", gid[128] = "";
        char* tab = strchr(line, '\t');
        if (!tab) continue;
        *tab = 0;
        strncpy(sid, line, 127);
        strncpy(gid, tab + 1, 127);

        // Strip BAM extensions from sid
        char* ext;
        char* slash = strrchr(sid, '/');
        if (slash) memmove(sid, slash + 1, strlen(slash));
        if ((ext = strstr(sid, ".sorted.markdup.bam"))) *ext = 0;
        else if ((ext = strstr(sid, ".bam"))) *ext = 0;

        // Find sample index
        for (int i = 0; i < n_samples; i++) {
            if (!strcmp(sample_names[i], sid)) {
                // Find or create group
                int gi = -1;
                for (int g = 0; g < n_groups; g++) {
                    if (!strcmp(group_names[g], gid)) { gi = g; break; }
                }
                if (gi < 0 && n_groups < 50) {
                    gi = n_groups;
                    strncpy(group_names[gi], gid, 63);
                    group_names[gi][63] = 0;
                    n_groups++;
                }
                sg[i].group_idx = gi;
                strncpy(sg[i].group, gid, 63);
                break;
            }
        }
    }
    fclose(f);
    return n_groups;
}

// ── Main ──

int main(int argc, char** argv) {
    const char *beagle_path = NULL, *sample_path = NULL, *groups_path = NULL;
    const char *filter_chr = NULL, *out_prefix = "rare_sfs";
    int bins[MAX_BINS], n_bins = 0;
    double het_threshold = 0.5;
    int range_start = -1, range_end = -1;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--beagle") && i+1<argc)      beagle_path = argv[++i];
        else if (!strcmp(argv[i], "--sample_list") && i+1<argc) sample_path = argv[++i];
        else if (!strcmp(argv[i], "--groups") && i+1<argc)      groups_path = argv[++i];
        else if (!strcmp(argv[i], "--chr") && i+1<argc)         filter_chr = argv[++i];
        else if (!strcmp(argv[i], "--out_prefix") && i+1<argc)  out_prefix = argv[++i];
        else if (!strcmp(argv[i], "--het_threshold") && i+1<argc) het_threshold = atof(argv[++i]);
        else if (!strcmp(argv[i], "--range") && i+1<argc)       sscanf(argv[++i], "%d:%d", &range_start, &range_end);
        else if (!strcmp(argv[i], "--bins") && i+1<argc) {
            char buf[256];
            strncpy(buf, argv[++i], 255); buf[255] = 0;
            char* tok = strtok(buf, ",");
            while (tok && n_bins < MAX_BINS) {
                bins[n_bins++] = atoi(tok);
                tok = strtok(NULL, ",");
            }
        }
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            fprintf(stderr,
                "rare_sfs_pairwise — Pairwise rare allele sharing from BEAGLE dosage\n\n"
                "  --beagle <f>           BEAGLE GL file (.beagle.gz)\n"
                "  --sample_list <f>      Sample IDs (BAM list order)\n"
                "  --groups <f>           Group file: sample_id<tab>group_id\n"
                "  --chr <chr>            Filter chromosome\n"
                "  --bins 2,3,4,5         MAC bins (default: 2,3,4,5)\n"
                "  --het_threshold <f>    Dosage threshold for carrier (default: 0.5)\n"
                "  --range START:END      Genomic range filter\n"
                "  --out_prefix <prefix>  Output prefix\n");
            return 0;
        }
    }

    if (!beagle_path || !sample_path) {
        fprintf(stderr, "Need --beagle and --sample_list\n");
        return 1;
    }

    // Default bins
    if (n_bins == 0) {
        bins[0] = 2; bins[1] = 3; bins[2] = 4; bins[3] = 5;
        n_bins = 4;
    }

    // Load samples
    static char names[MAX_IND][64];
    int n_ind = load_samples(sample_path, names);
    fprintf(stderr, "[rare_sfs] %d samples loaded\n", n_ind);

    // Load groups (optional)
    static SampleGroup sg[MAX_IND];
    static char group_names[50][64];
    int n_groups = 0;
    if (groups_path) {
        n_groups = load_groups(groups_path, names, n_ind, sg, group_names);
        fprintf(stderr, "[rare_sfs] %d groups loaded\n", n_groups);
    }

    // Allocate sharing matrices: one per bin, each n_ind × n_ind
    // Use flat arrays for cache friendliness
    int** sharing = (int**)malloc(n_bins * sizeof(int*));
    long* site_counts = (long*)calloc(n_bins, sizeof(long));  // sites per bin
    for (int b = 0; b < n_bins; b++) {
        sharing[b] = (int*)calloc(n_ind * n_ind, sizeof(int));
    }

    // Carrier buffer
    int* carriers = (int*)malloc(n_ind * sizeof(int));

    // Read BEAGLE and process
    Gz gz;
    if (!gz_open(&gz, beagle_path)) {
        fprintf(stderr, "Cannot open %s\n", beagle_path);
        return 1;
    }

    char* line = (char*)malloc(1 << 20);
    gz_getline(&gz, line, 1 << 20); // skip header

    long n_sites_total = 0, n_sites_used = 0;

    while (gz_getline(&gz, line, 1 << 20)) {
        n_sites_total++;

        // Parse marker for chr and position
        char* tab1 = strchr(line, '\t');
        if (!tab1) continue;
        *tab1 = 0;
        char* us = strrchr(line, '_');
        if (!us) continue;
        *us = 0;
        char* chr = line;
        int pos = atoi(us + 1);

        if (filter_chr && strcmp(chr, filter_chr) != 0) continue;
        if (range_start >= 0 && (pos < range_start || pos > range_end)) continue;

        // Skip allele1, allele2
        char* p = tab1 + 1;
        char* tab2 = strchr(p, '\t'); if (!tab2) continue;
        char* tab3 = strchr(tab2 + 1, '\t'); if (!tab3) continue;
        p = tab3 + 1;

        // Parse dosages and determine carriers
        int mac = 0;
        int n_carriers = 0;
        for (int i = 0; i < n_ind; i++) {
            double gl0 = 0, gl1 = 0, gl2 = 0;
            char* endp;
            gl0 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;
            gl1 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;
            gl2 = strtod(p, &endp); p = endp; while (*p == '\t' || *p == ' ') p++;

            double s = gl0 + gl1 + gl2;
            if (s > 1e-15) { gl0 /= s; gl1 /= s; gl2 /= s; }
            double dos = gl1 + 2.0 * gl2;

            if (dos >= het_threshold) {
                carriers[n_carriers++] = i;
                mac++;
            }
        }

        // Check which bin this site falls in
        for (int b = 0; b < n_bins; b++) {
            if (mac == bins[b]) {
                site_counts[b]++;
                // Increment sharing for all carrier pairs
                for (int ci = 0; ci < n_carriers; ci++) {
                    int ii = carriers[ci];
                    // Self-sharing (diagonal)
                    sharing[b][ii * n_ind + ii]++;
                    for (int cj = ci + 1; cj < n_carriers; cj++) {
                        int jj = carriers[cj];
                        sharing[b][ii * n_ind + jj]++;
                        sharing[b][jj * n_ind + ii]++;
                    }
                }
                break; // each site in exactly one bin
            }
        }
        n_sites_used++;
    }

    gz_close(&gz);
    free(line);
    free(carriers);

    fprintf(stderr, "[rare_sfs] %ld sites total, %ld after chr/range filter\n",
            n_sites_total, n_sites_used);

    // ── Write output per bin ──
    for (int b = 0; b < n_bins; b++) {
        char path[512];
        snprintf(path, sizeof(path), "%s.bin%d.tsv", out_prefix, bins[b]);
        FILE* f = fopen(path, "w");
        if (!f) { fprintf(stderr, "Cannot write %s\n", path); continue; }

        fprintf(stderr, "[rare_sfs] Bin MAC=%d: %ld sites\n", bins[b], site_counts[b]);

        // Header with sample names
        fprintf(f, "sample");
        for (int j = 0; j < n_ind; j++) fprintf(f, "\t%s", names[j]);
        fprintf(f, "\n");

        // Matrix rows
        for (int i = 0; i < n_ind; i++) {
            fprintf(f, "%s", names[i]);
            for (int j = 0; j < n_ind; j++) {
                fprintf(f, "\t%d", sharing[b][i * n_ind + j]);
            }
            fprintf(f, "\n");
        }
        fclose(f);

        // ── Group-level summary ──
        if (n_groups > 0) {
            snprintf(path, sizeof(path), "%s.group_summary.bin%d.tsv", out_prefix, bins[b]);
            f = fopen(path, "w");
            if (!f) continue;

            // Group × group sum
            fprintf(f, "group");
            for (int g2 = 0; g2 < n_groups; g2++) fprintf(f, "\t%s", group_names[g2]);
            fprintf(f, "\n");

            for (int g1 = 0; g1 < n_groups; g1++) {
                fprintf(f, "%s", group_names[g1]);
                for (int g2 = 0; g2 < n_groups; g2++) {
                    long sum = 0;
                    int count = 0;
                    for (int i = 0; i < n_ind; i++) {
                        if (sg[i].group_idx != g1) continue;
                        for (int j = 0; j < n_ind; j++) {
                            if (sg[j].group_idx != g2) continue;
                            if (i == j) continue; // skip diagonal
                            sum += sharing[b][i * n_ind + j];
                            count++;
                        }
                    }
                    // Mean per pair
                    double mean_val = count > 0 ? (double)sum / count : 0;
                    fprintf(f, "\t%.2f", mean_val);
                }
                fprintf(f, "\n");
            }
            fclose(f);
        }
    }

    // Clean up
    for (int b = 0; b < n_bins; b++) free(sharing[b]);
    free(sharing);
    free(site_counts);

    fprintf(stderr, "[rare_sfs] Done\n");
    return 0;
}
