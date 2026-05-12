# ENGINES.md — inversion-popgen-toolkit engine inventory

*Last updated: 2026-05-12*
*Status key: **STABLE** = field-tested, **WIP** = in development, **NEW** = just-added, not yet field-tested, **ORPHAN** = maybe abandoned*

---

## Purpose of this document

The toolkit has several compiled engines and companion scripts under
`unified_ancestry/`. Some are stable and documented, some are
half-built, and I haven't kept a clear mental model of which is which.
This file is the single source of truth: for each engine, one paragraph
covering **what it does · where it lives · inputs · outputs · status**.

Whenever I touch an engine, I update this file. If I don't know something,
I write `???` — that's a signal that I need to look, not guess.

---

## Quick summary table

| # | Engine | Path | Language | Status |
|---|---|---|---|---|
| 1 | `instant_q` | `unified_ancestry/engines/instant_q` | C++ | STABLE (+ schema v2 anchoring metadata, 2026-05) |
| 2 | `region_popstats` | `unified_ancestry/engines/region_popstats` | C | STABLE (+ schema v3: karyotype HWE, per-SNP F_IS Wilcoxon, 2026-05) |
| 3 | `hobs_windower` | `unified_ancestry/engines/hobs_windower` | C | WIP (+ `--groups` per-arrangement Hobs, 2026-05) |
| 4 | `rare_sfs_pairwise` | `unified_ancestry/engines/rare_sfs_pairwise` | C | WIP |
| 5 | `export_q_residual_dosage` | `unified_ancestry/engines/export_q_residual_dosage` | C | WIP |
| 6 | `codon_stats` | `unified_ancestry/engines/codon_stats` | C | NEW (2026-05) |
| 7 | `pi_NS` | `unified_ancestry/engines/pi_NS` | C | NEW (2026-05) |
| 8 | `fdM` | `unified_ancestry/engines/fdM` | C | NEW (2026-05) |
| 9 | `region_test` | `unified_ancestry/engines/region_test` | C | NEW (2026-05) |
| 10 | `ultrabootstrap` | `unified_ancestry/engines/ultrabootstrap` | C | NEW (2026-05) |
| 11 | `fst_outlier_scan` | `unified_ancestry/engines/fst_outlier_scan` | C | NEW (2026-05) |
| 12 | `homolog_index_build` | `unified_ancestry/engines/homolog_index_build` | C | NEW (2026-05) |
| 13 | `homolog_index_query` | `unified_ancestry/engines/homolog_index_query` | C | NEW (2026-05) |
| 14 | `homolog_atlas_server` | `unified_ancestry/engines/homolog_atlas_server` | C | NEW (2026-05) |

All sources live flat in `unified_ancestry/engines/` alongside a single
`Makefile` that builds every target. Build: `make -C engines [target]`.

JSON schemas (input/output contracts) live in `unified_ancestry/engines/schemas/`:
`codon_stats.{input,output}.schema.json`, `homolog_index.{binary,query.output}.schema.json`,
`homolog_atlas_server.api.schema.json`, `homolog_atlas_config.schema.json`.

*The stable engines (1–2) are the ones used in the LG28 analysis
(2026-04-20). Engines 6–14 were added 2026-05 to support the broader
"Diversity Atlas + Inversion Atlas" plan (πN/πS, π0/π4 fold, fdM
introgression, FST outlier scans, homolog/paralog atlas, generic
bootstrap/permutation testing). NEW = passes smoke tests in this
session but not yet run on real cohort data.*

---

## 1. Engine B — `instant_q`  ·  STABLE

**One-line purpose**: Local ancestry inference — produces per-window,
per-sample ancestry assignments (maxQ label, max_q score, per-K posteriors)
across a chromosome from a thinned BEAGLE + a reference NGSadmix F matrix.

**Location**:
- Binary: `unified_ancestry/engines/instant_q`
- Source: `unified_ancestry/engines/instant_q.cpp`
- R wrapper: `unified_ancestry/wrappers/instant_q.R`
- Python wrapper: `unified_ancestry/wrappers/instant_q.py`

**Inputs**:
- BEAGLE-format genotype likelihoods (thin or dense)
- NGSadmix F matrix (`*.fopt.gz`) for a given K
- Optional: sample list, chromosome filter

**Outputs** (under `unified_ancestry/local_Q/`):
- `<chr>.local_Q_summary.tsv.gz`  (per-window ancestry summary: delta12, entropy, ENA)
- `<chr>.local_Q_samples.tsv.gz`  (per-window-per-sample maxQ label + score)
- `<chr>.local_Q_meta.tsv.gz`     (per-window positions / SNP counts)

**Used by**:
- `phase_qc_shelf/STEP_Q06_precompute.sh` (multi-K × multi-scale caches)
- `phase_qc_shelf/STEP_Q06_ancestry_tracks.sh`
- `phase_qc_shelf/STEP_Q06_multiscale.sh`

**Known quirks**:
- `00_ancestry_config.sh` clobbers `INSTANT_Q_BIN` and `LOCAL_Q_DIR` when
  sourced. Fix: stash-and-restore before/after sourcing (done in Q06 precompute
  v2).

---

## 2. Engine F — `region_popstats`  ·  STABLE

**One-line purpose**: Cohort-wide and between-group popgen statistics over
sliding windows — theta_pi, theta_W, Tajima's D, Hp, Hudson Fst, dXY, dA, MI,
MInorm. Computes these for all-samples, per-group (Hom1/Het/Hom2), and
between-group pairs in a single pass.

**Location**:
- Binary: `unified_ancestry/engines/region_popstats`
- Source: `unified_ancestry/engines/region_popstats.c`
- R dispatcher: `unified_ancestry/dispatchers/region_stats_dispatcher.R`
- R track plotter: `unified_ancestry/plots/plot_fst_dxy_tracks.R`

**Inputs**:
- BEAGLE (dense)
- Position file (`*.pos.fixed`)
- Sample list
- Per-group sample files (e.g. one line per sample per group → `groups_invgt/`)
- Window size + step

**Outputs**:
- `popstats_<grouping>.<chr>.tsv`  with 29 columns: window_id, chrom, start, end,
  n_sites, n_sites_used, S, theta_pi_all, theta_W_all, Tajima_D, Hp,
  theta_pi_{Het,Hom1,Hom2}, Fst/dXY/dA/MI/MInorm_{Het_Hom1, Het_Hom2, Hom1_Hom2}

**Used by**:
- `phase_qc_shelf/STEP_Q07_popstats.sh` (main consumer)
- `phase_qc_shelf/R/q04_compose_plot.R` (reads Fst/pi panels from popstats TSV)

**Known quirks**:
- Writes **exact zeros** in rows where `n_sites_used < ~20`. These render as
  fake drop-to-baseline spikes in Fst if you don't filter. Q04 filters on
  `n_sites_used >= Q04_MIN_SITES` (default 20). Document this in the paper
  methods so reviewers don't ask.

---

## 3. Engine H — `hobs_windower`  ·  WIP

**One-line purpose**: Site-level observed heterozygosity (Hobs) + multi-scale
windowed aggregation. Wraps a patched ANGSD (`-doHWE 1` from Isoris fork) to
get per-site Hobs, then windows it. Used for secondary confirmation of
inversion breakpoints (Hobs elevation in Het group + depression in both Hom
groups is the classic inversion signature).

**Location**:
- C binary source: `unified_ancestry/engines/hobs_windower.c`
- Compiled binary: `unified_ancestry/engines/hobs_windower`
- Driver: `unified_ancestry/hobs_hwe/run_hobs_hwe.sh`
- Config: `unified_ancestry/hobs_hwe/00_hobs_hwe_config.sh`
- Step scripts: `unified_ancestry/hobs_hwe/0[1-4]_*.sh|.py`
- Plot script: `unified_ancestry/plots/plot_hobs_hwe.R`
- Phase_qc_shelf wrappers: `STEP_Q07b_hobs_per_group.sh`, `STEP_Q07c_hobs_windower.sh`
- Merge: `R/q07c_merge_hobs.R`

**Inputs**:
- Per-group BAM lists (01_build_subset_bamlists.sh)
- Reference FASTA
- Site list (optional)

**Outputs** (per group, per chromosome):
- `*.hwe.gz` (raw ANGSD output with Hobs column)
- Windowed Hobs tracks at multiple scales

**Used by**:
- `phase_qc_shelf/STEP_Q07b_*` and `STEP_Q07c_*`

**Status notes / WIP flags**:
- ??? — is the patched ANGSD binary path correctly set in every launcher?
- ??? — does the merge step (q07c_merge_hobs.R) handle missing groups gracefully?
- ??? — has Engine H been run on LG28? If yes, where are the outputs?
- Not yet wired into Q04 diagnostic PDF (should add a Hobs panel per invgt group)

---

## 4. `rare_sfs_pairwise`  ·  WIP

**One-line purpose**: Bergström 2020-style pairwise rare-allele sharing matrix
— per minor-allele-count bin (doubletons, tripletons, ...), counts how many
times each pair of samples both carry ≥1 copy of the minor allele. Signal of
cryptic relatedness / fine-scale population structure.

**Location**:
- Source: `unified_ancestry/engines/rare_sfs_pairwise.c`
- Binary: `unified_ancestry/engines/rare_sfs_pairwise`
- Plotter: `unified_ancestry/plots/plot_rare_sfs_heatmap.R`

**Inputs**:
- `--beagle <file.beagle.gz>` BEAGLE GL file
- `--sample_list <file>` sample IDs (BAM-list order)
- `--groups <file>` sample_id<TAB>group_id (for group-level aggregates)
- `--bins 2,3,4,5` MAC bins (default 2,3,4,5)
- `--het_threshold F` dosage threshold for "carries allele" (default 0.5)
- `--chr <chr>`, `--range START:END` filters

**Outputs**:
- `<prefix>.bin<N>.tsv` — symmetric N×N sample matrix per MAC bin
- `<prefix>.group_summary.bin<N>.tsv` — group-level aggregates

**Citation**: Bergström A et al. (2020) Science 367:eaay5012.

**Status notes**:
- NOT orphan — clear, well-documented purpose in the source header.
- Useful for the inversion paper IF a relatedness / structure section is
  needed. Otherwise downstream for population-structure QC.
- Not wired into phase_qc_shelf yet.

---

## 5. `export_q_residual_dosage`  ·  WIP

**One-line purpose**: Q-corrected residual dosage from BEAGLE — for each
sample × site, computes `observed − Σ_k Q_ik · 2 · f_jk`. The residual
removes the ancestry component for cohort-wide structure correction before
downstream association / drift tests.

**Location**:
- Source: `unified_ancestry/engines/export_q_residual_dosage.c`
- Binary: `unified_ancestry/engines/export_q_residual_dosage`

**Inputs**:
- `--beagle <file.beagle.gz>` BEAGLE GL file
- `--fopt <file.fopt.gz>` NGSadmix F matrix (n_sites × K)
- `--local_q <file.tsv.gz>` window-level local Q cache (output of `instant_q`)
  OR a single global Q file
- `--sample_list <file>`, `--chr <chr>`, `--window_size N` (default 100, must
  match the Engine B run that produced local_q)

**Outputs** (under `--out_prefix`):
- `*.residual.beagle.gz` — BEAGLE-format residual GLs for compatibility with
  downstream tools
- `*.residual.dosage.tsv.gz` — marker × sample matrix of residual dosages
- `*.sites.tsv.gz` — chrom, pos, marker

**Status notes**:
- NOT orphan — clear, well-documented purpose in the source header.
- Use case: pre-conditioning step before running anything sensitive to
  population structure on the cohort.
- Not yet wired into the inversion-paper pipeline; would slot in between
  Engine B (Q) and any selection / association scan.

---

## 6. `codon_stats`  ·  NEW (2026-05)

**One-line purpose**: Pairwise dN/dS / Ka/Ks + 4-fold degenerate sites from
aligned CDS pairs. NG86 (Nei-Gojobori 1986, JC distance) and YN00
(Yang-Nielsen 2000, κ-corrected K80 distance) are both supported in one binary.

**Location**:
- Source: `engines/codon_stats.c`
- Binary: `engines/codon_stats`
- Schemas: `engines/schemas/codon_stats.{input,output}.schema.json`

**Inputs**: `--pairs <tsv>` (pair_id<TAB>seqA<TAB>seqB) OR `--fasta <multi.fa>`
(consecutive records form pairs). `--method ng86|yn00`, `--kappa F`, `--ncores N`.

**Outputs**: TSV with 25 columns per pair (S, N, sd, nd, sd_ts/tv, nd_ts/tv,
n_4d_sites/diffs/ts/tv, pS, pN, dS, dN, omega, p4d, d4d). Schema v2.

**Status**: NEW. Passes hand-verified smoke tests on 4 synthetic pair cases
(identical / syn-only / mixed / gappy). Not yet run on real data.

---

## 7. `pi_NS`  ·  NEW (2026-05)

**One-line purpose**: Population-level synonymous (πS) and non-synonymous (πN)
diversity per CDS locus, plus the per-codon-position simpler proxies π4-fold
(4-fold degenerate sites, neutral) and π0-fold (0-fold degenerate sites,
always nonsyn). Per-group splits (e.g. HOM_STD / HET / HOM_INV) supported via
`--groups` for inversion-burden comparisons.

**Location**:
- Source: `engines/pi_NS.c`
- Binary: `engines/pi_NS`

**Inputs**: `--fasta <multi.fa>` (one locus) OR `--fasta_list <tsv>` (many),
optional `--groups seq→label.tsv`, `--coverage_factor F` or
`--coverage_tsv f.tsv`, `--bootstrap N` (codon resampling, 95% CI), `--seed K`.

**Outputs**: TSV with πS, πN, πN/πS, π4-fold, π0-fold, π0/π4, coverage-corrected
versions, and (with --bootstrap) lower/upper 2.5%/97.5% percentiles. Schema v2.

**Atlas placement**:
- Diversity Atlas → "Functional burden / selection efficacy" panel (main).
- Inversion Atlas → per-candidate consequence panel (with `--groups karyo.tsv`).

**Status**: NEW. Numbers hand-verified on a 4-hap × 7-codon synthetic.

---

## 8. `fdM`  ·  NEW (2026-05)

**One-line purpose**: Patterson's D + Malinsky 2015 fdM introgression statistic
per genomic window from BEAGLE GLs. Genealogy `((P1, (P2, P3)), O)`. fdM > 0 →
P3 → P2 gene flow; fdM < 0 → P3 → P1. Optional block-jackknife for genome-wide
SE / Z / p.

**Location**:
- Source: `engines/fdM.c`  (in-source Patterson's-D explainer header)
- Binary: `engines/fdM`

**Inputs**: `--beagle <f.beagle.gz>`, `--sample_list`, `--pops
P1:f,P2:f,P3:f,O:f`, `--chr`, `--windows <bed>` or `--fixed_win W:S`,
`--jackknife_blocks N`.

**Outputs**: per-window TSV (D, fdM, sum_ABBA, sum_BABA, sum_num, sum_denom,
plus jackknife_SE / Z / p in a final GENOME row when --jackknife_blocks is set).

**Use case**: works for cross-species AND within-species (P1=family A,
P2=family B, P3=candidate donor family, O=distant outgroup family). Useful
inside-vs-flanking inversion candidates to test "does this inversion carry
introgressed ancestry?"

**Status**: NEW. Sign-convention verified on deliberate ABBA-rich and
BABA-rich synthetic data (caught a sign bug in the denom branch during
testing, fixed before commit).

---

## 9. `region_test`  ·  NEW (2026-05)

**One-line purpose**: Generic Wilcoxon + permutation test for arbitrary TSV
value columns, partitioned by a region column. Use for "region vs rest of
genome" comparisons on pi_NS / fdM / region_popstats outputs.

**Location**:
- Source: `engines/region_test.c`
- Binary: `engines/region_test`

**Inputs**: `--input <tsv>`, `--region_col <name>`, `--value_cols a,b,c`,
optional `--rest_label <name>`, `--skip_regions a,b`, `--permutations N`
(default 1000), `--seed K`, `--ncores N`.

**Outputs**: per (region × value_col): n_in, n_out, mean_in/out, median_in/out,
effect_diff, W, z, wilcoxon_p, perm_p.

**Status**: NEW. Smoke-tested on 13-locus synthetic: regA (high πS) → p=0.003
(wilcoxon) / 0.001 (permutation); regB (matches rest) → p=1.0.

**Driver**: paper Bonhomme et al. methods excerpt — "Wilcoxon tests or
resampling procedures (1000 tests) to compare diversity estimates between
candidate regions and the rest of the genome."

---

## 10. `ultrabootstrap`  ·  NEW (2026-05)

**One-line purpose**: Generic bootstrap CI for arbitrary TSV columns. Row mode
(default) or block mode (`--block_col` → resample whole blocks of rows as a
unit; keeps within-block correlation for codons-within-locus / sites-within-chrom).
Stats: mean / median / sum / sd. Optional per-group output.

**Location**:
- Source: `engines/ultrabootstrap.c`
- Binary: `engines/ultrabootstrap`

**Inputs**: `--input <tsv>`, `--value_cols a,b,c`, `--statistic mean[,median,sum,sd]`,
`--group_col`, `--block_col`, `--n_boot N` (default 1000), `--ci_lo/hi`,
`--seed K`, `--ncores N`.

**Outputs**: per (group × value_col × statistic): n, observed, boot_mean,
boot_sd, boot_lo, boot_hi.

**Status**: NEW. Row + block bootstrap modes verified on synthetic. For
multi-FASTA inputs, pre-process into a TSV with locus_id as block_col.

---

## 11. `fst_outlier_scan`  ·  NEW (2026-05)

**One-line purpose**: Genome-wide FST outlier scan, cichlid-paper style.
Threshold + gap-aware merge of per-window FST into OUTLIER regions; HDR
escalation if region overlaps a stronger 10-kb-window FST entry (or by
region max_FST fallback). Distinct from inversion-candidate-confirmation:
this scans the whole genome de novo for small differentiated regions.

**Location**:
- Source: `engines/fst_outlier_scan.c`
- Binary: `engines/fst_outlier_scan`

**Inputs**: `--input <fst_windows.tsv>` (auto-detects chrom/start/end/fst
columns; can pass `--fst_col Fst_HOM_A_HOM_B`), `--fst_threshold` or
`--top_pct P`, `--merge_gap_bp` (default 10000), `--stronger_threshold`
(default 0.30), optional `--stronger_windows <tsv>` for HDR check, optional
BED annotations (`--genes_tsv`, `--inversion_bed`, `--breakpoint_bed`,
`--te_bed`).

**Outputs**: `out_windows.tsv` (flagged input rows), `out_regions.tsv`
(merged), `out_hdr.tsv` (HDR-class only), `out_summary.json`.

**Caveat baked into the binary**: phrased as "empirical FST outliers" /
"highly differentiated regions" — NOT "selected regions". No neutral
calibration here.

**Status**: NEW. Full smoke-test passes (16-window synthetic with stronger
10-kb windows + inversion BED + genes TSV).

---

## 12. `homolog_index_build`  ·  NEW (2026-05)

**One-line purpose**: Builds a sorted, mmap-friendly flat-binary homolog
index from DIAMOND tabular (`-outfmt 6`) and/or miniprot PAF files. Manifest
TSV input also supported. Powers the "click a gene → see all
chains/paralogs/orthologs instantly" atlas use case.

**Location**:
- Source: `engines/homolog_index_build.c`
- Header: `engines/homolog_index.h` (shared with query + server)
- Binary: `engines/homolog_index_build`
- Schema: `engines/schemas/homolog_index.binary.schema.json`

**Inputs**: any combination of `--diamond <path>[=<label>]`,
`--miniprot <path>[=<label>]`, `--manifest <tsv>`; `--out <file.holindx>`.

**Outputs**: single binary `.holindx` file. Layout (little-endian, packed):
IndexHeader → GeneRecord[sorted by name] → HitRecord[clustered by gene,
bitscore desc] → string pool.

**Status**: NEW. End-to-end tested with synthetic DIAMOND + miniprot inputs.

---

## 13. `homolog_index_query`  ·  NEW (2026-05)

**One-line purpose**: Microsecond gene→hits lookup over a `.holindx`. mmap +
binary search. CLI only; the HTTP server (#14) is the typical atlas hook.

**Location**:
- Source: `engines/homolog_index_query.c`
- Binary: `engines/homolog_index_query`
- Schema: `engines/schemas/homolog_index.query.output.schema.json`

**Inputs**: `--index <f.holindx>`, `--gene <id>`, `--format tsv|json`,
`--source diamond|miniprot`, `--min_identity F`, `--min_bitscore F`, `--limit N`.

**Outputs**: stdout. TSV (default) or JSON with hits ordered by bitscore desc.

**Status**: NEW. Tested.

---

## 14. `homolog_atlas_server`  ·  NEW (2026-05)

**One-line purpose**: Tiny single-binary HTTP server over a `.holindx`. mmaps
the index once, serves `/lookup` + `/health` + `/` as JSON with CORS open.
Pure C + libc + POSIX sockets + pthread; thread-per-connection.

**Location**:
- Source: `engines/homolog_atlas_server.c`
- Binary: `engines/homolog_atlas_server`
- Schema: `engines/schemas/homolog_atlas_server.api.schema.json`

**Inputs**: `--index <f.holindx>`, `--port N` (default 8765), `--bind <addr>`
(default 127.0.0.1).

**Endpoints**: `GET /lookup?gene=<id>[&source=...][&min_identity=F]
[&min_bitscore=F][&limit=N]`; `GET /health`; `GET /`.

**Status**: NEW. End-to-end tested with curl against a fresh index — `/health`,
filtered lookups, missing-gene (`found:false`), CORS preflight all behave.

**Atlas integration**:
- `dispatchers/build_homolog_atlas.sh` orchestrates miniprot (strict/chain/loose
  presets) + optional DIAMOND blastp + index build in one command.
- `dispatchers/run_miniprot_chain.sh` is the lower-level single-target miniprot
  wrapper.
- `dispatchers/homolog_atlas.R` is the R wrapper exposing
  `atlas_lookup()`, `atlas_orthologs()`, `atlas_paralogs()`,
  `codon_stats_run()`, `atlas_build()`, `atlas_health_check()` to atlas pages.

---

## Related engines (not under `unified_ancestry/engines/` but worth flagging)

### `STEP_UA_C_snp_q_support.py`

**Location**: `unified_ancestry/steps/STEP_UA_C_snp_q_support.py`
**One-line (from header)**: Per-SNP soft support vectors to Q components.
For each SNP, computes Pearson |corr| of dosage vs each Q column, normalizes
to soft weights, classifies markers as `dominant_single_Q` /
`dominant_plus_secondary` / `mixed_Q_support` / `diffuse_Q_support` /
`low_information`. Aggregates into windows with coherence score, switch counts,
dominant Q.
**Inputs**: `--beagle`, `--q_cache_dir` (Engine B output) OR `--qopt` (global),
`--sample_list`, `--K`, `--chr`.
**Outputs** (in `--outdir`): `q_support_per_snp.tsv`, `q_support_windows.tsv`.
**Status**: WIP. Origin: adapted from `MODULE_2C/helpers/compute_snp_support.py +
summarize_windows.py`. Reads Q from Engine B cache (instant_q) instead of old
Module 2B tables.

### `STEP_UA_D_internal_ancestry_composition.py` (formerly `nested_composition.py`)

**Location**: `unified_ancestry/steps/STEP_UA_D_internal_ancestry_composition.py`
**One-line (from header)**: Generic interval-internal ancestry structure classifier.
**Status**: WIRED via Phase 7 wrapper `STEP_C01i_c_nested_composition.py`.
**Notes**: Renamed 2026-04-24 to match schema/block_type (`internal_ancestry_composition`). The engine is ALSO called directly as a CLI by `unified_ancestry/run_full_pipeline.sh` Step 4 (chromosome-level composition pass). See `docs/NESTED_VS_COMPOSITE.md` for the role split between the engine and the Phase 7 inversion-candidate wrapper.

### `STEP_UA_E_candidate_classifier.py`

**Location**: `unified_ancestry/steps/STEP_UA_E_candidate_classifier.py`
**One-line (from header)**: Classify inversion candidates from Q support
windows + Engine B per-sample Q. For each candidate, overlaps Q-support windows
and produces `coherent_q_loaded_core`, `likely_mixed_background`,
`likely_parasite_windows`, plus `ancestry_class` (`ancestry_driven` /
`structural` / `mixed`) and `internal_structure_type` (from
`STEP_UA_D` nested-composition output).
**Inputs**: `--candidates`, `--q_windows` (output of `STEP_UA_C`),
`--q_cache_dir`, `--nested` (output of `STEP_UA_D`), `--K`.
**Outputs**: `candidate_classification.tsv` under `--outdir`.
**Status**: WIP. Origin: `MODULE_2C/helpers/summarize_candidates.py +` new logic.
Important downstream of phase_9_classification.

### `STEP_UA_F_export_module5b.py`

**One-line (from header)**: Export marker/window annotation strips for
Module 5B (heatmap + regime track + candidate segmentation views).
**Location**: `unified_ancestry/steps/STEP_UA_F_export_module5b.py`
**Inputs**: `--snp_support` (output of `STEP_UA_C`),
`--window_support` (also from `STEP_UA_C`), `--K`.
**Outputs** (in `exports_for_module5b/`):
- `marker_strip_dominant_Q.tsv` — per-SNP dominant ancestry + colour
- `marker_strip_entropy.tsv` — per-SNP support entropy
- `marker_strip_confidence.tsv` — per-SNP class + confidence
- `window_strip_coherence.tsv` — per-window coherence + dominant Q
- `module5b_integration_table.tsv` — full per-SNP table for Module 5B
**Status**: WIP. Origin: `MODULE_2C/helpers/export_for_5b.py`.

---

## Registries — producer / consumer split (NOT competing)

Two halves of one registry system live in two repos:

### A. `inversion-popgen-toolkit/registries/api/` (READER side)

- `registry_loader.R` — "the unified entry point for the three registries"
  (sample / interval / results). Looked up by `dispatchers/region_stats_dispatcher.R`
  via a search path (`registries/api/R/registry_loader.R` first, with
  fallbacks).
- Has integration tests (`test_results_registry.py`,
  `test_interval_registry_extensions.R`, `test_sample_registry_extensions.R`).
- **Status**: partially working, tests exist but WIP.

### B. `unified_ancestry/registries/build_registries.py` (WRITER side)

- Single Python builder script. Origin: `MODULE_2B_Step2/helpers/make_intervals.py +
  make_subsets.py + make_cov_registry.py`.
- WRITES TO:
  - `registries/interval_registry.tsv` — all intervals at all scales
  - `registries/sample_subsets.tsv` — named sample subsets with paths
  - `registries/cov_registry.tsv` — PCAngsd covariance file registry
- Inputs: `--config 00_ancestry_config.sh`, optional `--candidates`.
- **Status**: builder for the TSVs A consumes. Resolved 2026-05-12: **A and B
  are producer/consumer, NOT parallel systems.** B writes the registry TSVs
  in this repo; A's loader reads them from the sister inversion-popgen-toolkit
  repo (via the dispatcher's search path).
- No further unification needed; one is the API surface (R/Python), the other
  is the population step.

---

## How I use this document

1. **Before I reach for an engine**: read its section to confirm status.
2. **After I touch an engine**: update the "Known quirks" or "Status notes".
3. **Before writing methods**: copy the one-line purpose straight into the manuscript.
4. **Periodic review** (monthly): re-walk `unified_ancestry/engines/` with
   `build_inventory.py` and compare to this file. Any new file that isn't
   listed here is a pending engine to document.

---

## What's NOT in this document (on purpose)

- **Phase modules** (`inversion_modules/phase_*`) — those are pipelines,
  not engines. They'll get their own `PHASES.md` if needed.
- **DELLY/Manta SV calling** (`Modules/MODULE_4*`) — established third-party
  tools, documented by their upstream projects.
- **Deprecated code** (`_archive/`, `_archive_superseded/`) — explicitly
  not engines I should reach for.

---

## Open questions for future-me

1. Should `rare_sfs_pairwise` and `export_q_residual_dosage` be promoted to
   "STABLE"? Both are NOT orphan (clear purpose, documented headers).
   Both have non-zero startup cost to wire into the inversion paper — decide
   at lab whether the relatedness/structure track (`rare_sfs_pairwise`) and
   the Q-residual pre-conditioning step (`export_q_residual_dosage`) earn
   their keep for this paper or get parked until next.
2. ~~Is registry system A or B canonical?~~ Resolved 2026-05-12. **They are
   producer/consumer.** `unified_ancestry/registries/build_registries.py`
   writes the three TSVs; `inversion-popgen-toolkit/registries/api/R/
   registry_loader.R` reads them. No conflict, no consolidation needed.
3. ~~Is there a sixth engine I haven't found?~~ Resolved 2026-05-12. **Yes —
   nine new ones** (codon_stats, pi_NS, fdM, region_test, ultrabootstrap,
   fst_outlier_scan, homolog_index_build/_query, homolog_atlas_server),
   all documented above as #6–14.
4. Should Engine H be added as a Q04 diagnostic panel? Still open — clean
   v2.2 feature for phase_qc_shelf when Engine H goes STABLE.
5. **New 2026-05**: should the two-mode architecture (candidate-vs-flanks
   AND genome-scan) be generalized beyond `fst_outlier_scan`? `pi`, `dXY`,
   `dA`, `Tajima's D` would all benefit from the same pattern. Plan: rename
   `fst_outlier_scan` → `outlier_scan` (statistic-agnostic, `--stat_col`,
   `--direction high|low|abs`) + add a sibling `candidate_vs_flanks` binary.
   Scope explicitly excludes LD / ROH / VESM for now (per 2026-05 decision).
6. **New 2026-05**: the 9 new engines all pass smoke tests on synthetic
   data but NONE have been run on the real cohort yet. Once any is run
   in anger, promote NEW → WIP, then to STABLE on second use.
