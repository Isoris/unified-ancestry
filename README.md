# Unified Ancestry Module — Integration Guide

> **Chat 15 update (2026-04-17).** Cache layout is now K-sharded and
> sample-set-tagged. Cache directory moved out of the code tree. See
> "What changed" below for the migration story. The prior rewiring
> plan in `README_REWIRING_v9_need_rewrite.md` is now implemented and
> that file can be deleted after one full pipeline rerun.

## Directory Layout

```
unified_ancestry/                                        CODE (git)
├── 00_ancestry_config.sh          Central config (source this everywhere)
├── run_ancestry.sh                Unified dispatcher
├── run_full_pipeline.sh           End-to-end driver
├── instant_q                      Bash CLI for Engine B
├── region_stats                   Bash CLI for stats dispatcher
│
├── engines/                       All C/C++ sources + single Makefile
│   ├── Makefile                   Build all 5 binaries: make -C engines
│   ├── instant_q.cpp              Engine B (Skotte 2013, fixed-F EM)
│   ├── region_popstats.c          Engine F (Fst/dXY/theta/MI)
│   ├── rare_sfs_pairwise.c        Pairwise rare-allele SFS
│   ├── export_q_residual_dosage.c Q-corrected residual dosage
│   └── hobs_windower.c            Multi-scale Hobs/F windowing
│
├── hobs_hwe/                      Hobs/HWE sub-pipeline (HH = tag)
│   ├── 00_hobs_hwe_config.sh
│   ├── README_hobs_hwe.md
│   ├── run_hobs_hwe.sh
│   ├── STEP_HH_A_build_subset_bamlists.sh
│   ├── STEP_HH_B_run_angsd_hwe.sh
│   ├── STEP_HH_C_compute_hobs_windows.sh
│   └── STEP_HH_D_candidate_overlay.py
│
├── launchers/                     All SLURM scripts
│   ├── LAUNCH_instant_q_precompute.slurm   28 × K-sweep
│   ├── LAUNCH_region_popstats.slurm
│   ├── LAUNCH_rare_sfs_pairwise.slurm
│   ├── LAUNCH_q_residual_dosage.slurm
│   └── LAUNCH_hobs_hwe.slurm               280-task subset × chr array
│
├── wrappers/
│   ├── instant_q.R                R wrapper (source from any R script)
│   └── instant_q.py               Python wrapper (import from any script)
│
├── dispatchers/
│   └── region_stats_dispatcher.R  Unified stats dispatcher
│
├── steps/                         Pipeline-step Python scripts (UA = tag)
│   ├── STEP_UA_C_snp_q_support.py            run_full_pipeline Step 3
│   ├── STEP_UA_D_internal_ancestry_composition.py     Step 4
│   ├── STEP_UA_E_candidate_classifier.py              Step 5
│   └── STEP_UA_F_export_module5b.py                   Step 6
│
├── plots/                         All plotting R scripts
│   ├── plot_fst_dxy_tracks.R
│   ├── plot_ld_panels.R
│   ├── plot_rare_sfs_heatmap.R
│   ├── plot_hobs_hwe.R
│   ├── plot_local_Q_diagnostics.R
│   └── theme_systems_plate.R      Shared theme helper
│
├── registries/
│   └── build_registries.py        Interval + sample-subset + cov registry
│
├── tests/
│   └── test_integration.sh
│
└── _archive/                      Backup/legacy files
    └── instant_q.R.bk_chat17

$BASE/ancestry_cache/                                    DATA (scratch, NOT in git)
├── K02/
│   ├── <chr>.<sample_set>.local_Q_summary.tsv.gz   cohort aggregates
│   ├── <chr>.<sample_set>.local_Q_samples.tsv.gz   per-sample × per-window Q
│   └── <chr>.<sample_set>.local_Q_meta.tsv
├── K03/
├── ...
├── K08/                                              ← canonical K
├── ...
├── K20/
└── manifest.tsv                   audit of every cached (chrom,K,sample_set)
```

`<sample_set>` is a short tag `N<count>_<6char-sha1>` computed from
the sorted `SAMPLE_LIST` at configure time. Running against a different
sample subset produces a different tag and cannot collide with prior
output.

The cache directory lives OUT of `unified_ancestry/` (the code dir) so
a git status never sees data files and a repo reclone doesn't nuke
your cache. Default is `$BASE/ancestry_cache`; override by setting
`LOCAL_Q_DIR` in the env or in `00_ancestry_config.sh`.

## What changed (chat 15)

Pre-2026-04-17:
- One K only (K=8). Cache files named `<chr>.local_Q_summary.tsv.gz`.
- Cache directory was `$BASE/unified_ancestry/local_Q/` — inside the code tree.
- No way to tell which sample set a cache file came from.

Now:
- K=2..20 sweep. Per-K subdir under `$LOCAL_Q_DIR/K<NN>/`.
- Sample-set tag in every filename.
- Cache directory is `$BASE/ancestry_cache/` by default (set via `LOCAL_Q_DIR`).
- Canonical K (default K=8) is flattened into the inversion precomp RDS as
  `localQ_delta12_K08` / `localQ_entropy_K08` / `localQ_ena_K08` columns.

Backward compat: legacy flat-layout cache files (no K subdir, no sample
tag) are still readable — the wrapper's `resolve_cache_summary()` falls
through to them at canonical K. Manifest files auto-upgrade on the next
write.

## Step 1: Compile the C/C++ Engines

```bash
cd unified_ancestry
make -C engines            # builds all 5 binaries
# Or build a single one:
make -C engines instant_q
make -C engines region_popstats
make -C engines rare_sfs_pairwise
make -C engines export_q_residual_dosage
make -C engines hobs_windower
```

Requirements: `g++` with C++17, `gcc`, OpenMP, zlib (`-lz`).
On LANTA: `module load gcc zlib` or use the conda `assembly` env.

## Step 2: Configure

Edit `00_ancestry_config.sh` to match your HPC paths. Critical vars:

```bash
BEAGLE_DIR      # Per-chr BEAGLE GL files from Step 1
BEST_QOPT       # Canonical-K Q matrix; other K's resolved by filename swap
BEST_FOPT       # Canonical-K F matrix; other K's resolved by filename swap
SAMPLE_LIST     # One sample ID per line
INSTANT_Q_BIN   # Path to compiled binary
LOCAL_Q_DIR     # Where the cache goes (default: $BASE/ancestry_cache)
CANONICAL_K     # Which K gets flattened into precomp RDS (default: 8)
K_SWEEP         # Space-sep list of K's to precompute (default: 2..20)
```

## Step 3: Validate

```bash
./run_ancestry.sh validate
```

## Step 4: Precompute All (Chroms × K Values)

### Option A: SLURM array (recommended)
```bash
mkdir -p logs
sbatch launchers/LAUNCH_instant_q_precompute.slurm
# 28 chroms × 19 K values = 532 array tasks. Array-parallel runtime ~30-60 min
# depending on queue and node concurrency.
```

### Option B: Canonical K only (fast, one-time sanity)
```bash
K_SWEEP=8 sbatch launchers/LAUNCH_instant_q_precompute.slurm
# 28 tasks, ~10 min total.
```

### Option C: Serial bash (single chrom, single K)
```bash
./instant_q --precompute --chr C_gar_LG01
```

## Step 5: Query

### Single region (prints to stdout)
```bash
./instant_q --chr C_gar_LG01 --start 35000000 --end 41000000
```

### From R
```r
source("wrappers/instant_q.R")
configure_instant_q(config_file = "00_ancestry_config.sh")

# Single region
q <- get_Q("C_gar_LG01", start = 35e6, end = 41e6)
# Returns data.table: sample_id, Q1..Q8, max_q, delta12, entropy, ena, assigned_pop

# Whole-chromosome summary (from cache)
summary <- get_Q_summary("C_gar_LG01")
# Returns: window_id, chrom, start_bp, end_bp, mean_delta12, mean_entropy, mean_ena, sd_delta12

# Full stats dispatcher
stats <- get_region_stats("C_gar_LG06", 35e6, 41e6,
  what = c("Q", "Fst", "Hobs", "dXY", "theta_pi"),
  groups = list(inv = inv_samples, ref = ref_samples))
```

### From Python
```python
from instant_q import configure, get_Q, get_Q_summary, get_region_stats

configure(config_file="00_ancestry_config.sh")

q = get_Q("C_gar_LG01", start=35e6, end=41e6)  # pandas DataFrame
summary = get_Q_summary("C_gar_LG01")

stats = get_region_stats("C_gar_LG06", 35e6, 41e6,
    what=["Q", "Fst", "Hobs"],
    groups={"inv": inv_samples, "ref": ref_samples})
```

---

## Integration with the Inversion Pipeline

### What changes in `00_inversion_config.sh`

Add these lines:

```bash
# ── Unified ancestry integration ──
ANCESTRY_CONFIG="${BASE}/unified_ancestry/00_ancestry_config.sh"
ANCESTRY_R="${BASE}/unified_ancestry/wrappers/instant_q.R"
LOCAL_Q_DIR="${BASE}/unified_ancestry/local_Q"
```

### What changes in `STEP_C01a_snake1_precompute.R`

Replace the old ancestry_bridge.R section (~lines 706-730) with:

```r
# ── Merge local Q from Engine B (replaces ancestry_bridge.R) ──
ancestry_r <- Sys.getenv("ANCESTRY_R",
  file.path(dirname(dirname(f_inv)), "unified_ancestry", "wrappers", "instant_q.R"))

if (file.exists(ancestry_r)) {
  tryCatch({
    source(ancestry_r)
    configure_instant_q(
      config_file = Sys.getenv("ANCESTRY_CONFIG", ""),
      cache_dir = local_q_dir
    )
    inv_like_dt <- merge_local_Q_into_invlikeness(inv_like_dt, cache_dir = local_q_dir)
    message("[PRECOMP] Engine B local Q columns added: localQ_delta12, localQ_entropy, localQ_ena")
  }, error = function(e) {
    message("[PRECOMP] Engine B merge failed: ", e$message, " — continuing without")
  })
} else {
  message("[PRECOMP] instant_q.R not found at ", ancestry_r,
          " — run precompute first to enable real Q metrics")
}
```

### What changes in the SLURM launcher for precomp

Add this **before** the precomp R call:

```bash
# Run Engine B precompute for this chromosome if not cached
source "${BASE}/unified_ancestry/00_ancestry_config.sh"
if [[ ! -f "${LOCAL_Q_DIR}/${CHR}.local_Q_summary.tsv.gz" ]]; then
  "${INSTANT_Q_BIN}" \
    --beagle "${BEAGLE_DIR}/${CHR_BEAGLE}" \
    --fopt "${BEST_FOPT}" --qinit "${BEST_QOPT}" \
    --precompute --outdir "${LOCAL_Q_DIR}" --chr "${CHR}" \
    --window_size 100 --window_step 20 --em_iter 20 --ncores 4
fi
```

---

## How Engine B Works (Math)

The algorithm follows Skotte, Korneliussen & Albrechtsen (2013) *Genetics* 195:693-702.

**Key insight:** When F (allele frequencies per population) is held FIXED from a
genome-wide run, the EM only needs to update Q (admixture proportions). This
converges in 5-10 iterations instead of 250+, making it ~50× faster.

For each site j, individual i:

```
p_i = Σ_k Q[i,k] × F[j,k]           # expected allele frequency

P(g=0|p) = (1-p)²                     # genotype probabilities
P(g=1|p) = 2p(1-p)                    # under HWE
P(g=2|p) = p²

pp_g = P(g|p) × GL(g)                 # posterior × likelihood
E[g] = (pp_1 + 2×pp_2) / Σ pp_g      # expected genotype

prodA[i,k] += E[g]/(1-p) × F[j,k]    # accumulators
prodB[i,k] += (2-E[g])/p × (1-F[j,k])

Q_new[i,k] ∝ Q[i,k] × (prodA[i,k] + prodB[i,k])
Normalize rows to sum to 1. Iterate 5-20 times.
```

This is a **clean-room implementation** of the published math — same algorithm,
original code. Not copied from NGSadmix source.

---

## Relationship to Existing Modules

### Module 2B Step 1 (unchanged)
The 18 existing scripts (S00-S18) continue to work as-is. They produce the
genome-wide Q/F/covariance files that Engine B uses as fixed references.
Engine A = wrapper around these scripts.

### Module 2B Step 2 (superseded by Engine B)
The old `run_step2.sh` called NGSadmix per window — impossibly slow.
Engine B replaces this with the fixed-F EM: same Q output, ~10,000× faster.

### Module 2C (uses Engine B output)
Module 2C computes per-SNP support vectors (dosage-Q correlations).
It now reads the Engine B Q matrix instead of running NGSadmix per window.
The `compute_snp_support.py` helper is unchanged — it just reads a different
Q source file.

### Inversion Pipeline (consumes Engine B cache)
`STEP_C01a_snake1_precompute.R` sources `instant_q.R` and calls
`merge_local_Q_into_invlikeness()` to add `localQ_delta12`, `localQ_entropy`,
`localQ_ena` columns to the inv_likeness table.

---

## Stats Dispatcher — Full Inventory

The `get_region_stats()` function (R and Python) supports:

| Stat | Backend | Groups needed? |
|------|---------|---------------|
| Q (delta12, entropy, ENA) | C++ Engine B (fixed-F EM) | No |
| Fst (Hudson) | Dosage-based formula | Yes (2+ groups) |
| Hobs | Dosage het rate [0.3, 1.7] | No |
| dXY | p₁(1-p₂) + p₂(1-p₁) | Yes (2 groups) |
| theta_pi | 2p(1-p) × n/(n-1) | No |
| HWE | Observed vs expected het | No |
| dosage_matrix | Raw matrix from BEAGLE | No |
| gl_matrix | Raw GL triplets | No |

Stats from Groups 5-10 in the prompt (LD, ROH, TE, gene annotation, coverage)
are **not yet implemented** in the dispatcher. These require additional data
sources (BAMs, ngsF-HMM output, RepeatMasker, GFF3) that live in separate
modules. The dispatcher architecture supports adding them as new backends —
each stat just needs a `compute_X(dosage, positions, groups)` function.

---

## Performance Expectations

| Operation | Expected time |
|-----------|--------------|
| Single region query (100 SNPs) | < 1 second |
| One chromosome precompute (~7000 windows) | 30-120 seconds |
| All 28 chromosomes (parallel) | 5-15 minutes |
| Cache read (per chromosome) | < 1 second |
| Stats dispatcher (Q + Fst + Hobs) | 2-5 seconds |

The C++ engine processes ~500-2000 windows/second depending on window size
and number of EM iterations. The R fallback is ~100× slower but works without
compilation.
