# ENGINES.md — inversion-popgen-toolkit engine inventory

*Last updated: 2026-04-20*
*Status key: **STABLE** = field-tested, **WIP** = in development, **ORPHAN** = maybe abandoned*

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
| 1 | `instant_q` | `unified_ancestry/engines/instant_q` | C++ | STABLE |
| 2 | `region_popstats` | `unified_ancestry/engines/region_popstats` | C | STABLE |
| 3 | `hobs_windower` | `unified_ancestry/engines/hobs_windower` | C | WIP |
| 4 | `rare_sfs_pairwise` | `unified_ancestry/engines/rare_sfs_pairwise` | C | WIP |
| 5 | `export_q_residual_dosage` | `unified_ancestry/engines/export_q_residual_dosage` | C | WIP |

All sources live flat in `unified_ancestry/engines/` alongside a single
`Makefile` that builds every target. Build: `make -C engines [target]`.

*The stable engines are the ones I actually used in the LG28 analysis
(2026-04-20). The others exist but haven't been wired into a complete
end-to-end pipeline yet.*

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

**One-line purpose**: ???  (from directory name: pairwise rare-allele site
frequency spectrum — measures rare allele sharing between sample pairs as a
signal of cryptic relatedness or population structure).

**Location**:
- Source: `unified_ancestry/engines/rare_sfs_pairwise.c`
- Plotter: `unified_ancestry/plots/plot_rare_sfs_heatmap.R`

**Inputs**: ???
**Outputs**: ??? (heatmap implied by plotter name)
**Used by**: not wired into phase_qc_shelf yet
**Status notes**:
- Binary compiled but I don't remember when I last ran it
- Is this needed for the inversion paper or is it for a separate analysis?
- **Action at lab**: read `plots/plot_rare_sfs_heatmap.R` header to find out what it expects

---

## 5. `export_q_residual_dosage`  ·  WIP

**One-line purpose**: ???  (from filename: residual allele dosage after
removing Q-inferred ancestry component. Useful for cohort-wide structure
correction before downstream tests).

**Location**:
- Source: `unified_ancestry/engines/export_q_residual_dosage.c`

**Inputs**: ???
**Outputs**: ???
**Used by**: ??? (no obvious downstream consumer in current pipeline)
**Status notes**:
- Compiled binary exists; don't remember using it in any analysis
- **Action at lab**: open the `.c` source, read the first 50 lines to
  understand args. If I still don't know why I wrote it, consider it ORPHAN
  and move to `_archive_superseded/` so it stops showing up in the engine list.

---

## Related engines (not under `unified_ancestry/engines/` but worth flagging)

### `STEP_UA_C_snp_q_support.py`

**Location**: `unified_ancestry/steps/STEP_UA_C_snp_q_support.py`
**Status**: ??? — header says "MARKER CLASSIFICATION" but I need to read more.
**Action**: re-read when reviewing the ancestry pipeline.

### `STEP_UA_D_internal_ancestry_composition.py` (formerly `nested_composition.py`)

**Location**: `unified_ancestry/steps/STEP_UA_D_internal_ancestry_composition.py`
**One-line (from header)**: Generic interval-internal ancestry structure classifier.
**Status**: WIRED via Phase 7 wrapper `STEP_C01i_c_nested_composition.py`.
**Notes**: Renamed 2026-04-24 to match schema/block_type (`internal_ancestry_composition`). The engine is ALSO called directly as a CLI by `unified_ancestry/run_full_pipeline.sh` Step 4 (chromosome-level composition pass). See `docs/NESTED_VS_COMPOSITE.md` for the role split between the engine and the Phase 7 inversion-candidate wrapper.

### `STEP_UA_E_candidate_classifier.py`

**Location**: `unified_ancestry/steps/STEP_UA_E_candidate_classifier.py`
**One-line (from header)**: Classify inversion candidates from ???.
**Status**: ??? — sounds like it could be important for phase_9_classification.

### `STEP_UA_F_export_module5b.py`

**Location**: `unified_ancestry/steps/STEP_UA_F_export_module5b.py`
**One-line (from header)**: Export marker/window annotation strips.
**Status**: ??? — references phase_5, presumably annotation export.

---

## Registries — related but separate

Two registry systems coexist (likely needs unification):

### A. `inversion-popgen-toolkit/registries/` (the "chat 16" generation)

- Three registries (sample, interval, results) + a loader in bash/R/python
- Has integration tests (`test_results_registry.py`, `test_interval_registry_extensions.R`, `test_sample_registry_extensions.R`)
- `registry_loader.R` is "the unified entry point for the three registries"
- **Status**: partially working, tests exist but WIP

### B. `unified_ancestry/registries/build_registries.py`

- A single Python builder script
- **Status**: ??? — is this the builder that feeds A's registries, or a parallel system?
- **Action at lab**: open this file, read first 30 lines, figure out what it emits and where it writes to. If it writes to the same tables A reads from, it's the builder → mark clearly. If it's parallel, decide which wins.

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
   "STABLE" after being wired into a real analysis, or deprecated if they
   don't find a use in the inversion paper?
2. Is registry system A or B canonical? Pick one and move forward.
3. Is there a sixth engine I haven't found yet? (If the inventory walker
   reports something new under `unified_ancestry/engines/` next month, add it
   here immediately.)
4. Should Engine H be added as a Q04 diagnostic panel? If yes, that's a
   clean v2.2 feature for phase_qc_shelf.
