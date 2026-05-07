#!/usr/bin/env bash
# =============================================================================
# 00_ancestry_config.sh — Central config for the unified ancestry module (v12.2)
#
# v12.2 REWIRED (repo reorg):
#   - All MODULE_* directories moved to Modules/ at repo root
#   - 12 files moved out of unified_ancestry/engines and utils to
#     Modules/MODULE_5A_inversion_discovery/ and Modules/MODULE_5B_inversion_followup/
#   - Added QRES_BIN, RARE_SFS_BIN, MODULE_5B_ENGINES_DIR for moved C engines
#   - Added SNP_Q_SUPPORT_PY, NESTED_COMPOSITION_PY, CANDIDATE_CLASSIFIER_PY,
#     EXPORT_MODULE5B_PY, HOBS_CANDIDATE_OVERLAY_PY, HOBS_PLOT_R for moved scripts
#
# v12.1 REWIRED:
#   - Added REGISTRY_DIR, LOAD_BRIDGE, SAMPLES_IND for cross-module wiring
#   - Added SAMPLE_MAP_R path
#   - init_registry_if_needed() helper for auto-initialization
#   - Palette unchanged
#
# Consumed by:
#   Engine A (full survey), Engine B (instant Q), Module 2C (SNP support),
#   Inversion pipeline (via instant_q.R), Stats dispatcher, Region popstats,
#   Hobs/HWE engine, pipeline_bridge.sh
# =============================================================================

# ── Project root ─────────────────────────────────────────────────────────────
export BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ── Reference genome ─────────────────────────────────────────────────────────
export REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
export REF_FAI="${REF}.fai"

# ── Programs ─────────────────────────────────────────────────────────────────
export RSCRIPT_BIN="${RSCRIPT_BIN:-/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript}"
export CONDA_ENV="assembly"

# ── Binaries (compile from src/) ─────────────────────────────────────────────
export INSTANT_Q_BIN="${BASE}/unified_ancestry/src/instant_q"
export POPSTATS_BIN="${BASE}/unified_ancestry/engines/fst_dxy/region_popstats"
export HOBS_WINDOWER_BIN="${BASE}/unified_ancestry/engines/hobs_hwe/scripts/hobs_windower"

# ── Binaries moved to Modules/MODULE_5B_inversion_followup/engines (v12.2) ──
# These C engines feed 5B R plots (plot_ld_panels.R, 07_plot_rare_sfs_heatmap.R)
# and were relocated from unified_ancestry/engines/fst_dxy/ to live next to
# their consumers. See Modules/MODULE_5B_inversion_followup/engines/Makefile.
export QRES_BIN="${BASE}/Modules/MODULE_5B_inversion_followup/engines/export_q_residual_dosage"
export RARE_SFS_BIN="${BASE}/Modules/MODULE_5B_inversion_followup/engines/rare_sfs_pairwise"
export MODULE_5B_ENGINES_DIR="${BASE}/Modules/MODULE_5B_inversion_followup/engines"

# ── Cross-module Python / R scripts (v12.2) ─────────────────────────────────
# Ancestry-driven candidate analysis/classification moved to MODULE_5A:
export SNP_Q_SUPPORT_PY="${BASE}/Modules/MODULE_5A_inversion_discovery/analysis/snp_q_support.py"
# Ancestry composition engine: single source of truth at the repo path below.
# (Previously pointed at MODULE_5A_inversion_discovery which no longer carries
#  this file; dedup pass 2026-04-24, see docs/NESTED_VS_COMPOSITE.md.)
# The engine file is also imported as a library by
# inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_c_nested_composition.py
# via ANCESTRY_ENGINE_DIR — this variable is the CLI-standalone path.
export NESTED_COMPOSITION_PY="${BASE}/unified_ancestry/engines/nested_composition/internal_ancestry_composition.py"
export CANDIDATE_CLASSIFIER_PY="${BASE}/Modules/MODULE_5A_inversion_discovery/classification/candidate_classifier.py"
# Inversion followup utilities / figures moved to MODULE_5B:
export EXPORT_MODULE5B_PY="${BASE}/Modules/MODULE_5B_inversion_followup/utils/export_module5b.py"
export HOBS_CANDIDATE_OVERLAY_PY="${BASE}/Modules/MODULE_5B_inversion_followup/analysis/04_candidate_overlay.py"
export HOBS_PLOT_R="${BASE}/Modules/MODULE_5B_inversion_followup/figures/05_plot_hobs_hwe.R"

# ── Cross-module wiring (v12.1) ─────────────────────────────────────────────
export SAMPLES_IND="${BASE}/het_roh/01_inputs_check/samples.ind"
export REGISTRY_DIR="${BASE}/sample_registry"
# LOAD_BRIDGE path: prefer the flattened layout, fall back to legacy v8.5.
# The bridge itself auto-detects its location and doesn't strictly need this
# env var; it's set so downstream scripts that just `source "$LOAD_BRIDGE"`
# pick it up correctly.
if [[ -z "${LOAD_BRIDGE:-}" ]]; then
  for _try in \
      "${BASE}/utils/load_bridge.R" \
      "${BASE}/inversion_modules/utils/load_bridge.R" \
      "${BASE}/inversion_codebase_v8.5/utils/load_bridge.R"; do
    if [[ -f "$_try" ]]; then
      export LOAD_BRIDGE="$_try"; break
    fi
  done
  unset _try
  : "${LOAD_BRIDGE:=${BASE}/utils/load_bridge.R}"  # default even if nothing found
fi
# sample_map.R / sample_registry.R co-live with load_bridge.R
_LB_DIR=$(dirname "${LOAD_BRIDGE}" 2>/dev/null || echo "${BASE}/utils")
export SAMPLE_MAP_R="${SAMPLE_MAP_R:-${_LB_DIR}/sample_map.R}"
export SAMPLE_REGISTRY_R="${SAMPLE_REGISTRY_R:-${_LB_DIR}/sample_registry.R}"
unset _LB_DIR

# ── Step 1 outputs (Engine A → Engine B) ─────────────────────────────────────
export STEP1_DIR="${BASE}/popstruct_thin"
export STEP1_BEAGLE_DIR="${STEP1_DIR}/04_beagle_byRF_majmin"
export STEP1_SITES_DIR="${STEP1_DIR}/03_sites"
export STEP1_SAMPLE_LIST="${STEP1_DIR}/list_of_samples_one_per_line_same_bamfile_list.tsv"
export STEP1_PRUNED_LIST="${STEP1_DIR}/06_relatedness/pruned_samples.txt"

# ── BEAGLE directory for Engine B ────────────────────────────────────────────
export BEAGLE_DIR="${STEP1_BEAGLE_DIR}"

# ── Best NGSadmix run ────────────────────────────────────────────────────────
export BEST_SEED_TABLE="${STEP1_DIR}/05_ngsadmix_global/runs_thin500/best_seed_by_K.tsv"

# Default K: AUTO-DETECT from best_seed_by_K.tsv
if [[ -f "${BEST_SEED_TABLE}" && -z "${DEFAULT_K:-}" ]]; then
  _AUTO_K=$(awk -F'\t' 'NR>1 {print $1, $3}' "${BEST_SEED_TABLE}" 2>/dev/null \
            | sort -k2 -rn | head -1 | awk '{print $1}')
  export DEFAULT_K="${_AUTO_K:-8}"
  unset _AUTO_K
else
  export DEFAULT_K="${DEFAULT_K:-8}"
fi

# Best Q and F
_BEST_SEED=""
if [[ -f "${BEST_SEED_TABLE}" ]]; then
  _BEST_SEED=$(awk -F'\t' -v k="${DEFAULT_K}" 'NR>1 && $1==k {print $2; exit}' "${BEST_SEED_TABLE}" 2>/dev/null)
fi
_BEST_SEED="${_BEST_SEED:-1}"
_K_PAD=$(printf "%02d" "${DEFAULT_K}")
export BEST_QOPT="${BEST_QOPT:-${STEP1_DIR}/05_ngsadmix_global/runs_thin500/thin500_K${_K_PAD}_seed${_BEST_SEED}.qopt}"
export BEST_FOPT="${BEST_FOPT:-${STEP1_DIR}/05_ngsadmix_global/runs_thin500/thin500_K${_K_PAD}_seed${_BEST_SEED}.fopt.gz}"
unset _BEST_SEED _K_PAD

# Aliases
export QINIT_FILE="${BEST_QOPT}"
export FOPT_FILE="${BEST_FOPT}"
export SAMPLE_LIST="${STEP1_SAMPLE_LIST}"

# ── Cache directory (Engine B output) ────────────────────────────────────────
# 2026-04-17: moved OUT of the code tree. Cache files are mutable run-state,
# code is versioned under $BASE/unified_ancestry/. Keeping them separate means
# a git status or a repo reclone doesn't see data files as untracked.
# Legacy flat layout <LOCAL_Q_DIR>/<chr>.local_Q_*.tsv.gz is tolerated on read;
# new writes go under <LOCAL_Q_DIR>/K<NN>/<chr>.local_Q_*.tsv.gz so the K
# sweep has one sharded cache per K level.
export LOCAL_Q_DIR="${LOCAL_Q_DIR:-${BASE}/ancestry_cache}"

# Canonical K — what the precomp RDS flattens in as localQ_*_K<NN> columns.
# Mirrors DEFAULT_K above; kept as a separate alias so the purpose at each
# reference site is obvious.
export CANONICAL_K="${CANONICAL_K:-${DEFAULT_K:-8}}"

# K-sweep — which K levels to precompute during the per-chrom precompute
# launcher. Default: full 2..20. Set K_SWEEP=8 to skip the sweep and only
# precompute canonical.
export K_SWEEP="${K_SWEEP:-2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20}"

# Sample group — which registered sample_registry group this run targets.
# Filenames + results_registry manifest rows use this as the group_id FK.
# 2026-04-17 chat-16: replaces the chat-15 content-hash scheme (see
# registries/DATABASE_DESIGN.md). Must be a group registered in
# sample_registry (default "all_226" — registered by load_bridge.R STEP 6.5).
export SAMPLE_GROUP="${SAMPLE_GROUP:-all_226}"

# ── Window registry ──────────────────────────────────────────────────────────
export WINDOW_REGISTRY="${BASE}/inversion_localpca_v7/06_mds_candidates/snake_regions_multiscale/precomp/windows_master.tsv.gz"

# ── Genome ───────────────────────────────────────────────────────────────────
export GENOME_SIZE=963905721
export N_SAMPLES=226
export N_CHROMOSOMES=28

# ── Engine A: NGSadmix parameters ────────────────────────────────────────────
export K_MIN=2
export K_MAX=20
export SEEDS=(1 2 3 4 5)
export NGSADMIX_MINMAF=0.05
export NGSADMIX_MAXITER=250

# ── Engine B: Instant Q parameters ──────────────────────────────────────────
export IQ_WINDOW_SIZE=100
export IQ_WINDOW_STEP=20
export IQ_EM_ITER=20
export IQ_TOL="1e-4"
export IQ_NCORES=4

# ── ANGSD parameters ─────────────────────────────────────────────────────────
export ANGSD_GL=1
export ANGSD_MINQ=25
export ANGSD_MINMAPQ=25
export ANGSD_BAQ=1
export ANGSD_C=50
export ANGSD_SETMINDEPTHIND=3
export ANGSD_SETMAXDEPTHIND=57
export ANGSD_MININD=200
export ANGSD_SNP_PVAL="1e-6"
export ANGSD_MINMAF=0.05

# ── Module 2C: SNP Q Support ────────────────────────────────────────────────
export MODULE_2C_DIR="${BASE}/module02c_snp_q_support"
export MODULE_2C_WINDOW_COUNT=50
export MODULE_2C_WINDOW_STEP=25

# ── Stats dispatcher + C popstats ────────────────────────────────────────────
export STATS_CACHE_DIR="${BASE}/unified_ancestry/region_stats_cache"

# ── Region popstats defaults ────────────────────────────────────────────────
export POPSTATS_FIXED_WIN="50000:10000"
export POPSTATS_DOWNSAMPLE=1
export POPSTATS_TYPE=2

# ── Hobs/HWE engine ─────────────────────────────────────────────────────────
export HOBS_OUTDIR="${BASE}/hobs_hwe_confirmation"

# ── Inversion pipeline integration ──────────────────────────────────────────
export INVDIR="${BASE}/inversion_localpca_v7"
export SNAKE1_DIR="${INVDIR}/06_mds_candidates/snake_regions_multiscale"
export SNAKE_CAND_FILE="${SNAKE1_DIR}/snake_candidate_regions.tsv.gz"

# ── SLURM defaults ───────────────────────────────────────────────────────────
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-lt200308}"
export SLURM_PARTITION="${SLURM_PARTITION:-compute}"
export SLURM_DEFAULT_TIME="02:00:00"
export SLURM_DEFAULT_MEM="32G"
export SLURM_DEFAULT_CPUS=8

# ── Palette ──────────────────────────────────────────────────────────────────
export ANCESTRY_PALETTE=(
  "#4E79A7" "#F28E2B" "#E15759" "#76B7B2" "#59A14F"
  "#EDC948" "#B07AA1" "#FF9DA7" "#9C755F" "#BAB0AC"
  "#86BCB6" "#D37295"
)

# ── Convenience ──────────────────────────────────────────────────────────────
timestamp(){ date '+%F %T'; }
export -f timestamp

anc_log()  { echo "[$(timestamp)] [ancestry] $*"; }
anc_die()  { echo "[$(timestamp)] [ancestry] FATAL: $*" >&2; exit 1; }
anc_check_file() { [[ -s "$1" ]] || anc_die "Missing or empty: $1"; }

anc_init_dirs() {
  mkdir -p "${LOCAL_Q_DIR}" "${STATS_CACHE_DIR}" "${HOBS_OUTDIR}" \
           "${REGISTRY_DIR}/groups" "${REGISTRY_DIR}/backups"
}
export -f anc_log anc_die anc_check_file anc_init_dirs
