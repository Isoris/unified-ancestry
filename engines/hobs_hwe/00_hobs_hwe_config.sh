#!/usr/bin/env bash
# =============================================================================
# 00_hobs_hwe_config.sh — Hobs/HWE Secondary Confirmation Module
#
# Claire Mérot-style Hobs/HWE sliding-window concept adapted into a
# subset-aware, multiscale, robust summary framework for secondary
# inversion confirmation in a structured hatchery dataset.
#
# Citation note: Hobs/HWE sliding-window logic adapted from
#   Claire Mérot (github.com/clairemerot/angsd_pipeline)
#   ANGSD HWE from Korneliussen et al. (2014)
#   Patched ANGSD from github.com/Isoris/angsd_fixed_HWE
#
# Source this at the top of all module scripts.
# =============================================================================

# ── Inherit from unified ancestry config ──
ANCESTRY_CONFIG="${ANCESTRY_CONFIG:-$(dirname "${BASH_SOURCE[0]}")/../00_ancestry_config.sh}"
[[ -f "$ANCESTRY_CONFIG" ]] && source "$ANCESTRY_CONFIG"

# ── Module paths ──
export HOBS_MODULE_DIR="${BASE}/unified_ancestry/engines/hobs_hwe"
export HOBS_OUTDIR="${BASE}/hobs_hwe_confirmation"

# ── ANGSD binary (use patched version with -doHetFreq) ──
# This must be the Isoris/angsd_fixed_HWE patched build, NOT stock ANGSD.
# The patch adds -doHetFreq 1 (output-only, no filtering) and fixes
# abcFilterSNP null deref + abcFilter bounds check.
export ANGSD_BIN="${ANGSD_BIN:-angsd}"
export ANGSD_PATCHED="${ANGSD_PATCHED:-${BASE}/programs/angsd_fixed_HWE/angsd}"

# ── Input files ──
export BAMLIST="${BAMLIST:-${BASE}/het_roh/01_inputs_check/bamlist_qcpass.txt}"
export SAMPLE_LIST="${SAMPLE_LIST:-${STEP1_SAMPLE_LIST}}"

# ── Subset definitions ──
# Groups come from NGSadmix K=8 ancestry clusters.
# Each group gets its own Hobs/HWE scan.
export BEST_QOPT="${BEST_QOPT}"
export DEFAULT_K="${DEFAULT_K:-8}"
export PRUNED_LIST="${STEP1_PRUNED_LIST}"

# ── ANGSD parameters for -doHWE ──
export HWE_GL=1
export HWE_DOMAJMIN=1
export HWE_DOMAF=1
export HWE_DOHWE=1
export HWE_DOHETFREQ=1          # patched ANGSD: output hetFreq without filtering
export HWE_MINQ=25
export HWE_MINMAPQ=25
export HWE_BAQ=1
export HWE_C=50
export HWE_SETMINDEPTHIND=3
export HWE_SETMAXDEPTHIND=57
export HWE_MININD_FRAC=0.8      # minInd = floor(n_samples_in_subset * frac)
export HWE_SNP_PVAL="1e-6"
export HWE_MINMAF=0.05

# ── Window scales ──
# Format: "label:window_bp:step_bp"
export WINDOW_SCALES=(
  "5kb:5000:1000"
  "10kb:10000:2000"
  "50kb:50000:10000"
  "100kb:100000:20000"
  "250kb:250000:50000"
  "500kb:500000:100000"
  "1Mb:1000000:200000"
)

# ── Outlier thresholds ──
# Outliers defined relative to subset-level distributions (NOT pooled panel).
# Uses MAD-based robust thresholds: median ± N_MAD * MAD
export OUTLIER_N_MAD=3.0

# ── SLURM defaults ──
export HOBS_SLURM_ACCOUNT="${SLURM_ACCOUNT:-lt200308}"
export HOBS_SLURM_PARTITION="${SLURM_PARTITION:-compute}"
export HOBS_SLURM_TIME="04:00:00"
export HOBS_SLURM_MEM="16G"
export HOBS_SLURM_CPUS=8

# ── Convenience ──
hobs_log() { echo "[$(date '+%F %T')] [hobs] $*"; }
hobs_die() { echo "[$(date '+%F %T')] [hobs] FATAL: $*" >&2; exit 1; }
export -f hobs_log hobs_die
