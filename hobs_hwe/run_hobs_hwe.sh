#!/usr/bin/env bash
# =============================================================================
# run_hobs_hwe.sh — Hobs/HWE Secondary Confirmation Module
#
# Subset-aware, multiscale Hobs/HWE confirmation for inversion candidates.
# NOT a primary detection layer — secondary confirmation only.
#
# Steps:
#   1. Build subset BAM lists from ancestry clusters
#   2. Run ANGSD -doHWE per subset × chromosome
#   3. Compute site-level Hobs + multi-scale windowed summaries
#   4. Candidate interval overlay
#
# Usage:
#   bash run_hobs_hwe.sh --step 1          # just build subsets
#   bash run_hobs_hwe.sh --from 2 --to 3   # ANGSD + windowing
#   bash run_hobs_hwe.sh --all             # everything (serial, slow)
#   sbatch ../launchers/LAUNCH_hobs_hwe.slurm   # parallel (recommended)
#
# Citation:
#   Hobs sliding-window concept: Claire Mérot (angsd_pipeline)
#   ANGSD: Korneliussen et al. (2014)
#   Patched ANGSD: github.com/Isoris/angsd_fixed_HWE
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_hobs_hwe_config.sh"

START_STEP=1; END_STEP=99
while [[ $# -gt 0 ]]; do
  case "$1" in
    --step) START_STEP="$2"; END_STEP="$2"; shift 2;;
    --from) START_STEP="$2"; shift 2;;
    --to)   END_STEP="$2"; shift 2;;
    --all)  START_STEP=1; END_STEP=4; shift;;
    --force) export FORCE=1; shift;;
    *) shift;;
  esac
done

hobs_log "=== Hobs/HWE Confirmation Module ==="
hobs_log "Steps: $START_STEP to $END_STEP"

# Step 1
if [[ $START_STEP -le 1 && $END_STEP -ge 1 ]]; then
  hobs_log "Step 1: Build subset BAM lists..."
  bash "${SCRIPT_DIR}/STEP_HH_A_build_subset_bamlists.sh"
fi

# Step 2+3 (serial — use SLURM for parallel)
if [[ $START_STEP -le 3 && $END_STEP -ge 2 ]]; then
  hobs_log "Step 2+3: ANGSD + windowing (serial)..."
  MANIFEST="${HOBS_OUTDIR}/subsets/subset_manifest.tsv"
  mapfile -t CHROMS < <(awk '{print $1}' "$REF_FAI")

  while IFS=$'\t' read -r sid stype nsamp bam samp notes; do
    [[ "$sid" == "subset_id" ]] && continue
    for chr in "${CHROMS[@]}"; do
      if [[ $START_STEP -le 2 && $END_STEP -ge 2 ]]; then
        bash "${SCRIPT_DIR}/STEP_HH_B_run_angsd_hwe.sh" "$sid" "$chr" || true
      fi
      if [[ $START_STEP -le 3 && $END_STEP -ge 3 ]]; then
        bash "${SCRIPT_DIR}/STEP_HH_C_compute_hobs_windows.sh" "$sid" "$chr" || true
      fi
    done
  done < "$MANIFEST"
fi

# Step 4
if [[ $START_STEP -le 4 && $END_STEP -ge 4 ]]; then
  hobs_log "Step 4: Candidate overlay..."
  CAND="${SNAKE_CAND_FILE:-}"
  if [[ -n "$CAND" && -f "$CAND" ]]; then
    OVERLAY_DIR="${HOBS_OUTDIR}/04_candidate_overlay"
    mkdir -p "$OVERLAY_DIR"
    for scale_def in "${WINDOW_SCALES[@]}"; do
      label=$(echo "$scale_def" | cut -d: -f1)
      python3 "$HOBS_CANDIDATE_OVERLAY_PY" \
        --candidates "$CAND" \
        --windows_dir "${HOBS_OUTDIR}/03_hobs_windows" \
        --subsets "${HOBS_OUTDIR}/subsets/subset_manifest.tsv" \
        --scale "$label" \
        --outdir "$OVERLAY_DIR"
    done
  else
    hobs_log "  No candidate file — skip overlay"
  fi
fi

hobs_log "=== Hobs/HWE module complete ==="
