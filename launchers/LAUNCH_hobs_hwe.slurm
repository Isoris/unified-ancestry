#!/usr/bin/env bash
# =============================================================================
# LAUNCH_hobs_hwe.slurm — SLURM array for ANGSD -doHWE + windowing
#
# Two-stage design:
#   Stage 1 (this script): ANGSD -doHWE per subset × chromosome
#   Stage 2 (--dependency): hobs_windower per subset × chromosome
#
# Array index maps to (subset_idx * n_chroms + chrom_idx).
# With 10 subsets × 28 chromosomes = 280 tasks.
#
# Usage:
#   # Build job list first
#   bash hobs_hwe/STEP_HH_A_build_subset_bamlists.sh
#   # Then submit
#   sbatch launchers/LAUNCH_hobs_hwe.slurm
# =============================================================================
#SBATCH --job-name=hobs_hwe
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-280
#SBATCH --output=logs/hobs_%a.out
#SBATCH --error=logs/hobs_%a.err

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HOBS_DIR="${SCRIPT_DIR}/hobs_hwe"
source "${HOBS_DIR}/00_hobs_hwe_config.sh"

eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-assembly}"

mkdir -p logs

# ── Build job list: subset_id × chromosome ──
MANIFEST="${HOBS_OUTDIR}/subsets/subset_manifest.tsv"
[[ -f "$MANIFEST" ]] || hobs_die "Run 01_build_subset_bamlists.sh first"

# Read subsets (skip header)
mapfile -t SUBSETS < <(tail -n +2 "$MANIFEST" | awk -F'\t' '{print $1}')
# Read chromosomes
mapfile -t CHROMS < <(awk '{print $1}' "$REF_FAI")

N_SUBSETS=${#SUBSETS[@]}
N_CHROMS=${#CHROMS[@]}
N_TOTAL=$((N_SUBSETS * N_CHROMS))

IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [[ $IDX -ge $N_TOTAL ]]; then
  hobs_log "Task $SLURM_ARRAY_TASK_ID exceeds $N_TOTAL jobs — skip"
  exit 0
fi

SUB_IDX=$((IDX / N_CHROMS))
CHR_IDX=$((IDX % N_CHROMS))
SUBSET="${SUBSETS[$SUB_IDX]}"
CHR="${CHROMS[$CHR_IDX]}"

hobs_log "Task $SLURM_ARRAY_TASK_ID: subset=$SUBSET chr=$CHR"

# ── Stage 1: ANGSD -doHWE ──
bash "${HOBS_DIR}/STEP_HH_B_run_angsd_hwe.sh" "$SUBSET" "$CHR"

# ── Stage 2: Windowing ──
bash "${HOBS_DIR}/STEP_HH_C_compute_hobs_windows.sh" "$SUBSET" "$CHR"

hobs_log "Done: $SUBSET / $CHR"
