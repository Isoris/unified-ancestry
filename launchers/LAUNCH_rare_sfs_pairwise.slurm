#!/usr/bin/env bash
# =============================================================================
# LAUNCH_rare_sfs_pairwise.slurm — Pairwise rare allele sharing from BEAGLE
#
# Per-chromosome array (1-28). Computes doubleton/tripleton/quad sharing
# matrices per chromosome, then merges genome-wide in a final step.
#
# Usage:
#   sbatch launchers/LAUNCH_rare_sfs_pairwise.slurm
#   BINS="2,3,4,5,10" sbatch ...
# =============================================================================
#SBATCH --job-name=rare_sfs
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-29
#SBATCH --output=logs/rare_sfs_%a.out
#SBATCH --error=logs/rare_sfs_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${CONFIG:-${SCRIPT_DIR}/00_ancestry_config.sh}"
[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"

eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-assembly}"
mkdir -p logs

BINS="${BINS:-2,3,4,5}"
OUTDIR="${BASE}/unified_ancestry/rare_sfs_cache"
mkdir -p "$OUTDIR"

# Build binary
SFS_BIN="${SCRIPT_DIR}/engines/fst_dxy/rare_sfs_pairwise"
if [[ ! -x "$SFS_BIN" ]]; then
  anc_log "Compiling rare_sfs_pairwise..."
  cd "${SCRIPT_DIR}/engines/fst_dxy" && make rare_sfs_pairwise && cd -
fi

# Build group file from ancestry clusters
GROUPS_FILE="${OUTDIR}/sample_groups.tsv"
if [[ ! -f "$GROUPS_FILE" ]]; then
  SUBSETS_DIR="${HOBS_OUTDIR:-${BASE}/hobs_hwe_confirmation}/subsets"
  if [[ -d "$SUBSETS_DIR" ]]; then
    > "$GROUPS_FILE"
    for k in $(seq 1 "${DEFAULT_K}"); do
      SFILE="${SUBSETS_DIR}/Q${k}_full.samples"
      [[ -f "$SFILE" ]] || continue
      while read sid; do
        echo -e "${sid}\tQ${k}" >> "$GROUPS_FILE"
      done < "$SFILE"
    done
    anc_log "Built group file: $GROUPS_FILE"
  fi
fi
GROUPS_ARG=""
[[ -f "$GROUPS_FILE" ]] && GROUPS_ARG="--groups $GROUPS_FILE"

# ── Map array index ──
mapfile -t CHROMS < <(awk '{print $1}' "$REF_FAI" | sort -V)
N_CHROMS=${#CHROMS[@]}
IDX=$((SLURM_ARRAY_TASK_ID - 1))

# Task 29 = genome-wide merge
if [[ $IDX -ge $N_CHROMS ]]; then
  anc_log "Merging genome-wide rare SFS matrices..."

  IFS=',' read -ra BINS_ARR <<< "$BINS"
  for b in "${BINS_ARR[@]}"; do
    # Sum all per-chr matrices
    FIRST=1
    for CHR in "${CHROMS[@]}"; do
      CHR_FILE="${OUTDIR}/${CHR}.bin${b}.tsv"
      [[ -f "$CHR_FILE" ]] || continue

      if [[ $FIRST -eq 1 ]]; then
        cp "$CHR_FILE" "${OUTDIR}/genome.bin${b}.tsv"
        FIRST=0
      else
        # Add matrices (skip header, sum values)
        python3 -c "
import sys
import csv
f1 = '${OUTDIR}/genome.bin${b}.tsv'
f2 = '${CHR_FILE}'
with open(f1) as a, open(f2) as b_f:
    ra = list(csv.reader(a, delimiter='\t'))
    rb = list(csv.reader(b_f, delimiter='\t'))
    hdr = ra[0]
    with open(f1 + '.tmp', 'w') as out:
        out.write('\t'.join(hdr) + '\n')
        for i in range(1, len(ra)):
            row = [ra[i][0]]
            for j in range(1, len(ra[i])):
                v = int(ra[i][j]) + int(rb[i][j])
                row.append(str(v))
            out.write('\t'.join(row) + '\n')
import shutil; shutil.move(f1 + '.tmp', f1)
"
      fi
    done
    anc_log "Merged bin${b}: ${OUTDIR}/genome.bin${b}.tsv"
  done

  # Also merge group summaries
  for b in "${BINS_ARR[@]}"; do
    FIRST=1
    for CHR in "${CHROMS[@]}"; do
      GF="${OUTDIR}/${CHR}.group_summary.bin${b}.tsv"
      [[ -f "$GF" ]] || continue
      if [[ $FIRST -eq 1 ]]; then
        cp "$GF" "${OUTDIR}/genome.group_summary.bin${b}.tsv"
        FIRST=0
      else
        python3 -c "
import csv
f1 = '${OUTDIR}/genome.group_summary.bin${b}.tsv'
f2 = '${GF}'
with open(f1) as a, open(f2) as b_f:
    ra = list(csv.reader(a, delimiter='\t'))
    rb = list(csv.reader(b_f, delimiter='\t'))
    hdr = ra[0]
    with open(f1 + '.tmp', 'w') as out:
        out.write('\t'.join(hdr) + '\n')
        for i in range(1, len(ra)):
            row = [ra[i][0]]
            for j in range(1, len(ra[i])):
                v = float(ra[i][j]) + float(rb[i][j])
                row.append(f'{v:.2f}')
            out.write('\t'.join(row) + '\n')
import shutil; shutil.move(f1 + '.tmp', f1)
"
      fi
    done
  done

  anc_log "Genome-wide merge complete"
  exit 0
fi

# ── Per-chromosome computation ──
CHR="${CHROMS[$IDX]}"
BGL=$(ls "${BEAGLE_DIR}"/*${CHR}*.beagle.gz 2>/dev/null | head -1)
[[ -n "$BGL" && -f "$BGL" ]] || { anc_log "No BEAGLE for ${CHR}"; exit 0; }

OUT_FILE="${OUTDIR}/${CHR}"
if [[ -f "${OUT_FILE}.bin2.tsv" && "${FORCE:-0}" -ne 1 ]]; then
  anc_log "${CHR} exists — skip"
  exit 0
fi

anc_log "rare_sfs_pairwise: ${CHR}, bins=${BINS}"
"$SFS_BIN" \
  --beagle "$BGL" \
  --sample_list "$SAMPLE_LIST" \
  --chr "$CHR" \
  --bins "$BINS" \
  --out_prefix "$OUT_FILE" \
  $GROUPS_ARG

anc_log "Done: ${CHR}"
