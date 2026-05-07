#!/usr/bin/env bash
# =============================================================================
# 03_compute_hobs_windows.sh — Site-level Hobs + multi-scale windows
#
# For each subset × chromosome, reads ANGSD .hwe.gz and runs hobs_windower
# to produce site-level Hobs/F and 7-scale windowed summaries.
#
# Usage:
#   bash 03_compute_hobs_windows.sh <subset_id> <chromosome>
#   bash 03_compute_hobs_windows.sh Q1_full C_gar_LG01
#
# Or all subsets × all chromosomes:
#   bash 03_compute_hobs_windows.sh --all
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../00_hobs_hwe_config.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Compile C windower if needed ──
WINDOWER="${SCRIPT_DIR}/hobs_windower"
if [[ ! -x "$WINDOWER" ]]; then
  hobs_log "Compiling hobs_windower..."
  gcc -O3 -o "$WINDOWER" "${SCRIPT_DIR}/hobs_windower.c" -lz -lm
  hobs_log "  Done: $WINDOWER"
fi

# ── Build scales string for C binary ──
SCALES_STR=""
for s in "${WINDOW_SCALES[@]}"; do
  [[ -n "$SCALES_STR" ]] && SCALES_STR="${SCALES_STR},"
  SCALES_STR="${SCALES_STR}${s}"
done

# ── Parse args ──
if [[ "${1:-}" == "--all" ]]; then
  # Run all subsets × all chromosomes
  MANIFEST="${HOBS_OUTDIR}/subsets/subset_manifest.tsv"
  [[ -f "$MANIFEST" ]] || hobs_die "Run 01 first"

  mapfile -t CHROMS < <(awk '{print $1}' "$REF_FAI")

  while IFS=$'\t' read -r sid stype nsamp bam samp notes; do
    [[ "$sid" == "subset_id" ]] && continue
    for chr in "${CHROMS[@]}"; do
      HWE_FILE="${HOBS_OUTDIR}/02_angsd_hwe/${sid}/${chr}.hwe.gz"
      [[ -f "$HWE_FILE" ]] || continue

      OUTDIR="${HOBS_OUTDIR}/03_hobs_windows/${sid}"
      mkdir -p "$OUTDIR"
      PREFIX="${OUTDIR}/${chr}"

      # Get chromosome size
      CHR_SIZE=$(awk -v c="$chr" '$1==c {print $2}' "$REF_FAI")
      [[ -n "$CHR_SIZE" ]] || continue

      # Skip if already done
      if [[ -f "${PREFIX}.win1Mb.tsv" && "${FORCE:-0}" -ne 1 ]]; then
        continue
      fi

      hobs_log "Windowing: ${sid} / ${chr}"
      "$WINDOWER" "$HWE_FILE" "$PREFIX" "$CHR_SIZE" \
        --scales "$SCALES_STR" --mad_n "$OUTLIER_N_MAD"
    done
  done < "$MANIFEST"

  hobs_log "All windowing complete"
  exit 0
fi

# ── Single subset × chromosome ──
SUBSET_ID="${1:-}"
CHR="${2:-}"
[[ -n "$SUBSET_ID" && -n "$CHR" ]] || { echo "Usage: $0 <subset_id> <chr>"; exit 1; }

HWE_FILE="${HOBS_OUTDIR}/02_angsd_hwe/${SUBSET_ID}/${CHR}.hwe.gz"
[[ -f "$HWE_FILE" ]] || hobs_die "Missing: $HWE_FILE — run 02 first"

OUTDIR="${HOBS_OUTDIR}/03_hobs_windows/${SUBSET_ID}"
mkdir -p "$OUTDIR"
PREFIX="${OUTDIR}/${CHR}"

CHR_SIZE=$(awk -v c="$CHR" '$1==c {print $2}' "$REF_FAI")
[[ -n "$CHR_SIZE" ]] || hobs_die "Chromosome $CHR not in FAI"

hobs_log "Windowing: ${SUBSET_ID} / ${CHR} (${CHR_SIZE} bp)"
"$WINDOWER" "$HWE_FILE" "$PREFIX" "$CHR_SIZE" \
  --scales "$SCALES_STR" --mad_n "$OUTLIER_N_MAD"
