#!/usr/bin/env bash
# =============================================================================
# STEP_HH_B_run_angsd_hwe.sh — Run ANGSD -doHWE per subset per chromosome
#
# Uses the patched ANGSD binary (Isoris/angsd_fixed_HWE) with -doHetFreq 1
# to output hetFreq column without activating het-frequency filtering.
#
# Produces: <outdir>/<subset>/<chr>.hwe.gz
# Each file has columns: Chromo Position Major Minor hweFreq Freq F LRT p-value hetFreq
#
# Citation: ANGSD (Korneliussen et al. 2014), patched for -doHetFreq
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/00_hobs_hwe_config.sh"

# ── Parse args ──
SUBSET_ID="${1:-}"
CHR="${2:-}"

if [[ -z "$SUBSET_ID" || -z "$CHR" ]]; then
  echo "Usage: $0 <subset_id> <chromosome>"
  echo "  subset_id: Q1_full, Q1_pruned, allSamples_referenceOnly, etc."
  echo "  chromosome: C_gar_LG01, C_gar_LG02, ..."
  exit 1
fi

SUBSETS_DIR="${HOBS_OUTDIR}/subsets"
MANIFEST="${SUBSETS_DIR}/subset_manifest.tsv"
[[ -f "$MANIFEST" ]] || hobs_die "Run 01_build_subset_bamlists.sh first"

# ── Look up BAM list for this subset ──
BAMLIST_SUB=$(awk -F'\t' -v sid="$SUBSET_ID" '$1==sid {print $4}' "$MANIFEST")
[[ -n "$BAMLIST_SUB" && -f "$BAMLIST_SUB" ]] || hobs_die "Subset $SUBSET_ID not found in manifest"
N_SAMPLES=$(wc -l < "$BAMLIST_SUB")

# ── Compute minInd for this subset ──
MININD=$(python3 -c "import math; print(max(2, math.floor($N_SAMPLES * $HWE_MININD_FRAC)))")

# ── Output ──
OUTDIR="${HOBS_OUTDIR}/02_angsd_hwe/${SUBSET_ID}"
mkdir -p "$OUTDIR"
OUT_PREFIX="${OUTDIR}/${CHR}"

# ── Skip if exists ──
if [[ -f "${OUT_PREFIX}.hwe.gz" && "${FORCE:-0}" -ne 1 ]]; then
  hobs_log "Exists: ${OUT_PREFIX}.hwe.gz — skip"
  exit 0
fi

hobs_log "ANGSD -doHWE: subset=${SUBSET_ID} chr=${CHR} N=${N_SAMPLES} minInd=${MININD}"

# ── Select ANGSD binary ──
ANGSD="$ANGSD_BIN"
if [[ -x "$ANGSD_PATCHED" ]]; then
  ANGSD="$ANGSD_PATCHED"
  hobs_log "  Using patched ANGSD: $ANGSD"
fi

# ── Run ANGSD ──
$ANGSD \
  -bam "$BAMLIST_SUB" \
  -ref "$REF" \
  -r "$CHR": \
  -GL "$HWE_GL" \
  -doMajorMinor "$HWE_DOMAJMIN" \
  -doMaf "$HWE_DOMAF" \
  -doHWE "$HWE_DOHWE" \
  -doHetFreq "$HWE_DOHETFREQ" \
  -minQ "$HWE_MINQ" \
  -minMapQ "$HWE_MINMAPQ" \
  -baq "$HWE_BAQ" \
  -C "$HWE_C" \
  -setMinDepthInd "$HWE_SETMINDEPTHIND" \
  -setMaxDepthInd "$HWE_SETMAXDEPTHIND" \
  -minInd "$MININD" \
  -SNP_pval "$HWE_SNP_PVAL" \
  -minMaf "$HWE_MINMAF" \
  -nThreads "${SLURM_CPUS_PER_TASK:-$HOBS_SLURM_CPUS}" \
  -out "$OUT_PREFIX"

# ── Verify output ──
if [[ -f "${OUT_PREFIX}.hwe.gz" ]]; then
  N_SITES=$(zcat "${OUT_PREFIX}.hwe.gz" | tail -n +2 | wc -l)
  hobs_log "  Output: ${OUT_PREFIX}.hwe.gz ($N_SITES sites)"
else
  hobs_log "  WARNING: No output produced for ${SUBSET_ID}/${CHR}"
fi
