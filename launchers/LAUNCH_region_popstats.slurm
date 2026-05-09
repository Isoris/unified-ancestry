#!/usr/bin/env bash
# =============================================================================
# LAUNCH_region_popstats.slurm — SLURM array for C popstats (v12.1)
#
# v12.1: Sources pipeline_bridge.sh. GROUPS_FROM=registry pulls from
#         sample_registry instead of flat files.
#
# Usage:
#   sbatch launchers/LAUNCH_region_popstats.slurm
#   GROUPS_FROM=registry sbatch launchers/LAUNCH_region_popstats.slurm
# =============================================================================
#SBATCH --job-name=popstats
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-28
#SBATCH --output=logs/popstats_%a.out
#SBATCH --error=logs/popstats_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${CONFIG:-${SCRIPT_DIR}/00_ancestry_config.sh}"
[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"

# v12.1: Source pipeline bridge
BRIDGE_SH="${BASE}/inversion_codebase_v8.5/utils/pipeline_bridge.sh"
[[ -f "$BRIDGE_SH" ]] && source "$BRIDGE_SH"

eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-assembly}"
mkdir -p logs

GROUPS_FROM="${GROUPS_FROM:-ancestry}"
FIXED_WIN="${FIXED_WIN:-50000:10000}"
DOWNSAMPLE="${DOWNSAMPLE:-1}"
TYPE="${TYPE:-2}"
FORCE="${FORCE:-0}"

# ── Build binary if needed ──
if [[ ! -x "$POPSTATS_BIN" ]]; then
  anc_log "Compiling region_popstats..."
  make -C "${SCRIPT_DIR}/engines" region_popstats -j "${SLURM_CPUS_PER_TASK:-4}" 2>/dev/null || true
fi

# ── Map array index to chromosome ──
mapfile -t CHROMS < <(awk '{print $1}' "$REF_FAI" | sort -V)
N_CHROMS=${#CHROMS[@]}
IDX=$((SLURM_ARRAY_TASK_ID - 1))
[[ $IDX -ge $N_CHROMS ]] && exit 0

CHR="${CHROMS[$IDX]}"
anc_log "Task ${SLURM_ARRAY_TASK_ID}/${N_CHROMS}: chr=${CHR}"

BGL=$(ls "${BEAGLE_DIR}"/*${CHR}*.beagle.gz 2>/dev/null | head -1)
[[ -n "$BGL" && -f "$BGL" ]] || { anc_log "No BEAGLE for ${CHR}"; exit 0; }

OUTDIR="${STATS_CACHE_DIR:-${BASE}/unified_ancestry/region_stats_cache}"
mkdir -p "$OUTDIR"
OUT_FILE="${OUTDIR}/${CHR}.popstats.tsv"
[[ -f "$OUT_FILE" && "$FORCE" -ne 1 ]] && { anc_log "Exists: skip"; exit 0; }

# ── Build group spec ──
GROUPS_ARG=""

if [[ "$GROUPS_FROM" == "registry" && -f "${LOAD_BRIDGE}" ]]; then
  # v12.1: Pull groups from sample_registry
  anc_log "Resolving groups from registry..."
  GROUPS_SPEC=$("${RSCRIPT_BIN}" --vanilla -e "
    source('${LOAD_BRIDGE}')
    K <- ${DEFAULT_K}
    specs <- character(0)
    for (qi in 1:K) {
      gid <- paste0('ancestry_K', K, '_Q', qi)
      if (reg\$has_group(gid)) {
        members <- reg\$get_group(gid)
        if (length(members) > 0) {
          tf <- file.path(tempdir(), paste0('grp_Q', qi, '.txt'))
          writeLines(members, tf)
          specs <- c(specs, paste0('Q', qi, ':', tf))
        }
      }
    }
    if (length(specs) > 0) cat(paste(specs, collapse = ','))
  " 2>/dev/null || echo "")

  if [[ -n "$GROUPS_SPEC" ]]; then
    GROUPS_ARG="--groups ${GROUPS_SPEC}"
    anc_log "Registry groups resolved: ${GROUPS_SPEC%%,*}..."
  else
    anc_log "No registry groups found, falling back to ancestry subsets"
  fi
fi

if [[ -z "$GROUPS_ARG" && "$GROUPS_FROM" == "ancestry" ]]; then
  SUBSETS_DIR="${HOBS_OUTDIR:-${BASE}/hobs_hwe_confirmation}/subsets"
  if [[ -d "$SUBSETS_DIR" ]]; then
    GROUPS_SPEC=""
    for k in $(seq 1 "${DEFAULT_K}"); do
      SFILE="${SUBSETS_DIR}/Q${k}_full.samples"
      [[ -f "$SFILE" ]] || continue
      [[ -n "$GROUPS_SPEC" ]] && GROUPS_SPEC="${GROUPS_SPEC},"
      GROUPS_SPEC="${GROUPS_SPEC}Q${k}:${SFILE}"
    done
    [[ -n "$GROUPS_SPEC" ]] && GROUPS_ARG="--groups ${GROUPS_SPEC}"
  fi
fi

# ── Run ──
anc_log "region_popstats: chr=${CHR} win=${FIXED_WIN}"
T_START=$(date +%s)

"$POPSTATS_BIN" \
  --beagle "$BGL" \
  --sample_list "$SAMPLE_LIST" \
  --chr "$CHR" \
  --fixed_win "$FIXED_WIN" \
  --downsample "$DOWNSAMPLE" \
  --type "$TYPE" \
  --ncores "${SLURM_CPUS_PER_TASK:-4}" \
  --out "$OUT_FILE" \
  $GROUPS_ARG

ELAPSED=$(( $(date +%s) - T_START ))
anc_log "Done: ${CHR} (${ELAPSED}s)"
