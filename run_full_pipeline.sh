#!/usr/bin/env bash
# =============================================================================
# run_full_pipeline.sh — Orchestrate the complete unified ancestry pipeline
#
# v12.2 REWIRED (repo reorg):
#   - Uses moved-script variables (SNP_Q_SUPPORT_PY, NESTED_COMPOSITION_PY,
#     CANDIDATE_CLASSIFIER_PY, EXPORT_MODULE5B_PY, HOBS_PLOT_R) from config
#   - Moved scripts now live under Modules/MODULE_5A/ and Modules/MODULE_5B/
#
# v12.1 REWIRED:
#   - Sources pipeline_bridge.sh at top
#   - Step 1 (registries) uses sample_registry if available
#   - Step 8 (popstats) can pull groups from registry
#   - Step 9 (hobs) bamlist builder accepts --from-registry
#
# Usage:
#   bash run_full_pipeline.sh [--config 00_ancestry_config.sh] [--step N]
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/00_ancestry_config.sh"
START_STEP=1
END_STEP=99

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2;;
    --step)   START_STEP="$2"; END_STEP="$2"; shift 2;;
    --from)   START_STEP="$2"; shift 2;;
    --to)     END_STEP="$2"; shift 2;;
    *) shift;;
  esac
done

[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG"; exit 1; }
source "$CONFIG"

# v12.1: Source pipeline bridge
BRIDGE_SH="${BASE}/inversion_codebase_v8.5/utils/pipeline_bridge.sh"
if [[ -f "$BRIDGE_SH" ]]; then
  source "$BRIDGE_SH"
fi

RSCRIPT="${RSCRIPT_BIN:-Rscript}"
PYTHON="${PYTHON_BIN:-python3}"

anc_log "=== Unified Ancestry Pipeline (v12.2 REWIRED) ==="
anc_log "Config: $CONFIG"
anc_log "Steps: ${START_STEP} to ${END_STEP}"
anc_log "Registry: ${REGISTRY_DIR:-NOT SET}"

# ── Step 1: Build registries ──
if [[ $START_STEP -le 1 && $END_STEP -ge 1 ]]; then
  anc_log "Step 1: Building registries..."
  $PYTHON "${SCRIPT_DIR}/registries/build_registries.py" \
    --config "$CONFIG" \
    --outdir "${BASE}/unified_ancestry/registries/" \
    ${SNAKE_CAND_FILE:+--candidates "$SNAKE_CAND_FILE"}
fi

# ── Step 2: Engine B precompute ──
if [[ $START_STEP -le 2 && $END_STEP -ge 2 ]]; then
  anc_log "Step 2: Engine B precompute (instant Q)..."
  "${SCRIPT_DIR}/instant_q" --config "$CONFIG" --precompute --all --ncores "${IQ_NCORES:-4}"
fi

# ── Step 3: SNP Q support ──
if [[ $START_STEP -le 3 && $END_STEP -ge 3 ]]; then
  anc_log "Step 3: SNP Q support..."
  SNP_OUTDIR="${BASE}/unified_ancestry/snp_q_support"
  mkdir -p "$SNP_OUTDIR"

  for BGL in "${BEAGLE_DIR}"/*.beagle.gz; do
    [[ -f "$BGL" ]] || continue
    CHR=$(basename "$BGL" | grep -oP 'C_gar_LG\d+' | head -1)
    [[ -n "$CHR" ]] || continue

    if [[ -f "${SNP_OUTDIR}/q_support_per_snp_${CHR}.tsv" && "${FORCE:-0}" -ne 1 ]]; then
      continue
    fi

    anc_log "  SNP support: ${CHR}..."
    $PYTHON "$SNP_Q_SUPPORT_PY" \
      --beagle "$BGL" \
      --qopt "$BEST_QOPT" \
      --sample_list "$SAMPLE_LIST" \
      --outdir "${SNP_OUTDIR}/${CHR}" \
      --K "$DEFAULT_K" \
      --chr "$CHR" \
      ${STEP1_PRUNED_LIST:+--keep "$STEP1_PRUNED_LIST"}
  done
fi

# ── Step 4: Nested composition ──
if [[ $START_STEP -le 4 && $END_STEP -ge 4 ]]; then
  anc_log "Step 4: Nested composition..."
  NEST_OUTDIR="${BASE}/unified_ancestry/nested_composition"
  mkdir -p "$NEST_OUTDIR"

  PARENTS="${SNAKE_CAND_FILE:-}"
  if [[ -n "$PARENTS" && -f "$PARENTS" ]]; then
    $PYTHON "$NESTED_COMPOSITION_PY" \
      --q_cache_dir "$LOCAL_Q_DIR" \
      --parents "$PARENTS" \
      --outdir "$NEST_OUTDIR" \
      --K "$DEFAULT_K" \
      --sample_list "$SAMPLE_LIST"
  else
    anc_log "  No candidate file — skip"
  fi
fi

# ── Step 5: Candidate classification ──
if [[ $START_STEP -le 5 && $END_STEP -ge 5 ]]; then
  anc_log "Step 5: Candidate classification..."
  CLASS_OUTDIR="${BASE}/unified_ancestry/candidate_classification"
  mkdir -p "$CLASS_OUTDIR"

  PARENTS="${SNAKE_CAND_FILE:-}"
  if [[ -n "$PARENTS" && -f "$PARENTS" ]]; then
    Q_WIN=""
    for d in "${BASE}/unified_ancestry/snp_q_support"/*/; do
      qw="${d}q_support_windows.tsv"
      if [[ -f "$qw" ]]; then
        if [[ -z "$Q_WIN" ]]; then
          cp "$qw" "${CLASS_OUTDIR}/q_support_windows_merged.tsv"
          Q_WIN="${CLASS_OUTDIR}/q_support_windows_merged.tsv"
        else
          tail -n +2 "$qw" >> "$Q_WIN"
        fi
      fi
    done

    NEST_SUMM="${BASE}/unified_ancestry/nested_composition/nested_composition_summary.tsv"

    $PYTHON "$CANDIDATE_CLASSIFIER_PY" \
      --candidates "$PARENTS" \
      ${Q_WIN:+--q_windows "$Q_WIN"} \
      --q_cache_dir "$LOCAL_Q_DIR" \
      ${NEST_SUMM:+--nested "$NEST_SUMM"} \
      --outdir "$CLASS_OUTDIR" \
      --K "$DEFAULT_K"
  fi
fi

# ── Step 6: Module 5B exports ──
if [[ $START_STEP -le 6 && $END_STEP -ge 6 ]]; then
  anc_log "Step 6: Module 5B exports..."
  EXPORT_OUTDIR="${BASE}/unified_ancestry/exports_for_module5b"
  mkdir -p "$EXPORT_OUTDIR"

  SNP_FILE=""
  for d in "${BASE}/unified_ancestry/snp_q_support"/*/; do
    sf="${d}q_support_per_snp.tsv"
    if [[ -f "$sf" ]]; then
      if [[ -z "$SNP_FILE" ]]; then
        cp "$sf" "${EXPORT_OUTDIR}/q_support_per_snp_merged.tsv"
        SNP_FILE="${EXPORT_OUTDIR}/q_support_per_snp_merged.tsv"
      else
        tail -n +2 "$sf" >> "$SNP_FILE"
      fi
    fi
  done

  if [[ -n "$SNP_FILE" ]]; then
    $PYTHON "$EXPORT_MODULE5B_PY" \
      --snp_support "$SNP_FILE" \
      --outdir "$EXPORT_OUTDIR" \
      --K "$DEFAULT_K"
  fi
fi

# ── Step 7: Region stats precompute (Q overlay) ──
if [[ $START_STEP -le 7 && $END_STEP -ge 7 ]]; then
  anc_log "Step 7: Region stats precompute (Q overlay)..."
  $RSCRIPT -e '
    source("'"${LOAD_BRIDGE}"'")
    chroms <- readLines("'"$REF_FAI"'")
    chroms <- sub("\\t.*", "", chroms)
    for (chr in chroms) {
      message("[pipeline] Q overlay for ", chr)
      tryCatch({
        q <- get_Q_summary(chr)
        message("[pipeline]   ", nrow(q), " windows")
      }, error = function(e) {
        message("[pipeline] WARN: ", chr, " failed: ", e$message)
      })
    }
  '
fi

# ── Step 8: Region popstats (C binary) ──
if [[ $START_STEP -le 8 && $END_STEP -ge 8 ]]; then
  anc_log "Step 8: Region popstats (C binary)..."

  if [[ ! -x "$POPSTATS_BIN" ]]; then
    anc_log "  Compiling region_popstats..."
    make -C "${SCRIPT_DIR}/engines/fst_dxy" -j4 2>/dev/null || true
  fi

  POPSTATS_OUTDIR="${STATS_CACHE_DIR:-${BASE}/unified_ancestry/region_stats_cache}"
  mkdir -p "$POPSTATS_OUTDIR"

  # v12.1: Build group spec from registry if available
  GROUPS_SPEC=""
  if [[ -f "${LOAD_BRIDGE}" ]]; then
    # Try registry-based group resolution
    GROUPS_SPEC=$("${RSCRIPT_BIN}" --vanilla -e "
      source('${LOAD_BRIDGE}')
      K <- ${DEFAULT_K}
      specs <- character(0)
      for (qi in 1:K) {
        gid <- paste0('ancestry_K', K, '_Q', qi)
        if (reg\$has_group(gid)) {
          members <- reg\$get_group(gid)
          if (length(members) > 0) {
            tf <- tempfile(pattern = paste0('grp_Q', qi, '_'), fileext = '.txt')
            writeLines(members, tf)
            specs <- c(specs, paste0('Q', qi, ':', tf))
          }
        }
      }
      if (length(specs) > 0) cat(paste(specs, collapse = ','))
    " 2>/dev/null || echo "")
  fi

  # Fallback to hobs subsets
  if [[ -z "$GROUPS_SPEC" ]]; then
    SUBSETS_DIR="${HOBS_OUTDIR:-${BASE}/hobs_hwe_confirmation}/subsets"
    if [[ -d "$SUBSETS_DIR" ]]; then
      for k in $(seq 1 "${DEFAULT_K}"); do
        SFILE="${SUBSETS_DIR}/Q${k}_full.samples"
        [[ -f "$SFILE" ]] || continue
        [[ -n "$GROUPS_SPEC" ]] && GROUPS_SPEC="${GROUPS_SPEC},"
        GROUPS_SPEC="${GROUPS_SPEC}Q${k}:${SFILE}"
      done
    fi
  fi

  mapfile -t CHROMS < <(awk '{print $1}' "$REF_FAI")

  for CHR in "${CHROMS[@]}"; do
    OUT_FILE="${POPSTATS_OUTDIR}/${CHR}.popstats.tsv"
    if [[ -f "$OUT_FILE" && "${FORCE:-0}" -ne 1 ]]; then
      continue
    fi

    BGL=$(ls "${BEAGLE_DIR}"/*${CHR}*.beagle.gz 2>/dev/null | head -1)
    [[ -n "$BGL" && -f "$BGL" ]] || continue

    GROUPS_ARG=""
    [[ -n "$GROUPS_SPEC" ]] && GROUPS_ARG="--groups ${GROUPS_SPEC}"

    "$POPSTATS_BIN" \
      --beagle "$BGL" \
      --sample_list "$SAMPLE_LIST" \
      --chr "$CHR" \
      --fixed_win "50000:10000" \
      --downsample 1 \
      --type 2 \
      --ncores "${SLURM_CPUS_PER_TASK:-4}" \
      --out "$OUT_FILE" \
      $GROUPS_ARG || anc_log "  WARN: ${CHR} popstats failed"
  done
fi

# ── Step 9: Hobs/HWE plotting ──
if [[ $START_STEP -le 9 && $END_STEP -ge 9 ]]; then
  anc_log "Step 9: Hobs/HWE plotting suite..."
  PLOT_SCRIPT="$HOBS_PLOT_R"
  HOBS_CONFIG="${SCRIPT_DIR}/engines/hobs_hwe/00_hobs_hwe_config.sh"

  if [[ -f "$PLOT_SCRIPT" && -f "$HOBS_CONFIG" ]]; then
    $RSCRIPT "$PLOT_SCRIPT" \
      --config "$HOBS_CONFIG" \
      --mode all \
      --scale 50kb \
      ${SNAKE_CAND_FILE:+--candidates "$SNAKE_CAND_FILE"}
  fi
fi

anc_log "=== Pipeline complete (v12.2) ==="
