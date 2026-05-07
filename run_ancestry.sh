#!/usr/bin/env bash
# =============================================================================
# run_ancestry.sh — Unified ancestry dispatcher (v12.1 REWIRED)
#
# v12.1: Sources pipeline_bridge.sh for cross-module env vars + registry
#
# Subcommands:
#   full_survey, instant_q, snp_support, region_stats, status, validate
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/00_ancestry_config.sh"
TASK=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2;;
    -h|--help)
      sed -n '2,/^# ====/p' "$0" | head -30
      exit 0;;
    *)
      if [[ -z "$TASK" ]]; then
        TASK="$1"; shift
      else
        EXTRA_ARGS+=("$1"); shift
      fi;;
  esac
done

[[ -n "$TASK" ]] || {
  echo "Usage: $0 <subcommand> [options]"
  echo "Subcommands: full_survey instant_q snp_support region_stats status validate"
  exit 1
}

# ── Load config ──
[[ -f "$CONFIG" ]] || { echo "[run_ancestry] Config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"

# ── v12.1: Source pipeline bridge for cross-module wiring ──
BRIDGE_SH="${BASE}/inversion_codebase_v8.5/utils/pipeline_bridge.sh"
if [[ -f "$BRIDGE_SH" ]]; then
  source "$BRIDGE_SH"
fi

# ── Dispatch ──
case "$TASK" in

  full_survey)
    anc_log "Engine A: Full Survey"
    STEP1_RUNNER="${BASE}/popstruct_thin/run_step1.sh"
    if [[ -f "$STEP1_RUNNER" ]]; then
      bash "$STEP1_RUNNER" "${EXTRA_ARGS[@]}"
    else
      anc_log "Step 1 runner not found at $STEP1_RUNNER"
      anc_log "Run Module 2B Step 1 scripts manually."
      exit 1
    fi
    ;;

  instant_q)
    anc_log "Engine B: Instant Q"
    exec "${SCRIPT_DIR}/instant_q" --config "$CONFIG" "${EXTRA_ARGS[@]}"
    ;;

  snp_support)
    anc_log "Module 2C: SNP Q Support"
    M2C_DIR="${SCRIPT_DIR}/module_2c"
    if [[ -d "$M2C_DIR" && -f "$M2C_DIR/run_module2c.sh" ]]; then
      bash "$M2C_DIR/run_module2c.sh" --config "$CONFIG" "${EXTRA_ARGS[@]}"
    else
      anc_log "Module 2C not found at $M2C_DIR"
      exit 1
    fi
    ;;

  region_stats)
    anc_log "Stats Dispatcher"
    python3 "${SCRIPT_DIR}/wrappers/instant_q.py" "${EXTRA_ARGS[@]}" \
      --config "$CONFIG"
    ;;

  status)
    echo "=== Unified Ancestry Module Status (v12.1) ==="
    echo ""

    echo "── Engine A (Full Survey) ──"
    if [[ -f "$BEST_QOPT" ]]; then
      K=$(awk '{print NF; exit}' "$BEST_QOPT")
      N=$(wc -l < "$BEST_QOPT")
      echo "  Best Q:  $BEST_QOPT (K=$K, N=$N samples)"
    else
      echo "  Best Q:  NOT FOUND"
    fi

    echo ""
    echo "── Engine B (Instant Q Cache) ──"
    echo "  Cache dir: ${LOCAL_Q_DIR}"
    if [[ -d "$LOCAL_Q_DIR" ]]; then
      N_CACHED=$(ls "${LOCAL_Q_DIR}"/*.local_Q_summary.* 2>/dev/null | wc -l)
      echo "  Cached chromosomes: ${N_CACHED} / ${N_CHROMOSOMES}"
    else
      echo "  Cache dir does not exist"
    fi

    echo ""
    echo "── Binaries ──"
    [[ -x "${INSTANT_Q_BIN:-}" ]] && echo "  instant_q: OK" || echo "  instant_q: NOT FOUND"
    [[ -x "${POPSTATS_BIN:-}" ]] && echo "  popstats: OK" || echo "  popstats: NOT FOUND"

    echo ""
    echo "── Sample Registry ──"
    if [[ -d "${REGISTRY_DIR:-}" ]]; then
      N_GROUPS=$(wc -l < "${REGISTRY_DIR}/sample_groups.tsv" 2>/dev/null || echo 0)
      echo "  Registry: ${REGISTRY_DIR} ($((N_GROUPS - 1)) groups)"
    else
      echo "  Registry: NOT INITIALIZED"
    fi

    echo ""
    echo "── Cross-module bridge ──"
    [[ -f "${LOAD_BRIDGE:-}" ]] && echo "  load_bridge.R: OK" || echo "  load_bridge.R: NOT FOUND"
    ;;

  validate)
    anc_log "Validating inputs..."
    STATUS=0
    check() {
      local label="$1" path="$2"
      if [[ -f "$path" || -d "$path" ]]; then
        echo "  [OK]   $label"
      else
        echo "  [FAIL] $label: $path"
        STATUS=1
      fi
    }

    check "Reference FAI"     "$REF_FAI"
    check "Sample list"       "$STEP1_SAMPLE_LIST"
    check "Samples.ind"       "${SAMPLES_IND:-NOT SET}"
    check "BEAGLE dir"        "$BEAGLE_DIR"
    check "Best Q (qopt)"     "$BEST_QOPT"
    check "Best F (fopt)"     "$BEST_FOPT"
    check "Load bridge"       "${LOAD_BRIDGE:-NOT SET}"
    check "Registry dir"      "${REGISTRY_DIR:-NOT SET}"

    [[ -x "${INSTANT_Q_BIN:-}" ]] && echo "  [OK]   C++ binary" || echo "  [WARN] instant_q not compiled"
    [[ -x "${POPSTATS_BIN:-}" ]] && echo "  [OK]   popstats binary" || echo "  [WARN] popstats not compiled"

    [[ $STATUS -eq 0 ]] && anc_log "All checks passed" || anc_log "Some checks failed"
    exit $STATUS
    ;;

  *)
    echo "[run_ancestry] Unknown task: $TASK" >&2
    exit 1
    ;;
esac
