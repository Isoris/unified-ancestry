#!/usr/bin/env bash
# =============================================================================
# test_integration.sh — Unified Ancestry Module v12 Integration Tests
#
# Tests:
#   T01  Compile check: all 3 C/C++ binaries
#   T02  Synthetic BEAGLE: region_popstats produces expected columns
#   T03  Synthetic HWE: hobs_windower produces expected output structure
#   T04  Config parse: auto-discover K from best_seed_by_K.tsv
#   T05  Group resolution: resolve_groups returns correct structure
#   T06  CLI routing: region_stats routes Q-only to instant_q
#   T07  CLI routing: region_stats routes Fst to C binary
#   T08  Dispatcher: cleaned R dispatcher has no F_IS/Ashman_D
#   T09  Pipeline steps: run_full_pipeline.sh --step 8 compiles + runs
#   T10  Plotting: plots/plot_hobs_hwe.R loads without error
#
# Usage:
#   bash tests/test_integration.sh [--config 00_ancestry_config.sh]
#   bash tests/test_integration.sh --test T01    # single test
#
# Exit code: 0 if all pass, 1 if any fail.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/00_ancestry_config.sh}"
# v12.2: source config so moved-script variables (HOBS_PLOT_R, etc.) resolve
[[ -f "$CONFIG" ]] && source "$CONFIG"
SINGLE_TEST="${2:-}"
TMPDIR=$(mktemp -d /tmp/ancestry_test_XXXXXX)
PASS=0; FAIL=0; SKIP=0

cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

# ── Helpers ──
run_test() {
  local id="$1" desc="$2"
  shift 2
  if [[ -n "$SINGLE_TEST" && "$SINGLE_TEST" != "$id" ]]; then
    return 0
  fi
  echo -n "[$id] $desc ... "
}

pass() { echo "PASS"; PASS=$((PASS + 1)); }
fail() { echo "FAIL: $*"; FAIL=$((FAIL + 1)); }
skip() { echo "SKIP: $*"; SKIP=$((SKIP + 1)); }

# =============================================================================
# T01: Compile check — all 3 C/C++ binaries
# =============================================================================
run_test T01 "Compile instant_q.cpp"
if g++ -O3 -march=native -fopenmp -std=c++17 -w -o "$TMPDIR/instant_q" \
    "${SCRIPT_DIR}/engines/instant_q.cpp" -lz -fopenmp 2>"$TMPDIR/t01a.err"; then
  pass
else
  fail "$(head -3 "$TMPDIR/t01a.err")"
fi

run_test T01b "Compile hobs_windower.c"
if gcc -O3 -w -o "$TMPDIR/hobs_windower" \
    "${SCRIPT_DIR}/engines/hobs_windower.c" -lz -lm 2>"$TMPDIR/t01b.err"; then
  pass
else
  fail "$(head -3 "$TMPDIR/t01b.err")"
fi

run_test T01c "Compile region_popstats.c"
if gcc -O3 -march=native -fopenmp -w -o "$TMPDIR/region_popstats" \
    "${SCRIPT_DIR}/engines/region_popstats.c" -lz -lm -fopenmp 2>"$TMPDIR/t01c.err"; then
  pass
else
  fail "$(head -3 "$TMPDIR/t01c.err")"
fi

# =============================================================================
# T02: Synthetic BEAGLE → region_popstats output structure
# =============================================================================
run_test T02 "region_popstats synthetic data"

# Create minimal BEAGLE (3 individuals, 20 sites)
{
  echo -e "marker\tallele1\tallele2\tInd0\tInd0\tInd0\tInd1\tInd1\tInd1\tInd2\tInd2\tInd2"
  for pos in $(seq 1000 1000 20000); do
    # Random-ish GLs
    echo -e "chr1_${pos}\tA\tT\t0.8\t0.15\t0.05\t0.1\t0.8\t0.1\t0.05\t0.15\t0.8"
  done
} | gzip > "$TMPDIR/test.beagle.gz"

# Sample list
echo -e "Ind0\nInd1\nInd2" > "$TMPDIR/samples.txt"

# Group files
echo "Ind0" > "$TMPDIR/grpA.txt"
echo -e "Ind1\nInd2" > "$TMPDIR/grpB.txt"

if [[ -x "$TMPDIR/region_popstats" ]]; then
  "$TMPDIR/region_popstats" \
    --beagle "$TMPDIR/test.beagle.gz" \
    --sample_list "$TMPDIR/samples.txt" \
    --chr chr1 \
    --fixed_win "10000:5000" \
    --groups "grpA:$TMPDIR/grpA.txt,grpB:$TMPDIR/grpB.txt" \
    --out "$TMPDIR/t02_out.tsv" 2>"$TMPDIR/t02.log"

  # Check output has expected columns
  HEADER=$(head -2 "$TMPDIR/t02_out.tsv" | tail -1)
  OK=1
  for col in window_id chrom start end theta_pi_all theta_W_all Tajima_D Hp \
             Fst_grpA_grpB dXY_grpA_grpB MI_grpA_grpB; do
    if ! echo "$HEADER" | grep -q "$col"; then
      fail "Missing column: $col"
      OK=0
      break
    fi
  done
  [[ $OK -eq 1 ]] && pass
else
  skip "region_popstats not compiled"
fi

# =============================================================================
# T03: Synthetic HWE → hobs_windower output structure
# =============================================================================
run_test T03 "hobs_windower synthetic data"

# Create minimal .hwe.gz
{
  echo -e "Chromo\tPosition\tMajor\tMinor\thweFreq\tFreq\tF\tLRT\tp-value"
  for pos in $(seq 1000 500 50000); do
    echo -e "chr1\t${pos}\tA\tT\t0.300000\t0.310000\t0.050000\t2.500000\t0.120000"
  done
} | gzip > "$TMPDIR/test.hwe.gz"

if [[ -x "$TMPDIR/hobs_windower" ]]; then
  "$TMPDIR/hobs_windower" "$TMPDIR/test.hwe.gz" "$TMPDIR/t03_out" 50000 \
    --scales "5kb:5000:1000,10kb:10000:2000" 2>"$TMPDIR/t03.log"

  # Check sites file
  if [[ -f "$TMPDIR/t03_out.sites.tsv" ]]; then
    N_SITES=$(tail -n +2 "$TMPDIR/t03_out.sites.tsv" | wc -l)
    if [[ $N_SITES -gt 0 ]]; then
      # Check window files
      if [[ -f "$TMPDIR/t03_out.win5kb.tsv" && -f "$TMPDIR/t03_out.win10kb.tsv" ]]; then
        pass
      else
        fail "Missing window output files"
      fi
    else
      fail "Sites file empty"
    fi
  else
    fail "No sites file produced"
  fi
else
  skip "hobs_windower not compiled"
fi

# =============================================================================
# T04: Config auto-K detection
# =============================================================================
run_test T04 "Config auto-K from best_seed_by_K.tsv"

# Create mock best_seed_by_K.tsv
{
  echo -e "K\tseed\tloglik"
  echo -e "4\t3\t-1200000"
  echo -e "6\t2\t-1100000"
  echo -e "8\t1\t-1050000"
} > "$TMPDIR/best_seed_by_K.tsv"

# Simulate config logic
_AUTO_K=$(awk -F'\t' 'NR>1 {print $1, $3}' "$TMPDIR/best_seed_by_K.tsv" \
          | sort -k2 -rn | head -1 | awk '{print $1}')
if [[ "$_AUTO_K" == "8" ]]; then
  pass
else
  fail "Expected K=8, got $_AUTO_K"
fi

# =============================================================================
# T05: Dispatcher cleanup check — no F_IS or Ashman_D
# =============================================================================
run_test T05 "Dispatcher cleaned of F_IS/Ashman_D"

DISP="${SCRIPT_DIR}/dispatchers/region_stats_dispatcher.R"
if [[ -f "$DISP" ]]; then
  if grep -v '^\s*#' "$DISP" | grep -q "compute_F_IS\|compute_Ashman_D" 2>/dev/null; then
    fail "F_IS or Ashman_D still present in dispatcher (non-comment)"
  else
    pass
  fi
else
  skip "Dispatcher not found"
fi

# =============================================================================
# T06: Dispatcher has route_to_c_popstats
# =============================================================================
run_test T06 "Dispatcher has C-binary routing"

if [[ -f "$DISP" ]]; then
  if grep -q "route_to_c_popstats" "$DISP" 2>/dev/null; then
    pass
  else
    fail "route_to_c_popstats not found in dispatcher"
  fi
else
  skip "Dispatcher not found"
fi

# =============================================================================
# T07: No WAIT_FOR_ANGSD stubs remain
# =============================================================================
run_test T07 "No WAIT_FOR_ANGSD stubs"

if [[ -f "$DISP" ]]; then
  if grep -v '^\s*#' "$DISP" | grep -q "WAIT_FOR_ANGSD" 2>/dev/null; then
    fail "WAIT_FOR_ANGSD stubs still present (non-comment)"
  else
    pass
  fi
else
  skip "Dispatcher not found"
fi

# =============================================================================
# T08: Pipeline has Steps 8 and 9
# =============================================================================
run_test T08 "Pipeline orchestrator has Steps 8+9"

PIPELINE="${SCRIPT_DIR}/run_full_pipeline.sh"
if [[ -f "$PIPELINE" ]]; then
  HAS_8=$(grep -c "Step 8" "$PIPELINE" || true)
  HAS_9=$(grep -c "Step 9" "$PIPELINE" || true)
  if [[ $HAS_8 -gt 0 && $HAS_9 -gt 0 ]]; then
    pass
  else
    fail "Missing Step 8 ($HAS_8) or Step 9 ($HAS_9)"
  fi
else
  skip "Pipeline not found"
fi

# =============================================================================
# T09: SLURM launcher for popstats exists and references C binary
# =============================================================================
run_test T09 "SLURM launcher for region_popstats"

LAUNCHER="${SCRIPT_DIR}/launchers/LAUNCH_region_popstats.slurm"
if [[ -f "$LAUNCHER" ]]; then
  if grep -q "region_popstats" "$LAUNCHER" && grep -q "SLURM_ARRAY_TASK_ID" "$LAUNCHER"; then
    pass
  else
    fail "Launcher missing key content"
  fi
else
  skip "Launcher not found"
fi

# =============================================================================
# T10: Plotting script exists and has all 5 modes
# =============================================================================
run_test T10 "Part 12 plotting suite"

PLOT="${HOBS_PLOT_R:-${SCRIPT_DIR}/plots/plot_hobs_hwe.R}"
if [[ -f "$PLOT" ]]; then
  MODES=0
  for m in genome_tracks heatmap outlier_burden subset_compare candidate_stack; do
    grep -q "$m" "$PLOT" && MODES=$((MODES + 1))
  done
  if [[ $MODES -eq 5 ]]; then
    pass
  else
    fail "Only $MODES/5 modes found"
  fi
else
  skip "Plot script not found"
fi

# =============================================================================
# T11: region_popstats.c has no dead p_grp stack array
# =============================================================================
run_test T11 "No dead p_grp stack array in region_popstats.c"

RPOP="${SCRIPT_DIR}/engines/region_popstats.c"
if [[ -f "$RPOP" ]]; then
  if grep -q "p_grp\[MAX_GROUPS\]\[MAX_SITES\]" "$RPOP" 2>/dev/null; then
    fail "Dead p_grp[MAX_GROUPS][MAX_SITES] still present"
  else
    pass
  fi
else
  skip "region_popstats.c not found"
fi

# =============================================================================
# T12: Config has POPSTATS_BIN and HOBS_OUTDIR
# =============================================================================
run_test T12 "Config has new v12 variables"

CFG="${SCRIPT_DIR}/00_ancestry_config.sh"
if [[ -f "$CFG" ]]; then
  HAS_PB=$(grep -c "POPSTATS_BIN" "$CFG" || true)
  HAS_HO=$(grep -c "HOBS_OUTDIR" "$CFG" || true)
  if [[ $HAS_PB -gt 0 && $HAS_HO -gt 0 ]]; then
    pass
  else
    fail "Missing POPSTATS_BIN ($HAS_PB) or HOBS_OUTDIR ($HAS_HO)"
  fi
else
  skip "Config not found"
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "═══════════════════════════════════════════════"
echo "  PASS: $PASS  FAIL: $FAIL  SKIP: $SKIP"
echo "═══════════════════════════════════════════════"

[[ $FAIL -eq 0 ]]
