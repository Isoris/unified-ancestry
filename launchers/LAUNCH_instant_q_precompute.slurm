#!/usr/bin/env bash
# =============================================================================
# LAUNCH_instant_q_precompute.slurm — SLURM array for Engine B (v13 K-sweep)
#
# 2026-04-17 rewrite:
#   - Array tasks map to (chrom × K) tuples, not chrom alone.
#   - Cache written to ${LOCAL_Q_DIR}/K${KK}/ (sharded per K level).
#   - Manifest row gains a K column.
#   - Per-K qinit/fopt resolved by swapping K${CANONICAL_K} → K${K} in the
#     BEST_QOPT/BEST_FOPT filename (matches the naming produced by
#     NGSadmix: thin500_K<NN>_seed<S>.qopt / .fopt.gz).
#   - Default sweep K=2..20. Override via K_SWEEP env ("2 3 8 12").
#   - If K_SWEEP contains a single value, behaves like the old single-K launcher.
#
# SIZING:
#   K values = 19 (full), chroms = 28, total = 532 array tasks.
#   Submit with the --array=1-532 below, or recompute N if K_SWEEP is changed.
# =============================================================================
#SBATCH --job-name=iq_precomp_K
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-532
#SBATCH --output=logs/iq_precomp_K_%a.out
#SBATCH --error=logs/iq_precomp_K_%a.err

set -euo pipefail

# ── Load config ──
CONFIG="${CONFIG:-$(dirname "${BASH_SOURCE[0]}")/../00_ancestry_config.sh}"
[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"

# Bridge (optional, for smap etc.)
BRIDGE_SH="${BRIDGE_SH:-${BASE}/utils/pipeline_bridge.sh}"
[[ -f "$BRIDGE_SH" ]] || BRIDGE_SH="${BASE}/inversion_modules/utils/pipeline_bridge.sh"
[[ -f "$BRIDGE_SH" ]] || BRIDGE_SH="${BASE}/inversion_codebase_v8.5/utils/pipeline_bridge.sh"
[[ -f "$BRIDGE_SH" ]] && source "$BRIDGE_SH"

eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-assembly}"

# ── Map array index to (chrom, K) tuple ──
mapfile -t CHROMS < <(awk '{print $1}' "$REF_FAI" | sort -V)
N_CHROMS=${#CHROMS[@]}

# K_SWEEP from env/config. Default covers K=2..20.
K_LIST=(${K_SWEEP:-2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20})
N_K=${#K_LIST[@]}
N_TOTAL=$(( N_CHROMS * N_K ))

# Array index 1..N — decompose into (K_idx, chrom_idx). Outer loop over K,
# inner over chrom, so a partial sweep of consecutive K_LIST entries runs
# all 28 chroms before moving on.
IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [[ $IDX -ge $N_TOTAL ]]; then
  echo "[LAUNCH] Array index ${SLURM_ARRAY_TASK_ID} exceeds ${N_TOTAL} — skip"
  exit 0
fi
K_IDX=$((IDX / N_CHROMS))
CHR_IDX=$((IDX % N_CHROMS))
K=${K_LIST[$K_IDX]}
CHR="${CHROMS[$CHR_IDX]}"

KK=$(printf "K%02d" "$K")
anc_log "Task ${SLURM_ARRAY_TASK_ID}/${N_TOTAL}: ${CHR} ${KK}"

# ── Resolve per-K qinit / fopt by filename swap ──
CANONICAL_KK=$(printf "K%02d" "${CANONICAL_K:-${DEFAULT_K:-8}}")
if [[ "$KK" == "$CANONICAL_KK" ]]; then
  Q_THIS="$BEST_QOPT"
  F_THIS="$BEST_FOPT"
else
  Q_THIS="${BEST_QOPT/${CANONICAL_KK}/${KK}}"
  F_THIS="${BEST_FOPT/${CANONICAL_KK}/${KK}}"
  # Also swap the seed suffix to whatever BEST_SEED_TABLE says is best for this K
  if [[ -f "${BEST_SEED_TABLE:-}" ]]; then
    SEED_K=$(awk -F'\t' -v k="$K" 'NR>1 && $1==k {print $2; exit}' "$BEST_SEED_TABLE" 2>/dev/null)
    if [[ -n "$SEED_K" ]]; then
      Q_THIS=$(echo "$Q_THIS" | sed -E "s/seed[0-9]+/seed${SEED_K}/")
      F_THIS=$(echo "$F_THIS" | sed -E "s/seed[0-9]+/seed${SEED_K}/")
    fi
  fi
fi
if [[ ! -f "$Q_THIS" ]]; then
  anc_log "qinit for $KK not found: $Q_THIS — skip this (K,chr)"
  exit 0
fi
if [[ ! -f "$F_THIS" ]]; then
  anc_log "fopt for $KK not found: $F_THIS — skip this (K,chr)"
  exit 0
fi

# ── Resolve BEAGLE ──
BGL=$(ls "${BEAGLE_DIR}"/*${CHR}*.beagle.gz 2>/dev/null | head -1)
if [[ -z "$BGL" || ! -f "$BGL" ]]; then
  anc_log "No BEAGLE for ${CHR} — skip"
  exit 0
fi

# ── Check binary ──
if [[ ! -x "$INSTANT_Q_BIN" ]]; then
  anc_log "Compiling instant_q..."
  cd "$(dirname "$INSTANT_Q_BIN")"
  make -j "${SLURM_CPUS_PER_TASK:-4}"
  cd -
fi

# ── Sample group — registered group_id for this cohort ──────────────────────
# 2026-04-17 chat-16: filename identity is now a registered group_id from
# sample_registry (e.g. "all_226"). Default "all_226" — the full 226-sample
# cohort. Override by exporting SAMPLE_GROUP before sbatch, e.g.
#   SAMPLE_GROUP=unrelated_81 sbatch LAUNCH_instant_q_precompute.slurm
#
# This replaces the chat-15 content-hash scheme that produced R-vs-bash tag
# mismatches. No hashing in this launcher anymore — the group name is a
# single source of truth with the R wrapper's get_active_sample_group().
#
# If the name needs to resolve to a registered group (which it does for the
# results_registry FK check later in the pipeline), validate here before we
# write files under a name that will eventually be rejected.
SAMPLE_GROUP="${SAMPLE_GROUP:-all_226}"
if [[ -n "${SAMPLE_REGISTRY:-}" && -f "${SAMPLE_REGISTRY}/groups/${SAMPLE_GROUP}.txt" ]]; then
  anc_log "Sample group: ${SAMPLE_GROUP} (registered in sample_registry)"
elif [[ -n "${SAMPLE_REGISTRY:-}" && -d "${SAMPLE_REGISTRY}/groups" ]]; then
  anc_log "Sample group: ${SAMPLE_GROUP} (NOT registered — will fail FK on downstream put_*)"
else
  anc_log "Sample group: ${SAMPLE_GROUP} (sample_registry not accessible from launcher)"
fi
# Back-compat variable name for older downstream scripts that read SAMPLE_SET_TAG.
SAMPLE_SET_TAG="$SAMPLE_GROUP"
export SAMPLE_GROUP SAMPLE_SET_TAG

# ── Sharded cache dir ──
CACHE_K="${LOCAL_Q_DIR}/${KK}"
mkdir -p "${CACHE_K}" logs

SUMMARY="${CACHE_K}/${CHR}.${SAMPLE_SET_TAG}.local_Q_summary.tsv.gz"
if [[ -f "$SUMMARY" && "${FORCE:-0}" -ne 1 ]]; then
  anc_log "Cache exists: ${SUMMARY} — skip (set FORCE=1 to recompute)"
  exit 0
fi

# ── Run precompute ──
TMPDIR_IQ=$(mktemp -d /tmp/iq_${CHR}_${KK}_XXXXXX)
anc_log "Running: ${CHR} ${KK} (win=${IQ_WINDOW_SIZE}, step=${IQ_WINDOW_STEP}, iter=${IQ_EM_ITER})"
T_START=$(date +%s)

"$INSTANT_Q_BIN" \
  --beagle "$BGL" \
  --fopt "$F_THIS" \
  --qinit "$Q_THIS" \
  --precompute \
  --outdir "$TMPDIR_IQ" \
  --chr "$CHR" \
  --window_size "$IQ_WINDOW_SIZE" \
  --window_step "$IQ_WINDOW_STEP" \
  --em_iter "$IQ_EM_ITER" \
  --ncores "${SLURM_CPUS_PER_TASK:-4}" \
  --sample_output

T_END=$(date +%s)
ELAPSED=$((T_END - T_START))

# ── Move to K-sharded cache with sample-set tag in filename ──
# instant_q emits <chr>.local_Q_summary.tsv etc. We rename to include the
# sample-set tag so the filename uniquely identifies (chr, K, sample_set).
for f in "${TMPDIR_IQ}/${CHR}".*; do
  [[ -f "$f" ]] || continue
  bn=$(basename "$f")
  # Insert .${SAMPLE_SET_TAG} between the chrom and the .local_Q_* suffix
  new_bn=$(echo "$bn" | sed -E "s/^(${CHR})(\\.local_Q_)/\\1.${SAMPLE_SET_TAG}\\2/")
  if [[ "$new_bn" == *summary.tsv || "$new_bn" == *samples.tsv ]]; then
    gzip -c "$f" > "${CACHE_K}/${new_bn}.gz"
  else
    cp "$f" "${CACHE_K}/${new_bn}"
  fi
done

rm -rf "$TMPDIR_IQ"
anc_log "Done: ${CHR} ${KK} in ${ELAPSED}s"

# ── Update manifest (single file at LOCAL_Q_DIR root; includes K + sample_set) ──
MANIFEST="${LOCAL_Q_DIR}/manifest.tsv"
if [[ ! -f "$MANIFEST" ]]; then
  echo -e "chrom\tK\tsample_set\tsummary_file\tsamples_file\ttimestamp" > "$MANIFEST"
fi
# Idempotent manifest upgrade: if the header lacks K or sample_set columns,
# rewrite with defaults (canonical K, N226_legacy tag). Bash check is a
# whole-word grep so a stray "K" in timestamps doesn't fool it.
HEADER=$(head -1 "$MANIFEST" 2>/dev/null || echo "")
NEED_UPGRADE=0
echo "$HEADER" | grep -Pq "(^|\t)K(\t|$)" || NEED_UPGRADE=1
echo "$HEADER" | grep -Pq "(^|\t)sample_set(\t|$)" || NEED_UPGRADE=1
if [[ $NEED_UPGRADE -eq 1 ]]; then
  TMP_M=$(mktemp)
  echo -e "chrom\tK\tsample_set\tsummary_file\tsamples_file\ttimestamp" > "$TMP_M"
  awk -F'\t' -v k="${CANONICAL_K:-${DEFAULT_K:-8}}" -v st="N226_legacy" \
    'NR>1 {
       # Support both legacy-flat (4 cols: chrom summary samples ts) and
       # chat-15-K-only (5 cols: chrom K summary samples ts) layouts.
       if (NF == 5) print $1 "\t" $2 "\t" st "\t" $3 "\t" $4 "\t" $5;
       else         print $1 "\t"  k "\t" st "\t" $2 "\t" $3 "\t" $4;
     }' "$MANIFEST" >> "$TMP_M"
  mv "$TMP_M" "$MANIFEST"
fi
echo -e "${CHR}\t${K}\t${SAMPLE_SET_TAG}\t${CACHE_K}/${CHR}.${SAMPLE_SET_TAG}.local_Q_summary.tsv.gz\t${CACHE_K}/${CHR}.${SAMPLE_SET_TAG}.local_Q_samples.tsv.gz\t$(date -Iseconds)" >> "$MANIFEST"

anc_log "Manifest updated"
