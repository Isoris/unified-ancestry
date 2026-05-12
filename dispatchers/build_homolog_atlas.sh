#!/usr/bin/env bash
# =============================================================================
# build_homolog_atlas.sh — full "click a gene → see all chains/paralogs" pipeline.
#
# Steps:
#   1. For each --target genome=label: run miniprot in the chosen mode → PAF.
#   2. Optionally: run DIAMOND blastp on the query proteome against itself
#      (or a supplied subject) → TSV.
#   3. Build the combined homolog index from all PAFs + the DIAMOND TSV.
#
# Output: <out_dir>/homolog_atlas.holindx (mmap-friendly, microsecond lookups
# via engines/homolog_index_query).
#
# Configuration: command line, or --config <atlas_config.json>
# (schema: engines/schemas/homolog_atlas_config.schema.json).
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENGINES_DIR="$(cd "$SCRIPT_DIR/../engines" && pwd)"
RUN_MINIPROT="$SCRIPT_DIR/run_miniprot_chain.sh"
HIB_BIN="$ENGINES_DIR/homolog_index_build"

usage() {
  cat <<'EOF'
Usage: build_homolog_atlas.sh --query <protein.fa> --out_dir <dir>
                              --target <genome.fa>=<label> [--target ...]
                              [--diamond <hits.tsv>=<label>]   (may repeat)
                              [--run_diamond_blastp <subject.fa>=<label>]
                                  (runs DIAMOND blastp inside; needs `diamond` on PATH)
                              [--mode strict|chain|loose]      (default chain)
                              [--threads N]                    (default 8)
                              [--keep_intermediate]
                              [--dry_run]
                              [--config <atlas_config.json>]   (alternative to CLI)

Outputs (in --out_dir):
  homolog_atlas.holindx     — combined binary index
  miniprot/<label>.paf      — one PAF per --target
  diamond/<label>.tsv       — one TSV per --run_diamond_blastp
  build.log                 — stderr capture

After build, query like:
  engines/homolog_index_query --index out_dir/homolog_atlas.holindx \
                              --gene MY_GENE --format json
EOF
}

# ── State holders ───────────────────────────────────────────────────────────

QUERY=""; OUT_DIR=""; MODE="chain"; THREADS=8
KEEP_INTERMEDIATE=0; DRY_RUN=0; CONFIG=""
declare -a TARGETS=()         # entries: genome.fa=label
declare -a DIAMOND_FILES=()   # entries: hits.tsv=label
declare -a DIAMOND_SUBJECTS=() # entries: subject.fa=label (will run diamond blastp)

# ── Parse CLI ───────────────────────────────────────────────────────────────

while [[ $# -gt 0 ]]; do
  case "$1" in
    --query)              QUERY=$2; shift 2 ;;
    --out_dir)            OUT_DIR=$2; shift 2 ;;
    --target)             TARGETS+=("$2"); shift 2 ;;
    --diamond)            DIAMOND_FILES+=("$2"); shift 2 ;;
    --run_diamond_blastp) DIAMOND_SUBJECTS+=("$2"); shift 2 ;;
    --mode)               MODE=$2; shift 2 ;;
    --threads)            THREADS=$2; shift 2 ;;
    --keep_intermediate)  KEEP_INTERMEDIATE=1; shift ;;
    --dry_run)            DRY_RUN=1; shift ;;
    --config)             CONFIG=$2; shift 2 ;;
    -h|--help)            usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

# ── Optional config (JSON; parsed with python if available) ─────────────────

if [[ -n "$CONFIG" ]]; then
  command -v python3 >/dev/null || { echo "--config needs python3 on PATH" >&2; exit 2; }
  # Pull config into shell variables. The script's CLI args override config.
  eval "$(python3 - "$CONFIG" <<'PY'
import json, sys, shlex, os
cfg = json.load(open(sys.argv[1]))
out = []
if "query" in cfg and not os.environ.get("CFG_QUERY_OVERRIDE"):
    out.append(f"QUERY=${{QUERY:-{shlex.quote(cfg['query'])}}}")
if "out_dir" in cfg:
    out.append(f"OUT_DIR=${{OUT_DIR:-{shlex.quote(cfg['out_dir'])}}}")
if "mode" in cfg:
    out.append(f"MODE=${{MODE:-{shlex.quote(cfg['mode'])}}}")
if "threads" in cfg:
    out.append(f"THREADS=${{THREADS:-{cfg['threads']}}}")
for t in cfg.get("targets", []):
    out.append(f"TARGETS+=({shlex.quote(t['path']+'='+t['label'])})")
for d in cfg.get("diamond_files", []):
    out.append(f"DIAMOND_FILES+=({shlex.quote(d['path']+'='+d['label'])})")
for d in cfg.get("diamond_subjects", []):
    out.append(f"DIAMOND_SUBJECTS+=({shlex.quote(d['path']+'='+d['label'])})")
print("\n".join(out))
PY
)"
fi

# ── Validate ────────────────────────────────────────────────────────────────

[[ -z "$QUERY"   ]] && { echo "--query required"   >&2; usage; exit 1; }
[[ -z "$OUT_DIR" ]] && { echo "--out_dir required" >&2; usage; exit 1; }
[[ ${#TARGETS[@]} -eq 0 && ${#DIAMOND_FILES[@]} -eq 0 && ${#DIAMOND_SUBJECTS[@]} -eq 0 ]] && {
  echo "Need at least one --target, --diamond, or --run_diamond_blastp" >&2; exit 1;
}

[[ -x "$HIB_BIN" ]] || { echo "homolog_index_build not built (run 'make -C $ENGINES_DIR')" >&2; exit 3; }

mkdir -p "$OUT_DIR/miniprot" "$OUT_DIR/diamond"
LOG="$OUT_DIR/build.log"
: > "$LOG"

run() {
  echo "+ $*" | tee -a "$LOG" >&2
  [[ "$DRY_RUN" -eq 1 ]] && return 0
  "$@" 2>>"$LOG"
}

# ── miniprot passes ─────────────────────────────────────────────────────────

declare -a HIB_ARGS=()

for t in "${TARGETS[@]}"; do
  GENOME="${t%%=*}"; LABEL="${t##*=}"
  [[ "$LABEL" == "$GENOME" ]] && LABEL="$(basename "$GENOME" | sed -E 's/\.(fa|fasta|fna|gz)//g')"
  PAF="$OUT_DIR/miniprot/${LABEL}.paf"
  echo "[atlas] miniprot mode=$MODE label=$LABEL → $PAF" | tee -a "$LOG" >&2
  run "$RUN_MINIPROT" --query "$QUERY" --target "$GENOME" --out "$PAF" \
                      --mode "$MODE" --threads "$THREADS" --label "$LABEL"
  HIB_ARGS+=( --miniprot "$PAF=$LABEL" )
done

# ── DIAMOND passes ──────────────────────────────────────────────────────────

for s in "${DIAMOND_SUBJECTS[@]}"; do
  SUBJ="${s%%=*}"; LABEL="${s##*=}"
  [[ "$LABEL" == "$SUBJ" ]] && LABEL="$(basename "$SUBJ" | sed -E 's/\.(fa|fasta)//g')"
  DB="$OUT_DIR/diamond/${LABEL}.dmnd"
  TSV="$OUT_DIR/diamond/${LABEL}.tsv"
  command -v diamond >/dev/null || { echo "diamond not on PATH (needed for --run_diamond_blastp)" >&2; exit 4; }
  echo "[atlas] diamond makedb label=$LABEL → $DB" | tee -a "$LOG" >&2
  run diamond makedb --in "$SUBJ" -d "$DB" --threads "$THREADS" --quiet
  echo "[atlas] diamond blastp $QUERY → $TSV" | tee -a "$LOG" >&2
  run diamond blastp -q "$QUERY" -d "$DB" -o "$TSV" --threads "$THREADS" \
                     --outfmt 6 --quiet --evalue 1e-5
  HIB_ARGS+=( --diamond "$TSV=$LABEL" )
done

for d in "${DIAMOND_FILES[@]}"; do
  TSV="${d%%=*}"; LABEL="${d##*=}"
  [[ "$LABEL" == "$TSV" ]] && LABEL="$(basename "$TSV" | sed -E 's/\.(tsv|txt)//g')"
  HIB_ARGS+=( --diamond "$TSV=$LABEL" )
done

# ── Build index ─────────────────────────────────────────────────────────────

INDEX_OUT="$OUT_DIR/homolog_atlas.holindx"
echo "[atlas] building combined index → $INDEX_OUT" | tee -a "$LOG" >&2
run "$HIB_BIN" "${HIB_ARGS[@]}" --out "$INDEX_OUT"

# ── Cleanup ─────────────────────────────────────────────────────────────────

if [[ "$KEEP_INTERMEDIATE" -eq 0 && "$DRY_RUN" -eq 0 ]]; then
  echo "[atlas] removing intermediates (use --keep_intermediate to keep)" | tee -a "$LOG" >&2
  rm -rf "$OUT_DIR/miniprot" "$OUT_DIR/diamond"
fi

echo "[atlas] done. Query with:" | tee -a "$LOG" >&2
echo "  $ENGINES_DIR/homolog_index_query --index $INDEX_OUT --gene <ID> --format json" | tee -a "$LOG" >&2
