#!/usr/bin/env bash
# =============================================================================
# run_miniprot_chain.sh — single miniprot run with chain-context presets.
#
# Replaces the older strict-only dispatcher. Three modes:
#   strict  — single best hit per query, high confidence (orthologue picking).
#   chain   — primary + paralog chains, gene-context loose. Default.
#   loose   — wide net for paralog / homolog discovery.
#
# Output: PAF (standard miniprot format). Pipe into homolog_index_build, or
# pass --build_index <out.holindx> to do both in one shot.
#
# Required tools on PATH: miniprot.
# Optional: engines/homolog_index_build (resolved relative to this script).
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENGINES_DIR="$(cd "$SCRIPT_DIR/../engines" 2>/dev/null && pwd || echo "$SCRIPT_DIR/../engines")"
HIB_BIN="$ENGINES_DIR/homolog_index_build"

usage() {
  cat <<'EOF'
Usage: run_miniprot_chain.sh --query <protein.fa> --target <genome.fa>
                              --out <out.paf>
                              [--mode strict|chain|loose] (default chain)
                              [--threads N]              (default 4)
                              [--label <name>]           (default basename of --target)
                              [--extra '<extra miniprot flags>']
                              [--build_index <out.holindx>]
                              [--gff <out.gff>]
                              [--dry_run]

Mode presets (override with --extra):
  strict  : --outs=0.99  -N 1   -p 0.5  -G 100k
  chain   : --outs=0.95  -N 20  -p 0.3  -G 100k
  loose   : --outs=0.80  -N 100 -p 0.2  -G 200k

Notes:
  - For repeat-heavy fish/teleost genomes, use --extra '--no-frs' to disable
    repeat filtering on a per-query basis (rarely needed in chain mode).
  - --gff also runs a second miniprot pass writing GFF3 (slow; skip unless
    you need it for downstream gene-model code).
  - --build_index calls homolog_index_build on the PAF immediately.
EOF
}

QUERY=""; TARGET=""; OUT=""; MODE="chain"; THREADS=4; LABEL=""
EXTRA=""; BUILD_INDEX=""; GFF=""; DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --query)        QUERY=$2; shift 2 ;;
    --target)       TARGET=$2; shift 2 ;;
    --out)          OUT=$2; shift 2 ;;
    --mode)         MODE=$2; shift 2 ;;
    --threads)      THREADS=$2; shift 2 ;;
    --label)        LABEL=$2; shift 2 ;;
    --extra)        EXTRA=$2; shift 2 ;;
    --build_index)  BUILD_INDEX=$2; shift 2 ;;
    --gff)          GFF=$2; shift 2 ;;
    --dry_run)      DRY_RUN=1; shift ;;
    -h|--help)      usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "$QUERY"  ]] && { echo "--query required"  >&2; usage; exit 1; }
[[ -z "$TARGET" ]] && { echo "--target required" >&2; usage; exit 1; }
[[ -z "$OUT"    ]] && { echo "--out required"    >&2; usage; exit 1; }
[[ -z "$LABEL"  ]] && LABEL="$(basename "$TARGET" | sed -E 's/\.(fa|fasta|fna|gz)//g')"

case "$MODE" in
  strict) PRESET="--outs=0.99 -N 1   -p 0.5 -G 100k" ;;
  chain)  PRESET="--outs=0.95 -N 20  -p 0.3 -G 100k" ;;
  loose)  PRESET="--outs=0.80 -N 100 -p 0.2 -G 200k" ;;
  *) echo "Unknown --mode: $MODE (use strict|chain|loose)" >&2; exit 1 ;;
esac

mkdir -p "$(dirname "$OUT")"

CMD=(miniprot -t "$THREADS" $PRESET)
[[ -n "$EXTRA" ]] && CMD+=( $EXTRA )
CMD+=( "$TARGET" "$QUERY" )

echo "[run_miniprot_chain] mode=$MODE label=$LABEL threads=$THREADS" >&2
echo "[run_miniprot_chain] cmd: ${CMD[*]} > $OUT" >&2

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "[run_miniprot_chain] DRY RUN — exiting before execution" >&2
  exit 0
fi

command -v miniprot >/dev/null || { echo "miniprot not on PATH" >&2; exit 2; }

"${CMD[@]}" > "$OUT"
echo "[run_miniprot_chain] PAF written: $OUT ($(wc -l < "$OUT") lines)" >&2

if [[ -n "$GFF" ]]; then
  GFF_CMD=(miniprot -t "$THREADS" $PRESET --gff)
  [[ -n "$EXTRA" ]] && GFF_CMD+=( $EXTRA )
  GFF_CMD+=( "$TARGET" "$QUERY" )
  echo "[run_miniprot_chain] GFF pass: ${GFF_CMD[*]} > $GFF" >&2
  "${GFF_CMD[@]}" > "$GFF"
fi

if [[ -n "$BUILD_INDEX" ]]; then
  [[ -x "$HIB_BIN" ]] || { echo "homolog_index_build not found at $HIB_BIN (run 'make -C $ENGINES_DIR')" >&2; exit 3; }
  echo "[run_miniprot_chain] building index → $BUILD_INDEX" >&2
  "$HIB_BIN" --miniprot "$OUT=$LABEL" --out "$BUILD_INDEX"
fi
