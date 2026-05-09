#!/usr/bin/env bash
# =============================================================================
# STEP_HH_A_build_subset_bamlists.sh — Generate per-group BAM lists (v12.1 REWIRED)
#
# v12.1 changes:
#   1. Can accept --from-registry to build bamlists from registry groups
#   2. Registers generated subsets back to the registry
#   3. Fallback to NGSadmix Q matrix if no registry
#
# Usage:
#   # Default: build from NGSadmix ancestry clusters
#   bash STEP_HH_A_build_subset_bamlists.sh
#
#   # From registry group IDs (e.g., inversion genotype classes)
#   bash STEP_HH_A_build_subset_bamlists.sh --from-registry "inv_LG05_HOM_REF,inv_LG05_HET,inv_LG05_HOM_INV"
#
#   # From registry pattern match
#   bash STEP_HH_A_build_subset_bamlists.sh --from-registry-pattern "inv_.*_HET"
# =============================================================================
set -euo pipefail

# Source hobs config (which sources ancestry config)
HOBS_CONFIG="$(dirname "${BASH_SOURCE[0]}")/00_hobs_hwe_config.sh"
if [[ -f "$HOBS_CONFIG" ]]; then
  source "$HOBS_CONFIG"
else
  source "$(dirname "${BASH_SOURCE[0]}")/../00_ancestry_config.sh"
fi

# Source pipeline bridge for registry access
BRIDGE_SH="${BASE}/inversion_codebase_v8.5/utils/pipeline_bridge.sh"
[[ -f "$BRIDGE_SH" ]] && source "$BRIDGE_SH" 2>/dev/null || true

FROM_REGISTRY=""
FROM_REGISTRY_PATTERN=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from-registry) FROM_REGISTRY="$2"; shift 2;;
    --from-registry-pattern) FROM_REGISTRY_PATTERN="$2"; shift 2;;
    *) shift;;
  esac
done

OUTDIR="${HOBS_OUTDIR}/subsets"
mkdir -p "$OUTDIR"

hobs_log "Building subset BAM lists..."

# ── Load BAM list and sample names ──
mapfile -t BAMS < "$BAMLIST"
N=${#BAMS[@]}

# Build sample name array from samples.ind (canonical CGA names)
SAMPLES_FILE="${SAMPLES_IND:-${BASE}/het_roh/01_inputs_check/samples.ind}"
if [[ -f "$SAMPLES_FILE" ]]; then
  mapfile -t SAMPLES < "$SAMPLES_FILE"
else
  # Fallback: extract from BAM paths
  SAMPLES=()
  for bam in "${BAMS[@]}"; do
    bn=$(basename "$bam")
    bn="${bn%.sorted.markdup.bam}"
    bn="${bn%.bam}"
    SAMPLES+=("$bn")
  done
fi

hobs_log "  $N samples"

# ── Load pruned list ──
declare -A IS_PRUNED
if [[ -f "${PRUNED_LIST:-}" ]]; then
  while IFS= read -r sid; do
    sid_clean=$(echo "$sid" | xargs)
    [[ -n "$sid_clean" ]] && IS_PRUNED["$sid_clean"]=1
  done < "$PRUNED_LIST"
  hobs_log "  Pruned list: ${#IS_PRUNED[@]} samples"
fi

# =============================================================================
# MODE 1: Build from registry groups
# =============================================================================

if [[ -n "$FROM_REGISTRY" || -n "$FROM_REGISTRY_PATTERN" ]]; then
  hobs_log "Building subsets from registry..."

  # Use R to resolve registry groups → sample lists
  "${RSCRIPT_BIN}" -e "
    source('${LOAD_BRIDGE}')
    
    from_ids <- '${FROM_REGISTRY}'
    from_pat <- '${FROM_REGISTRY_PATTERN}'
    outdir   <- '${OUTDIR}'
    
    all_groups <- reg\$list_groups()
    
    if (nzchar(from_ids)) {
      gids <- trimws(strsplit(from_ids, ',')[[1]])
    } else if (nzchar(from_pat)) {
      gids <- all_groups[grepl(from_pat, group_id), group_id]
    } else {
      stop('No registry source specified')
    }
    
    cat('[registry→bamlists] Resolving', length(gids), 'groups\n')
    
    # Load BAM list for index mapping
    bams <- readLines('${BAMLIST}')
    samples <- readLines('${SAMPLES_FILE}')
    stopifnot(length(bams) == length(samples))
    
    manifest_rows <- list()
    
    for (gid in gids) {
      members <- reg\$get_group(gid)
      if (length(members) == 0) {
        cat('  SKIP:', gid, '(empty)\n')
        next
      }
      
      # Map CGA sample IDs → BAM paths
      idx <- match(members, samples)
      idx <- idx[!is.na(idx)]
      if (length(idx) == 0) next
      
      bam_subset <- bams[idx]
      
      # Write bamlist and samples file
      bamlist_file <- file.path(outdir, paste0(gid, '.bamlist'))
      samples_file <- file.path(outdir, paste0(gid, '.samples'))
      writeLines(bam_subset, bamlist_file)
      writeLines(members[members %in% samples], samples_file)
      
      cat('  ', gid, ':', length(idx), 'samples\n')
      
      manifest_rows[[length(manifest_rows) + 1]] <- data.table::data.table(
        subset_id = gid,
        subset_type = 'registry_group',
        n_samples = length(idx),
        bamlist_path = bamlist_file,
        samples_path = samples_file,
        notes = paste0('From registry: ', gid)
      )
    }
    
    if (length(manifest_rows) > 0) {
      manifest <- data.table::rbindlist(manifest_rows)
      manifest_file <- file.path(outdir, 'subset_manifest_registry.tsv')
      data.table::fwrite(manifest, manifest_file, sep = '\t')
      cat('[registry→bamlists] Manifest:', manifest_file, '\n')
    }
  "

  hobs_log "Registry-based subsets done"
  exit 0
fi

# =============================================================================
# MODE 2: Build from NGSadmix ancestry clusters (original behavior)
# =============================================================================

hobs_log "Building subsets from NGSadmix K=$DEFAULT_K ancestry clusters..."

# Assign each sample to its dominant cluster
python3 - "$BEST_QOPT" "$DEFAULT_K" "$N" << 'PYEOF' > "${OUTDIR}/sample_assignments.tsv"
import sys
qopt, K, N = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
with open(qopt) as f:
    for i, line in enumerate(f):
        if i >= N: break
        vals = [float(x) for x in line.strip().split()]
        if len(vals) != K: continue
        assigned = vals.index(max(vals)) + 1
        print(f"{i}\t{assigned}")
PYEOF

# Initialize files
for k in $(seq 1 "$DEFAULT_K"); do
  > "${OUTDIR}/Q${k}_full.bamlist"
  > "${OUTDIR}/Q${k}_pruned.bamlist"
  > "${OUTDIR}/Q${k}_full.samples"
  > "${OUTDIR}/Q${k}_pruned.samples"
done

declare -A GROUP_COUNTS_FULL
declare -A GROUP_COUNTS_PRUNED
for k in $(seq 1 "$DEFAULT_K"); do
  GROUP_COUNTS_FULL[$k]=0
  GROUP_COUNTS_PRUNED[$k]=0
done

while IFS=$'\t' read -r idx cluster; do
  [[ -z "$idx" ]] && continue
  i=$((idx))
  sid="${SAMPLES[$i]}"
  bam="${BAMS[$i]}"

  echo "$bam" >> "${OUTDIR}/Q${cluster}_full.bamlist"
  echo "$sid" >> "${OUTDIR}/Q${cluster}_full.samples"
  GROUP_COUNTS_FULL[$cluster]=$(( ${GROUP_COUNTS_FULL[$cluster]} + 1 ))

  if [[ -n "${IS_PRUNED[$sid]:-}" ]]; then
    echo "$bam" >> "${OUTDIR}/Q${cluster}_pruned.bamlist"
    echo "$sid" >> "${OUTDIR}/Q${cluster}_pruned.samples"
    GROUP_COUNTS_PRUNED[$cluster]=$(( ${GROUP_COUNTS_PRUNED[$cluster]} + 1 ))
  fi
done < "${OUTDIR}/sample_assignments.tsv"

# Full panel
cp "$BAMLIST" "${OUTDIR}/allSamples_referenceOnly.bamlist"

# Write manifest
MANIFEST="${OUTDIR}/subset_manifest.tsv"
{
  echo -e "subset_id\tsubset_type\tn_samples\tbamlist_path\tsamples_path\tnotes"
  for k in $(seq 1 "$DEFAULT_K"); do
    nf=${GROUP_COUNTS_FULL[$k]}
    np=${GROUP_COUNTS_PRUNED[$k]}
    echo -e "Q${k}_full\tancestry_cluster_full\t${nf}\t${OUTDIR}/Q${k}_full.bamlist\t${OUTDIR}/Q${k}_full.samples\tK=${DEFAULT_K} cluster ${k}"
    if [[ $np -gt 0 ]]; then
      echo -e "Q${k}_pruned\tancestry_cluster_pruned\t${np}\t${OUTDIR}/Q${k}_pruned.bamlist\t${OUTDIR}/Q${k}_pruned.samples\tK=${DEFAULT_K} cluster ${k} pruned"
    fi
  done
  echo -e "allSamples_referenceOnly\tfull_panel_confounded\t${N}\t${OUTDIR}/allSamples_referenceOnly.bamlist\t${SAMPLE_LIST}\tWARNING: confounded by structure+relatedness"
} > "$MANIFEST"

# v12.1: Register hobs subsets to registry
if [[ -f "${LOAD_BRIDGE}" ]]; then
  "${RSCRIPT_BIN}" -e "
    source('${LOAD_BRIDGE}')
    '%+%' <- function(a, b) paste0(a, b)
    for (k in 1:${DEFAULT_K}) {
      sf <- '${OUTDIR}/Q' %+% k %+% '_full.samples'
      if (file.exists(sf)) {
        members <- readLines(sf)
        members <- members[nzchar(members)]
        if (length(members) > 0) {
          register_group(
            paste0('hobs_Q', k, '_full'),
            members,
            description = paste0('Hobs/HWE subset: K=${DEFAULT_K} cluster ', k, ' (full)')
          )
        }
      }
    }
  " 2>/dev/null || true
fi

hobs_log "Manifest: $MANIFEST"
cat "$MANIFEST" | column -t
hobs_log "Done"
