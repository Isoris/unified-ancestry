#!/usr/bin/env python3
"""
internal_ancestry_composition.py — generic interval-internal ancestry classifier.

Scope
-----
Given a parent interval (chromosome, MDS candidate, inversion candidate, any
arbitrary region) and per-sample ancestry labels on windows inside it, score
whether the interval's internal ancestry composition is simple or structured.

Algorithm
---------
For each parent region, finds overlapping Engine B windows and characterizes
the internal Q composition per sample:
  - Dominant ancestry label blocks and switches
  - Fragmentation score
  - Internal entropy
  - Structure classification:
    homogeneous | dominant_plus_secondary | two_block_composite |
    continuous_gradient | multi_block_fragmented | diffuse_mixed

What "nested" used to mean (and why this file is no longer called that)
-----------------------------------------------------------------------
Previously `nested_composition.py`. The word "nested" described the data
structure (Q windows nested inside parent intervals) but in inversion
biology it implies "nested inversion" (inversion-inside-inversion), which
is a STRONGER biological claim than this analysis supports. This engine
measures internal ancestry composition only. Interpreting a "composite"
verdict as a nested inversion requires cross-checks with SV calls,
breakpoints, and karyotype evidence — that interpretation lives in the
pipeline wrapper (STEP_C01i_c_nested_composition.py in phase_7_karyotype_groups),
not here.

Consumes: Engine B per-sample local Q cache (*.local_Q_samples.tsv.gz)
Writes (standalone mode):
    nested_composition.tsv          (per-parent × per-sample)
    nested_composition_summary.tsv  (per-parent, sample-averaged)
    [output filenames kept for back-compat with run_full_pipeline.sh Step 4]

Importable API
--------------
This file is both a CLI tool and a library. Other scripts — notably the
Phase 7 inversion wrapper STEP_C01i_c_nested_composition.py — import:
    classify_structure, block_analysis, load_q_samples, load_parents,
    analyze_parent  (added 2026-04-24 chat C; returns the dict the Phase 7
                     wrapper needs, extracted from what used to live in
                     vendored nested_composition_core.py)

Usage (CLI, standalone)
-----------------------
  python3 internal_ancestry_composition.py \
    --q_cache_dir local_Q/ \
    --parents candidates.tsv \
    --outdir nested_composition/ \
    --K 8

Origin: Adapted from MODULE_2B_Step2/helpers/nested_composition.py.
Dedup pass 2026-04-24: this file was byte-identical across four locations;
all duplicates removed or converted to imports. See docs/NESTED_VS_COMPOSITE.md.

=============================================================================
REGISTRY_CONTRACT
  BLOCKS_WRITTEN: none
  # Pure compute engine. Does not touch the registry directly. The Phase 7
  # wrapper STEP_C01i_c_nested_composition.py is the registry writer; see
  # that file's contract for the internal_ancestry_composition block it
  # produces.
  KEYS_IN: none
=============================================================================
"""

import os, sys, csv, gzip, math, argparse
from collections import defaultdict, Counter

def classify_structure(labels):
    """Classify internal Q structure from a sequence of assigned_pop labels."""
    if not labels:
        return "too_sparse"

    unique = set(labels)
    n_unique = len(unique)
    n = len(labels)

    if n < 3:
        return "too_sparse"

    switches = sum(1 for i in range(1, n) if labels[i] != labels[i-1])
    mc = Counter(labels).most_common()
    dom_frac = mc[0][1] / n

    if n_unique == 1:
        return "homogeneous"
    elif dom_frac >= 0.8 and switches <= 2:
        return "dominant_plus_secondary"
    elif n_unique == 2 and 0.3 <= dom_frac <= 0.7:
        return "two_block_composite"
    elif switches >= n * 0.4:
        return "multi_block_fragmented"
    elif switches <= 2 and n_unique <= 3:
        # Check monotonic gradient
        first_half_dom = Counter(labels[:n//2]).most_common(1)[0][0]
        second_half_dom = Counter(labels[n//2:]).most_common(1)[0][0]
        if first_half_dom != second_half_dom:
            return "continuous_gradient"
        return "dominant_plus_secondary"
    else:
        return "diffuse_mixed"


def block_analysis(labels):
    """Decompose label sequence into contiguous blocks."""
    if not labels:
        return [], 0

    blocks = []
    current = labels[0]
    length = 1
    for i in range(1, len(labels)):
        if labels[i] == current:
            length += 1
        else:
            blocks.append((current, length))
            current = labels[i]
            length = 1
    blocks.append((current, length))

    switches = len(blocks) - 1
    return blocks, switches


def load_q_samples(q_cache_dir, chrom):
    """Load per-window × per-sample Q from Engine B cache."""
    for ext in [".local_Q_samples.tsv.gz", ".local_Q_samples.tsv"]:
        path = os.path.join(q_cache_dir, chrom + ext)
        if os.path.isfile(path):
            opener = gzip.open if path.endswith(".gz") else open
            rows = []
            with opener(path, "rt") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    rows.append(row)
            return rows
    return []


def load_parents(path):
    """Load parent intervals (candidates, chromosomes, etc.)."""
    parents = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom = row.get("chrom", row.get("chr", ""))
            start = int(row.get("start_bp", row.get("start", 0)))
            end = int(row.get("end_bp", row.get("end", 0)))
            pid = row.get("candidate_id", row.get("interval_id",
                  row.get("parent_id", f"{chrom}_{start}_{end}")))
            parents.append({"parent_id": pid, "chrom": chrom,
                           "start": start, "end": end})
    return parents


# =============================================================================
# Engine core: analyze one parent interval against its Q windows
# =============================================================================
# Extracted 2026-04-24 (chat C) from what used to be duplicated inline in
# main() here AND in the Phase 7 wrapper. Both callers now use this function.
#
# Returns a dict with:
#   parent_id, chrom, start_bp, end_bp, K_used, n_samples_analyzed,
#   n_windows_in_parent, dominant_structure_type, dominant_structure_count,
#   structure_breakdown, structure_counts, mean_fragmentation,
#   mean_internal_entropy, mean_switches, per_sample (list of dicts)
#
# Numeric fields are rounded floats (not pre-formatted strings). Callers
# that need a specific string format do the formatting themselves.
#
# Deliberately does NOT include composite_flag — that's an inversion-specific
# interpretation layer that belongs in the Phase 7 wrapper, not here.

def analyze_parent_composition(parent: dict, q_rows: list,
                                K: int = 8,
                                min_windows_per_sample: int = 3) -> dict | None:
    """Run composition analysis on one parent interval.

    Args:
        parent: dict with keys parent_id, chrom, start, end
        q_rows: list of per-window × per-sample dicts loaded from
                Engine B cache (output of load_q_samples).
        K: number of ancestry components (metadata only).
        min_windows_per_sample: drop samples with fewer than this many
                windows in the interval.

    Returns:
        A dict summarizing the interval's internal ancestry composition,
        or None if no overlapping windows or no samples meeting the
        min_windows_per_sample threshold.
    """
    pid = parent["parent_id"]
    chrom = parent["chrom"]
    pstart = parent["start"]
    pend = parent["end"]

    # Filter windows overlapping the parent region
    overlap = [r for r in q_rows
               if int(r.get("start_bp", 0)) < pend
               and int(r.get("end_bp", 0)) > pstart]
    if not overlap:
        return None

    # Group by sample
    by_sample = defaultdict(list)
    for row in overlap:
        sid = row.get("sample_id", row.get("sample_idx", ""))
        by_sample[sid].append(row)

    per_sample = []
    for sid, rows in by_sample.items():
        rows.sort(key=lambda r: int(r.get("start_bp", 0)))
        labels = [str(r.get("assigned_pop", "")) for r in rows]
        n_sub = len(rows)
        if n_sub < min_windows_per_sample:
            continue
        blocks, switches = block_analysis(labels)
        label_counts = Counter(labels)
        mc = label_counts.most_common()
        dom_label = mc[0][0]
        dom_frac = mc[0][1] / n_sub
        sec_label = mc[1][0] if len(mc) > 1 else ""
        sec_frac = mc[1][1] / n_sub if len(mc) > 1 else 0.0
        probs = [c / n_sub for c in label_counts.values()]
        H = -sum(p * math.log(p + 1e-15) for p in probs)
        frag = switches / max(1, n_sub - 1)
        struct = classify_structure(labels)
        sorted_blocks = sorted(blocks, key=lambda x: -x[1])
        largest = sorted_blocks[0][1] if sorted_blocks else 0
        second = sorted_blocks[1][1] if len(sorted_blocks) > 1 else 0

        mean_d12 = sum(float(r.get("delta12", 0)) for r in rows) / n_sub
        mean_H = sum(float(r.get("entropy", 0)) for r in rows) / n_sub
        mean_ena = sum(float(r.get("ena", 0)) for r in rows) / n_sub

        per_sample.append({
            "sample_id": sid,
            "n_windows": n_sub,
            "n_labels_unique": len(label_counts),
            "dominant_label": dom_label,
            "dominant_fraction": round(dom_frac, 4),
            "second_label": sec_label,
            "second_fraction": round(sec_frac, 4),
            "n_blocks": len(blocks),
            "n_switches": switches,
            "largest_block": largest,
            "second_block": second,
            "fragmentation_score": round(frag, 4),
            "internal_entropy": round(H, 4),
            "structure_type": struct,
            "mean_delta12": round(mean_d12, 4),
            "mean_entropy": round(mean_H, 4),
            "mean_ena": round(mean_ena, 4),
        })

    if not per_sample:
        return None

    struct_counts = Counter(r["structure_type"] for r in per_sample)
    struct_breakdown = "; ".join(f"{k}:{v}" for k, v in struct_counts.most_common())
    mean_frag = sum(r["fragmentation_score"] for r in per_sample) / len(per_sample)
    mean_int_H = sum(r["internal_entropy"] for r in per_sample) / len(per_sample)
    mean_switches = sum(r["n_switches"] for r in per_sample) / len(per_sample)

    return {
        "parent_id": pid,
        "chrom": chrom,
        "start_bp": pstart,
        "end_bp": pend,
        "K_used": K,
        "n_samples_analyzed": len(per_sample),
        "n_windows_in_parent": len({int(r["start_bp"]) for r in overlap}),
        "dominant_structure_type": struct_counts.most_common(1)[0][0],
        "dominant_structure_count": struct_counts.most_common(1)[0][1],
        "structure_breakdown": struct_breakdown,
        "structure_counts": dict(struct_counts),
        "mean_fragmentation": round(mean_frag, 4),
        "mean_internal_entropy": round(mean_int_H, 4),
        "mean_switches": round(mean_switches, 2),
        "per_sample": per_sample,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q_cache_dir", required=True, help="Engine B cache directory")
    ap.add_argument("--parents", required=True, help="Parent intervals TSV")
    ap.add_argument("--outdir", default="nested_composition/")
    ap.add_argument("--K", type=int, default=8)
    ap.add_argument("--sample_list", default=None, help="Sample IDs file")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    parents = load_parents(args.parents)

    # Load sample names (metadata only — not used downstream in this path)
    if args.sample_list and os.path.isfile(args.sample_list):
        with open(args.sample_list) as f:
            _ = [l.strip() for l in f if l.strip()]

    per_sample_rows = []
    summary_rows = []

    # Cache Q data per chromosome
    q_cache_per_chr = {}

    for parent in parents:
        chrom = parent["chrom"]
        if chrom not in q_cache_per_chr:
            q_cache_per_chr[chrom] = load_q_samples(args.q_cache_dir, chrom)
            if not q_cache_per_chr[chrom]:
                print(f"[SKIP] No Q data for {chrom}")

        q_rows = q_cache_per_chr[chrom]
        if not q_rows:
            continue

        result = analyze_parent_composition(parent, q_rows, K=args.K)
        if result is None:
            continue

        pid = result["parent_id"]
        pstart = result["start_bp"]
        pend = result["end_bp"]

        # CLI TSV rows preserve the legacy string-formatted numeric style
        # for back-compat with downstream readers of nested_composition.tsv.
        for ps in result["per_sample"]:
            per_sample_rows.append({
                "parent_id": pid,
                "chrom": chrom,
                "parent_start": pstart,
                "parent_end": pend,
                "sample_id": ps["sample_id"],
                "n_windows": ps["n_windows"],
                "n_labels_unique": ps["n_labels_unique"],
                "dominant_label": ps["dominant_label"],
                "dominant_fraction": f"{ps['dominant_fraction']:.4f}",
                "second_label": ps["second_label"],
                "second_fraction": f"{ps['second_fraction']:.4f}",
                "n_blocks": ps["n_blocks"],
                "n_switches": ps["n_switches"],
                "largest_block": ps["largest_block"],
                "second_block": ps["second_block"],
                "fragmentation_score": f"{ps['fragmentation_score']:.4f}",
                "internal_entropy": f"{ps['internal_entropy']:.4f}",
                "structure_type": ps["structure_type"],
                "mean_delta12": f"{ps['mean_delta12']:.4f}",
                "mean_entropy": f"{ps['mean_entropy']:.4f}",
                "mean_ena": f"{ps['mean_ena']:.4f}",
            })

        summary_rows.append({
            "parent_id": pid,
            "chrom": chrom,
            "parent_start": pstart,
            "parent_end": pend,
            "n_samples": result["n_samples_analyzed"],
            "dominant_structure_type": result["dominant_structure_type"],
            "dominant_structure_count": result["dominant_structure_count"],
            "mean_fragmentation": f"{result['mean_fragmentation']:.4f}",
            "mean_internal_entropy": f"{result['mean_internal_entropy']:.4f}",
            "mean_switches": f"{result['mean_switches']:.1f}",
            "structure_breakdown": result["structure_breakdown"],
        })

    # Write per-sample
    if per_sample_rows:
        out1 = os.path.join(args.outdir, "nested_composition.tsv")
        fields = list(per_sample_rows[0].keys())
        with open(out1, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(per_sample_rows)
        print(f"[OK] {out1} ({len(per_sample_rows)} rows)")

    # Write summary
    if summary_rows:
        out2 = os.path.join(args.outdir, "nested_composition_summary.tsv")
        fields = list(summary_rows[0].keys())
        with open(out2, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(summary_rows)
        print(f"[OK] {out2} ({len(summary_rows)} parents)")

    print(f"\n[DONE] {len(per_sample_rows)} per-sample rows, {len(summary_rows)} parent summaries")


if __name__ == "__main__":
    main()
