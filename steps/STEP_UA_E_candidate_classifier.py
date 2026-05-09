#!/usr/bin/env python3
"""
STEP_UA_E_candidate_classifier.py — Classify inversion candidates from Q support + Engine B.

For each candidate, overlaps Q support windows + Engine B per-sample Q and
produces:
  - coherent_q_loaded_core: high coherence + dominant Q fraction
  - likely_mixed_background: >50% mixed-support windows
  - likely_parasite_windows: high switch count + mixed
  - ancestry_class: "ancestry_driven" vs "structural" vs "mixed"
  - internal_structure_type (from nested_composition)

Usage:
  python3 STEP_UA_E_candidate_classifier.py \
    --candidates candidates.tsv \
    --q_windows snp_q_support/q_support_windows.tsv \
    --q_cache_dir local_Q/ \
    --nested nested_composition/nested_composition_summary.tsv \
    --outdir candidate_classification/ \
    --K 8

Origin: MODULE_2C/helpers/summarize_candidates.py + new logic
"""

import os, csv, argparse, math
from collections import Counter, defaultdict

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--q_windows", default=None, help="q_support_windows.tsv from SNP support")
    ap.add_argument("--q_cache_dir", default=None, help="Engine B local_Q cache")
    ap.add_argument("--nested", default=None, help="nested_composition_summary.tsv")
    ap.add_argument("--outdir", default="candidate_classification/")
    ap.add_argument("--K", type=int, default=8)
    args = ap.parse_args()

    K = args.K
    os.makedirs(args.outdir, exist_ok=True)

    # Load candidates
    with open(args.candidates) as f:
        candidates = list(csv.DictReader(f, delimiter="\t"))

    # Load Q support windows (optional)
    windows = []
    if args.q_windows and os.path.isfile(args.q_windows):
        with open(args.q_windows) as f:
            windows = list(csv.DictReader(f, delimiter="\t"))

    # Load nested composition (optional)
    nested = {}
    if args.nested and os.path.isfile(args.nested):
        with open(args.nested) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                nested[r["parent_id"]] = r

    # Load Engine B summary per chromosome (optional)
    q_summaries = {}
    if args.q_cache_dir:
        import glob, gzip
        for f in glob.glob(os.path.join(args.q_cache_dir, "*.local_Q_summary.tsv*")):
            opener = gzip.open if f.endswith(".gz") else open
            with opener(f, "rt") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    key = (r.get("chrom", ""), int(r.get("start_bp", 0)), int(r.get("end_bp", 0)))
                    q_summaries[key] = r

    mean_q_cols = [f"mean_support_Q{j}" for j in range(1, K+1)]

    fields = ["candidate_id", "chrom", "start", "end",
              "n_windows", "n_snps",
              "dominant_Q", "dominant_Q_fraction",
              "mixed_support_fraction", "support_entropy_mean",
              "q_switch_count", "mean_coherence",
              # Classifications
              "coherent_q_loaded_core", "likely_mixed_background",
              "likely_parasite_windows",
              # Engine B Q metrics (mean across region)
              "localQ_mean_delta12", "localQ_mean_entropy", "localQ_mean_ena",
              # Nested composition
              "nested_dominant_structure", "nested_mean_fragmentation",
              # Final ancestry class
              "ancestry_class",
              "notes"]

    out_path = os.path.join(args.outdir, "candidate_classification.tsv")
    with open(out_path, "w", newline="") as fout:
        w = csv.DictWriter(fout, fieldnames=fields, delimiter="\t")
        w.writeheader()

        for cand in candidates:
            c_chr = cand.get("chrom", cand.get("chr", ""))
            c_start = int(cand.get("start_bp", cand.get("start", 0)))
            c_end = int(cand.get("end_bp", cand.get("end", 0)))
            c_id = cand.get("candidate_id", cand.get("interval_id",
                   f"{c_chr}_{c_start}_{c_end}"))

            # Find overlapping Q support windows
            overlap = [ow for ow in windows
                       if ow["chrom"] == c_chr
                       and int(ow["start"]) < c_end
                       and int(ow["end"]) > c_start]

            n_win = len(overlap)
            n_snps = sum(int(ow["n_snps"]) for ow in overlap) if overlap else 0

            if not overlap:
                w.writerow({"candidate_id": c_id, "chrom": c_chr,
                            "start": c_start, "end": c_end,
                            "n_windows": 0, "n_snps": 0,
                            "ancestry_class": "no_data", "notes": "no overlapping windows"})
                continue

            # Aggregate Q support
            means = []
            for j in range(K):
                col = mean_q_cols[j]
                vals = [float(ow[col]) for ow in overlap if col in ow]
                means.append(sum(vals)/len(vals) if vals else 0)

            dom_idx = max(range(K), key=lambda j: means[j])
            dom_frac = means[dom_idx]

            # Mixed / coherence
            mixed_wins = sum(1 for ow in overlap if ow.get("mixed_support_yes_no") == "yes")
            mixed_frac = mixed_wins / n_win

            ent_vals = [float(ow["support_entropy_mean"]) for ow in overlap]
            ent_mean = sum(ent_vals) / len(ent_vals)

            coh_vals = [float(ow["coherence_score"]) for ow in overlap]
            mean_coh = sum(coh_vals) / len(coh_vals)

            dom_qs = [ow["dominant_Q"] for ow in overlap]
            switches = sum(1 for i in range(1, len(dom_qs)) if dom_qs[i] != dom_qs[i-1])

            # Engine B Q metrics
            lq_d12 = lq_H = lq_ena = ""
            # Scan q_summaries for overlapping windows
            matching_q = []
            for (qchr, qs, qe), qrow in q_summaries.items():
                if qchr == c_chr and qs < c_end and qe > c_start:
                    matching_q.append(qrow)
            if matching_q:
                lq_d12 = f"{sum(float(r.get('mean_delta12',0)) for r in matching_q)/len(matching_q):.4f}"
                lq_H = f"{sum(float(r.get('mean_entropy',0)) for r in matching_q)/len(matching_q):.4f}"
                lq_ena = f"{sum(float(r.get('mean_ena',0)) for r in matching_q)/len(matching_q):.4f}"

            # Nested composition
            nest = nested.get(c_id, {})
            nest_struct = nest.get("dominant_structure_type", "")
            nest_frag = nest.get("mean_fragmentation", "")

            # ── Classifications ──
            coherent = "yes" if mean_coh >= 0.70 and dom_frac >= 0.40 else "no"
            mixed_bg = "yes" if mixed_frac > 0.50 else "no"
            parasite = "yes" if switches > n_win * 0.3 and mixed_frac > 0.3 else "no"

            # Ancestry class: is this region driven by ancestry or structural?
            # High family_likeness + high Q coherence → ancestry-driven
            # Low family_likeness + high structure signal → structural
            if lq_d12:
                d12_val = float(lq_d12)
                if d12_val > 0.6 and dom_frac > 0.5:
                    anc_class = "ancestry_driven"
                elif d12_val < 0.3 and mixed_frac > 0.5:
                    anc_class = "mixed_uncertain"
                else:
                    anc_class = "structural_candidate"
            else:
                if coherent == "yes":
                    anc_class = "likely_ancestry"
                else:
                    anc_class = "unclassified"

            w.writerow({
                "candidate_id": c_id, "chrom": c_chr,
                "start": c_start, "end": c_end,
                "n_windows": n_win, "n_snps": n_snps,
                "dominant_Q": f"Q{dom_idx+1}",
                "dominant_Q_fraction": f"{dom_frac:.4f}",
                "mixed_support_fraction": f"{mixed_frac:.4f}",
                "support_entropy_mean": f"{ent_mean:.4f}",
                "q_switch_count": switches,
                "mean_coherence": f"{mean_coh:.4f}",
                "coherent_q_loaded_core": coherent,
                "likely_mixed_background": mixed_bg,
                "likely_parasite_windows": parasite,
                "localQ_mean_delta12": lq_d12,
                "localQ_mean_entropy": lq_H,
                "localQ_mean_ena": lq_ena,
                "nested_dominant_structure": nest_struct,
                "nested_mean_fragmentation": nest_frag,
                "ancestry_class": anc_class,
                "notes": "",
            })

    print(f"[DONE] {len(candidates)} candidates → {out_path}")

if __name__ == "__main__":
    main()
