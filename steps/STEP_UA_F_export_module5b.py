#!/usr/bin/env python3
"""
STEP_UA_F_export_module5b.py — Export marker/window annotation strips for Module 5B.

Reads:
  - snp_q_support/q_support_per_snp.tsv
  - snp_q_support/q_support_windows.tsv
  - candidate_classification/candidate_classification.tsv (optional)

Writes to exports_for_module5b/:
  1) marker_strip_dominant_Q.tsv     — per-SNP dominant ancestry + color
  2) marker_strip_entropy.tsv        — per-SNP support entropy
  3) marker_strip_confidence.tsv     — per-SNP class + confidence
  4) window_strip_coherence.tsv      — per-window coherence + dominant Q
  5) module5b_integration_table.tsv  — full per-SNP table for Module 5B

Module 5B consumes these as annotation strips alongside heatmaps, regime
tracks, and candidate segmentation views.

Origin: MODULE_2C/helpers/export_for_5b.py
"""

import os, csv, argparse

PALETTE = {f"Q{i+1}": c for i, c in enumerate([
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#86BCB6", "#D37295"
])}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snp_support", required=True, help="q_support_per_snp.tsv")
    ap.add_argument("--window_support", default=None, help="q_support_windows.tsv")
    ap.add_argument("--outdir", default="exports_for_module5b/")
    ap.add_argument("--K", type=int, default=8)
    args = ap.parse_args()

    K = args.K
    os.makedirs(args.outdir, exist_ok=True)

    # Load SNP support
    with open(args.snp_support) as f:
        snps = list(csv.DictReader(f, delimiter="\t"))
    print(f"[INFO] {len(snps)} SNPs loaded")

    # 1) Dominant Q strip
    p = os.path.join(args.outdir, "marker_strip_dominant_Q.tsv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["snp_id", "chrom", "pos", "major_allele", "minor_allele",
                     "best_Q", "best_Q_support", "color_hex",
                     "in_majmin", "in_thin500"])
        for r in snps:
            bq = r["best_Q"]
            w.writerow([r["snp_id"], r["chrom"], r["pos"],
                        r.get("major_allele", ""), r.get("minor_allele", ""),
                        bq, r["best_Q_support"], PALETTE.get(bq, "#999999"),
                        r.get("in_majmin", ""), r.get("in_thin500", "")])
    print(f"[OK] {p}")

    # 2) Entropy strip
    p = os.path.join(args.outdir, "marker_strip_entropy.tsv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["snp_id", "chrom", "pos", "support_entropy", "effective_num_Q"])
        for r in snps:
            w.writerow([r["snp_id"], r["chrom"], r["pos"],
                        r["support_entropy"], r["effective_num_Q"]])
    print(f"[OK] {p}")

    # 3) Confidence / class strip
    p = os.path.join(args.outdir, "marker_strip_confidence.tsv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["snp_id", "chrom", "pos", "marker_class",
                     "support_confidence", "delta12", "delta13"])
        for r in snps:
            w.writerow([r["snp_id"], r["chrom"], r["pos"],
                        r["marker_class"], r["support_confidence"],
                        r["delta12"], r.get("delta13", "")])
    print(f"[OK] {p}")

    # 4) Window coherence strip
    if args.window_support and os.path.isfile(args.window_support):
        p = os.path.join(args.outdir, "window_strip_coherence.tsv")
        with open(args.window_support) as fin, open(p, "w", newline="") as fout:
            reader = csv.DictReader(fin, delimiter="\t")
            w = csv.writer(fout, delimiter="\t")
            w.writerow(["window_id", "chrom", "start", "end", "coherence_score",
                        "dominant_Q", "mixed_support_yes_no", "bestQ_switch_count"])
            for r in reader:
                w.writerow([r["window_id"], r["chrom"], r["start"], r["end"],
                            r["coherence_score"], r["dominant_Q"],
                            r["mixed_support_yes_no"], r["bestQ_switch_count"]])
        print(f"[OK] {p}")

    # 5) Full integration table
    q_cols = [f"support_Q{j}" for j in range(1, K+1)]
    raw_cols = [f"raw_score_Q{j}" for j in range(1, K+1)]
    int_fields = (["snp_id", "chrom", "pos", "major_allele", "minor_allele", "K"]
                  + q_cols + raw_cols
                  + ["best_Q", "best_Q_support", "second_Q", "second_Q_support",
                     "delta12", "delta13", "support_entropy", "effective_num_Q",
                     "dosage_mean", "dosage_sd", "maf_approx",
                     "marker_class", "support_confidence",
                     "in_majmin", "in_thin500"])

    p = os.path.join(args.outdir, "module5b_integration_table.tsv")
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=int_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in snps:
            r["K"] = K
            w.writerow(r)
    print(f"[OK] {p} ({len(snps)} SNPs)")

    print(f"\n[DONE] Module 5B exports: {args.outdir}/")

if __name__ == "__main__":
    main()
