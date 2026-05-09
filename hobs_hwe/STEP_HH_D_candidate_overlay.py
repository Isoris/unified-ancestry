#!/usr/bin/env python3
"""
STEP_HH_D_candidate_overlay.py — Overlay Hobs/HWE windows on inversion candidates.

For each candidate × subset × scale, extracts overlapping windows and computes:
- mean/median F and Hobs within interval
- outlier burden within interval
- comparison to flanking regions
- distortion pattern classification

Output: candidate_hobs_overlay.tsv

Usage:
  python3 04_candidate_overlay.py \
    --candidates candidates.tsv \
    --windows_dir 03_hobs_windows/ \
    --subsets subset_manifest.tsv \
    --scale 50kb \
    --outdir 04_candidate_overlay/

Citation: Hobs/HWE confirmation concept from Claire Mérot
"""

import os, csv, gzip, argparse, math
from collections import defaultdict

def read_windows(path):
    """Read windowed Hobs/F TSV."""
    rows = []
    if not os.path.isfile(path):
        return rows
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)
    return rows


def classify_pattern(mean_F_int, median_F_int, mean_F_flank, frac_hi_F, frac_lo_Hobs):
    """Classify distortion pattern for manuscript interpretation.

    Patterns (from Part 8 of spec):
    1. broad localized genotype distortion (high mean+median F, many outliers)
    2. signal driven by few extreme loci (high mean, normal median)
    3. regional het deficit (low Hobs)
    4-5. need subset comparison (handled at report level)
    """
    if mean_F_int is None or mean_F_flank is None:
        return "insufficient_data"

    F_diff = mean_F_int - mean_F_flank

    if F_diff > 0.05 and frac_hi_F > 0.2:
        if abs(mean_F_int - median_F_int) < 0.02:
            return "broad_localized_distortion"
        else:
            return "few_extreme_loci"
    elif frac_lo_Hobs > 0.3:
        return "regional_het_deficit"
    elif abs(F_diff) < 0.02:
        return "no_signal"
    else:
        return "weak_signal"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--windows_dir", required=True, help="03_hobs_windows/ root")
    ap.add_argument("--subsets", required=True, help="subset_manifest.tsv")
    ap.add_argument("--scale", default="50kb", help="Window scale label")
    ap.add_argument("--flank_bp", type=int, default=500000, help="Flank size for comparison")
    ap.add_argument("--outdir", default="04_candidate_overlay/")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load candidates
    with open(args.candidates) as f:
        candidates = list(csv.DictReader(f, delimiter="\t"))

    # Load subset manifest
    subsets = []
    with open(args.subsets) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            subsets.append(r)

    fields = [
        "candidate_id", "chrom", "start", "end", "subset_id", "scale",
        "n_windows_interval", "n_windows_flank",
        "mean_Hobs_interval", "median_Hobs_interval",
        "mean_F_interval", "median_F_interval",
        "mean_Hobs_flank", "mean_F_flank",
        "frac_high_F_outlier_interval", "frac_low_Hobs_outlier_interval",
        "F_diff_vs_flank", "Hobs_diff_vs_flank",
        "pattern_class", "notes"
    ]

    out_path = os.path.join(args.outdir, f"candidate_hobs_overlay_{args.scale}.tsv")
    with open(out_path, "w", newline="") as fout:
        w = csv.DictWriter(fout, fieldnames=fields, delimiter="\t")
        w.writeheader()

        for cand in candidates:
            c_chr = cand.get("chrom", cand.get("chr", ""))
            c_start = int(cand.get("start_bp", cand.get("start", 0)))
            c_end = int(cand.get("end_bp", cand.get("end", 0)))
            c_id = cand.get("candidate_id", cand.get("interval_id",
                   f"{c_chr}_{c_start}_{c_end}"))

            for sub in subsets:
                sid = sub["subset_id"]

                # Load windows for this subset × chromosome × scale
                win_path = os.path.join(args.windows_dir, sid,
                                        f"{c_chr}.win{args.scale}.tsv")
                windows = read_windows(win_path)
                if not windows:
                    continue

                # Split into interval and flank windows
                interval_wins = []
                flank_wins = []
                flank_start = max(0, c_start - args.flank_bp)
                flank_end = c_end + args.flank_bp

                for win in windows:
                    ws = int(win["window_start"])
                    we = int(win["window_end"])
                    wc = (ws + we) // 2

                    if ws < c_end and we > c_start:
                        # Overlaps interval
                        interval_wins.append(win)
                    elif (ws >= flank_start and we <= c_start) or \
                         (ws >= c_end and we <= flank_end):
                        # In flanks
                        flank_wins.append(win)

                if not interval_wins:
                    continue

                # Compute stats
                def safe_mean(wins, col):
                    vals = [float(w[col]) for w in wins if w.get(col)]
                    return sum(vals)/len(vals) if vals else None

                def safe_median(wins, col):
                    vals = sorted([float(w[col]) for w in wins if w.get(col)])
                    if not vals: return None
                    n = len(vals)
                    return vals[n//2] if n%2==1 else (vals[n//2-1]+vals[n//2])/2

                mH_int = safe_mean(interval_wins, "mean_Hobs")
                mdH_int = safe_median(interval_wins, "median_Hobs")
                mF_int = safe_mean(interval_wins, "mean_F")
                mdF_int = safe_median(interval_wins, "median_F")
                mH_flk = safe_mean(flank_wins, "mean_Hobs")
                mF_flk = safe_mean(flank_wins, "mean_F")

                # Outlier fractions in interval
                frac_hi_F = safe_mean(interval_wins, "frac_high_F_outlier") or 0
                frac_lo_H = safe_mean(interval_wins, "frac_low_Hobs_outlier") or 0

                F_diff = (mF_int - mF_flk) if mF_int is not None and mF_flk is not None else None
                H_diff = (mH_int - mH_flk) if mH_int is not None and mH_flk is not None else None

                pattern = classify_pattern(mF_int, mdF_int, mF_flk, frac_hi_F, frac_lo_H)

                row = {
                    "candidate_id": c_id, "chrom": c_chr,
                    "start": c_start, "end": c_end,
                    "subset_id": sid, "scale": args.scale,
                    "n_windows_interval": len(interval_wins),
                    "n_windows_flank": len(flank_wins),
                    "mean_Hobs_interval": f"{mH_int:.6f}" if mH_int else "",
                    "median_Hobs_interval": f"{mdH_int:.6f}" if mdH_int else "",
                    "mean_F_interval": f"{mF_int:.6f}" if mF_int else "",
                    "median_F_interval": f"{mdF_int:.6f}" if mdF_int else "",
                    "mean_Hobs_flank": f"{mH_flk:.6f}" if mH_flk else "",
                    "mean_F_flank": f"{mF_flk:.6f}" if mF_flk else "",
                    "frac_high_F_outlier_interval": f"{frac_hi_F:.4f}",
                    "frac_low_Hobs_outlier_interval": f"{frac_lo_H:.4f}",
                    "F_diff_vs_flank": f"{F_diff:.6f}" if F_diff is not None else "",
                    "Hobs_diff_vs_flank": f"{H_diff:.6f}" if H_diff is not None else "",
                    "pattern_class": pattern,
                    "notes": "",
                }
                w.writerow(row)

    print(f"[DONE] {out_path}")


if __name__ == "__main__":
    main()
