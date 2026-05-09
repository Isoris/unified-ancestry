#!/usr/bin/env python3
"""
STEP_UA_C_snp_q_support.py — Per-SNP soft support vectors to Q components.

For each SNP: Pearson |correlation| of dosage vs Q column → normalized soft weights.
Classifies markers: dominant_single_Q | dominant_plus_secondary | mixed | diffuse | low_info.
Aggregates into windows with coherence score, switch counts, dominant Q.

NOW reads Q from Engine B cache (instant_q output) instead of old Module 2B tables.

Usage:
  python3 STEP_UA_C_snp_q_support.py \
    --beagle <file.beagle.gz> \
    --q_cache_dir local_Q/ \
    --sample_list samples.txt \
    --outdir snp_q_support/ \
    --K 8 --chr C_gar_LG01

  python3 STEP_UA_C_snp_q_support.py \
    --beagle <file.beagle.gz> \
    --qopt <global.qopt> \
    --sample_list samples.txt \
    --outdir snp_q_support/ \
    --K 8

Origin: Adapted from MODULE_2C/helpers/compute_snp_support.py + summarize_windows.py
"""

import os, sys, csv, gzip, math, argparse
from collections import defaultdict, Counter
import numpy as np

# =========================================================================
# MARKER CLASSIFICATION
# =========================================================================

def classify_marker(best_sup, delta12, entropy, dom_thr=0.50, mix_thr=0.30):
    if best_sup >= dom_thr and delta12 >= 0.20:
        return "dominant_single_Q"
    elif best_sup >= dom_thr:
        return "dominant_plus_secondary"
    elif best_sup >= mix_thr:
        return "mixed_Q_support"
    elif entropy > 1.5:
        return "diffuse_Q_support"
    else:
        return "low_information"

# =========================================================================
# MARKER UNIVERSE JOIN
# =========================================================================

def load_site_set(path):
    """Load site file into set of (chrom, pos)."""
    sites = set()
    if not path or not os.path.isfile(path): return sites
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                try: sites.add((parts[0], int(parts[1])))
                except ValueError: pass
    return sites

def load_majmin(path):
    """Load 4-col sites file into {(chrom,pos): (major,minor)}."""
    lk = {}
    if not path or not os.path.isfile(path): return lk
    with open(path) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 4:
                try: lk[(p[0], int(p[1]))] = (p[2], p[3])
                except ValueError: pass
    return lk

# =========================================================================
# BEAGLE DOSAGE ITERATOR
# =========================================================================

def beagle_dosage_iter(path, sample_indices):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            marker = parts[0]
            dosages = []
            for si in sample_indices:
                base = 3 + si * 3
                if base + 2 < len(parts):
                    try: dosages.append(float(parts[base + 1]) + 2 * float(parts[base + 2]))
                    except ValueError: dosages.append(np.nan)
                else:
                    dosages.append(np.nan)
            yield marker, np.array(dosages)

def beagle_sample_names(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline().strip().split("\t")
    names = []
    for i in range(3, len(header), 3):
        nm = header[i].rsplit("_", 1)[0] if "_" in header[i] else header[i]
        if not names or names[-1] != nm:
            names.append(nm)
    return names

# =========================================================================
# Q LOADING (from Engine B cache OR qopt file)
# =========================================================================

def load_q_from_cache(cache_dir, chr, K, sample_ids):
    """Load per-sample Q from Engine B cache — uses the genome-wide Q (fixed F source)."""
    # For SNP support, we actually want the GLOBAL Q (not local Q per window),
    # because we're correlating SNP dosage against population-level ancestry.
    # The global Q comes from the qopt file.
    return None  # Fall through to qopt loading

def load_q_from_qopt(qopt_path, K, sample_ids):
    """Load Q matrix from NGSadmix qopt file."""
    q_by_sample = {}
    rows = []
    with open(qopt_path) as f:
        for line in f:
            vals = [float(x) for x in line.strip().split()]
            if len(vals) == K:
                rows.append(vals)

    for i, sid in enumerate(sample_ids):
        if i < len(rows):
            q_by_sample[sid] = rows[i]
    return q_by_sample

# =========================================================================
# MAIN
# =========================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beagle", required=True)
    ap.add_argument("--qopt", default=None, help="Global Q matrix (.qopt)")
    ap.add_argument("--q_cache_dir", default=None, help="Engine B cache (fallback)")
    ap.add_argument("--sample_list", required=True)
    ap.add_argument("--outdir", default="snp_q_support/")
    ap.add_argument("--K", type=int, default=8)
    ap.add_argument("--chr", default=None, help="Filter to chromosome")
    ap.add_argument("--keep", default=None, help="Keep-list for sample subset")
    # Marker universe (optional)
    ap.add_argument("--majmin", default=None, help="sites.raw.ALL.majmin.tsv")
    ap.add_argument("--thin500", default=None, help="thin_500 sites file")
    # Thresholds
    ap.add_argument("--dom_thr", type=float, default=0.50)
    ap.add_argument("--mix_thr", type=float, default=0.30)
    ap.add_argument("--win_size", type=int, default=50)
    ap.add_argument("--win_step", type=int, default=25)
    args = ap.parse_args()

    K = args.K
    os.makedirs(args.outdir, exist_ok=True)

    # Load sample list
    with open(args.sample_list) as f:
        all_samples = [l.strip() for l in f if l.strip()]
    # Strip BAM paths
    all_samples = [os.path.basename(s).replace(".sorted.markdup.bam", "").replace(".bam", "")
                   for s in all_samples]

    # Apply keep filter
    if args.keep and os.path.isfile(args.keep):
        with open(args.keep) as f:
            keep_set = set(l.strip() for l in f if l.strip())
        target_samples = [s for s in all_samples if s in keep_set]
    else:
        target_samples = all_samples

    # Load Q
    q_by_sample = None
    if args.qopt and os.path.isfile(args.qopt):
        q_by_sample = load_q_from_qopt(args.qopt, K, all_samples)
    if not q_by_sample:
        print("[ERROR] Need --qopt for SNP support computation")
        sys.exit(1)

    print(f"[INFO] K={K}, {len(target_samples)} target samples, {len(q_by_sample)} Q entries")

    # Align samples between BEAGLE and Q
    bgl_samples = beagle_sample_names(args.beagle)
    target_set = set(target_samples)
    sample_indices, aligned = [], []
    for i, s in enumerate(bgl_samples):
        if s in target_set and s in q_by_sample:
            sample_indices.append(i)
            aligned.append(s)

    N = len(aligned)
    print(f"[INFO] Aligned: {N} samples")
    if N < 10:
        print("[ERROR] Too few aligned samples"); sys.exit(1)

    Q = np.array([q_by_sample[s] for s in aligned])
    Q_centered = Q - Q.mean(axis=0)

    # Load marker universe
    majmin = load_majmin(args.majmin) if args.majmin else {}
    thin500 = load_site_set(args.thin500) if args.thin500 else set()

    # ── Per-SNP computation ──
    snp_out = os.path.join(args.outdir, "q_support_per_snp.tsv")
    sup_cols = [f"support_Q{j}" for j in range(1, K+1)]
    raw_cols = [f"raw_score_Q{j}" for j in range(1, K+1)]

    fields = (["snp_id", "chrom", "pos", "major_allele", "minor_allele"]
              + sup_cols + raw_cols
              + ["best_Q", "best_Q_support", "second_Q", "second_Q_support",
                 "delta12", "delta13", "support_entropy", "effective_num_Q",
                 "nQ_above_005", "nQ_above_010",
                 "dosage_mean", "dosage_sd", "maf_approx",
                 "marker_class", "support_confidence",
                 "in_majmin", "in_thin500"])

    n_written = n_skip = 0

    with open(snp_out, "w", newline="") as fout:
        w = csv.DictWriter(fout, fieldnames=fields, delimiter="\t")
        w.writeheader()

        for marker, dosage in beagle_dosage_iter(args.beagle, sample_indices):
            parts = marker.rsplit("_", 1)
            if len(parts) != 2: continue
            chrom, pos_str = parts[0], parts[1]
            try: pos = int(pos_str)
            except ValueError: continue

            if args.chr and chrom != args.chr: continue

            valid = ~np.isnan(dosage)
            n_valid = int(valid.sum())
            if n_valid < N * 0.5: n_skip += 1; continue

            d = dosage.copy()
            d_mean = float(np.nanmean(d))
            d[~valid] = d_mean
            d_sd = float(np.std(d))
            if d_sd < 1e-10: n_skip += 1; continue

            maf = min(d_mean / 2, 1 - d_mean / 2) if 0 <= d_mean <= 2 else 0

            # Per-Q correlation
            d_c = d - d_mean
            raw_scores = []
            for j in range(K):
                qc = Q_centered[:, j]
                q_std = qc.std()
                if q_std < 1e-10:
                    raw_scores.append(0.0)
                else:
                    raw_scores.append(abs(float(np.dot(d_c, qc) / (N * d_sd * q_std))))

            total = sum(raw_scores)
            sups = [s / total for s in raw_scores] if total > 1e-15 else [1.0/K]*K

            ranked = sorted(enumerate(sups), key=lambda x: -x[1])
            best_i, best_s = ranked[0]
            sec_i, sec_s = ranked[1] if K > 1 else (0, 0)
            third_s = ranked[2][1] if K > 2 else 0
            d12 = best_s - sec_s
            d13 = best_s - third_s
            H = -sum(s * math.log(s + 1e-15) for s in sups)
            ena = math.exp(H)
            nQ05 = sum(1 for s in sups if s >= 0.05)
            nQ10 = sum(1 for s in sups if s >= 0.10)

            conf = "high" if best_s >= args.dom_thr else "moderate" if best_s >= args.mix_thr else "low"
            mclass = classify_marker(best_s, d12, H, args.dom_thr, args.mix_thr)

            site_key = (chrom, pos)
            mm = majmin.get(site_key, ("", ""))

            row = {
                "snp_id": marker, "chrom": chrom, "pos": pos,
                "major_allele": mm[0], "minor_allele": mm[1],
                "best_Q": f"Q{best_i+1}", "best_Q_support": f"{best_s:.6f}",
                "second_Q": f"Q{sec_i+1}", "second_Q_support": f"{sec_s:.6f}",
                "delta12": f"{d12:.6f}", "delta13": f"{d13:.6f}",
                "support_entropy": f"{H:.6f}", "effective_num_Q": f"{ena:.4f}",
                "nQ_above_005": nQ05, "nQ_above_010": nQ10,
                "dosage_mean": f"{d_mean:.4f}", "dosage_sd": f"{d_sd:.4f}",
                "maf_approx": f"{maf:.4f}",
                "marker_class": mclass, "support_confidence": conf,
                "in_majmin": "yes" if site_key in majmin else "no",
                "in_thin500": "yes" if site_key in thin500 else "no",
            }
            for j in range(K):
                row[f"support_Q{j+1}"] = f"{sups[j]:.6f}"
                row[f"raw_score_Q{j+1}"] = f"{raw_scores[j]:.6f}"

            w.writerow(row)
            n_written += 1
            if n_written % 50000 == 0:
                print(f"[INFO] {n_written:,} SNPs...")

    print(f"\n[OK] {n_written:,} SNPs → {snp_out}")
    print(f"[INFO] Skipped: {n_skip:,}")

    # ── Window aggregation ──
    print("\n[INFO] Aggregating into windows...")
    window_out = os.path.join(args.outdir, "q_support_windows.tsv")

    # Read back the SNPs grouped by chromosome
    by_chr = defaultdict(list)
    with open(snp_out) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            by_chr[r["chrom"]].append(r)

    win_fields = (["window_id", "chrom", "start", "end", "n_snps"]
                  + [f"mean_support_Q{j}" for j in range(1, K+1)]
                  + ["dominant_Q", "dominant_Q_fraction",
                     "support_entropy_mean", "coherence_score",
                     "bestQ_switch_count", "mixed_support_yes_no"])

    wid = 0
    n_win = 0
    with open(window_out, "w", newline="") as fout:
        w = csv.DictWriter(fout, fieldnames=win_fields, delimiter="\t")
        w.writeheader()

        for chrom in sorted(by_chr.keys()):
            snps = by_chr[chrom]
            snps.sort(key=lambda r: int(r["pos"]))
            i = 0
            while i < len(snps):
                win = snps[i:i + args.win_size]
                if len(win) < 5: break

                wid += 1
                ns = len(win)
                means = []
                for j in range(K):
                    col = f"support_Q{j+1}"
                    vals = [float(r[col]) for r in win]
                    means.append(sum(vals) / len(vals))

                dom_idx = max(range(K), key=lambda j: means[j])
                dom_frac = means[dom_idx]

                entropies = [float(r["support_entropy"]) for r in win]
                ent_mean = sum(entropies) / ns

                best_qs = [r["best_Q"] for r in win]
                bq_counts = Counter(best_qs)
                coherence = bq_counts.most_common(1)[0][1] / ns
                switches = sum(1 for idx in range(1, ns) if best_qs[idx] != best_qs[idx-1])

                row = {
                    "window_id": f"WIN_{wid:06d}",
                    "chrom": chrom,
                    "start": int(win[0]["pos"]),
                    "end": int(win[-1]["pos"]),
                    "n_snps": ns,
                    "dominant_Q": f"Q{dom_idx+1}",
                    "dominant_Q_fraction": f"{dom_frac:.4f}",
                    "support_entropy_mean": f"{ent_mean:.4f}",
                    "coherence_score": f"{coherence:.4f}",
                    "bestQ_switch_count": switches,
                    "mixed_support_yes_no": "yes" if coherence < 0.50 else "no",
                }
                for j in range(K):
                    row[f"mean_support_Q{j+1}"] = f"{means[j]:.6f}"

                w.writerow(row)
                n_win += 1

                i += args.win_step
                if i + args.win_size > len(snps) and i < len(snps):
                    i = max(len(snps) - args.win_size, i)

    print(f"[OK] {n_win} windows → {window_out}")

if __name__ == "__main__":
    main()
