#!/usr/bin/env python3
"""
build_registries.py — Build interval / sample-subset / set registries.

Generates:
  registries/interval_registry.tsv    — all intervals at all scales
  registries/sample_subsets.tsv       — named sample subsets with paths
  registries/set_registry.tsv         — generic registry of file 'sets'
                                          (cov / sites / mask / ... — discriminated
                                          by the set_type column)

Usage:
  python3 build_registries.py \
    --config 00_ancestry_config.sh \
    --outdir registries/ \
    [--candidates candidates.tsv]

Origin: MODULE_2B_Step2/helpers/make_intervals.py + make_subsets.py + make_cov_registry.py
"""

import os, sys, csv, gzip, glob, argparse

def parse_config(path):
    cfg = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or "=" not in line: continue
            line = line.lstrip("export").strip()
            eq = line.find("=")
            if eq > 0:
                k = line[:eq].strip()
                v = line[eq+1:].strip().strip('"').strip("'")
                # Expand ${VAR} references
                while "${" in v:
                    start = v.find("${")
                    end = v.find("}", start)
                    if end < 0: break
                    ref_key = v[start+2:end]
                    ref_val = cfg.get(ref_key, os.environ.get(ref_key, ""))
                    v = v[:start] + ref_val + v[end+1:]
                cfg[k] = v
    return cfg


def read_fai(path):
    chroms = []
    with open(path) as f:
        for line in f:
            p = line.strip().split("\t")
            chroms.append((p[0], int(p[1])))
    return chroms


def read_beagle_markers(path):
    markers = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        next(f)
        for line in f:
            m = line.split("\t", 1)[0]
            parts = m.rsplit("_", 1)
            if len(parts) == 2:
                try: markers.append((parts[0], int(parts[1])))
                except ValueError: pass
    return markers


def build_intervals(cfg, candidates_path=None):
    fai_path = cfg.get("REF_FAI", cfg.get("FAI", ""))
    beagle_dir = cfg.get("BEAGLE_DIR", cfg.get("STEP1_BEAGLE_DIR", ""))

    chroms = read_fai(fai_path) if os.path.isfile(fai_path) else []
    chrom_sizes = dict(chroms)

    # Read markers from first BEAGLE file for SNP positions
    by_chr = {}
    if beagle_dir and os.path.isdir(beagle_dir):
        for bgl in sorted(glob.glob(os.path.join(beagle_dir, "*.beagle.gz"))):
            for c, p in read_beagle_markers(bgl):
                by_chr.setdefault(c, []).append(p)
        for c in by_chr:
            by_chr[c].sort()

    intervals = []
    iid = 0

    def add(scale, chrom, start, end, mode, wfam="", note=""):
        nonlocal iid
        iid += 1
        nm = sum(1 for p in by_chr.get(chrom, []) if start <= p < end)
        intervals.append({
            "interval_id": f"INT_{iid:06d}",
            "scale_type": scale, "chrom": chrom,
            "start": start, "end": end, "interval_bp": end - start,
            "window_family": wfam, "mode": mode,
            "n_markers": nm, "note": note,
        })

    # Genome
    total_bp = sum(s for _, s in chroms)
    total_mk = sum(len(v) for v in by_chr.values())
    intervals.append({
        "interval_id": "INT_GENOME", "scale_type": "genome",
        "chrom": "genome", "start": 0, "end": total_bp,
        "interval_bp": total_bp, "window_family": "genome",
        "mode": "genome_direct", "n_markers": total_mk, "note": "",
    })

    # Chromosomes
    for c, sz in chroms:
        nm = len(by_chr.get(c, []))
        intervals.append({
            "interval_id": f"INT_CHR_{c}", "scale_type": "chromosome",
            "chrom": c, "start": 0, "end": sz, "interval_bp": sz,
            "window_family": "chromosome", "mode": "chromosome_direct",
            "n_markers": nm, "note": "",
        })

    # SNP windows (100/20 from Engine B default)
    for c in sorted(by_chr.keys()):
        positions = by_chr[c]
        wsz, wst = 100, 20
        i = 0
        while i < len(positions):
            j = min(i + wsz, len(positions))
            add("window", c, positions[i], positions[j-1] + 1,
                "snp_window", f"snp_{wsz}", f"snps={j-i}")
            i += wst
            if j == len(positions): break

    # Candidate intervals
    if candidates_path and os.path.isfile(candidates_path):
        with open(candidates_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                c = row.get("chrom", row.get("chr", ""))
                s = int(row.get("start_bp", row.get("start", 0)))
                e = int(row.get("end_bp", row.get("end", 0)))
                cid = row.get("candidate_id", row.get("interval_id", f"CAND_{c}_{s}_{e}"))
                add("region", c, s, e, "candidate_direct", "candidate", cid)

    return intervals


def build_subsets(cfg):
    subsets = []

    slist = cfg.get("SAMPLE_LIST", cfg.get("STEP1_SAMPLE_LIST", ""))
    if os.path.isfile(slist):
        with open(slist) as f:
            n = sum(1 for l in f if l.strip())
        subsets.append({
            "subset_id": "all", "subset_type": "full_cohort",
            "description": f"All {n} samples", "n_samples": n,
            "sample_list_path": slist,
        })

    pruned = cfg.get("STEP1_PRUNED_LIST", "")
    if os.path.isfile(pruned):
        with open(pruned) as f:
            n = sum(1 for l in f if l.strip())
        subsets.append({
            "subset_id": "pruned", "subset_type": "unrelated_pruned",
            "description": f"Pruned unrelated {n} samples", "n_samples": n,
            "sample_list_path": pruned,
        })

    return subsets


def _detect_chrom_in_path(path):
    """Best-effort chromosome detection from a filesystem path."""
    parts = path.split(os.sep)
    for p in parts:
        if p.startswith("C_gar_LG") or (p.startswith("LG") and p[2:].isdigit()):
            return p
    return ""


def build_set_registry(cfg):
    """Generic file-set registry: cov / sites / mask / ... All discriminated
    by the set_type column. Extend the per-type scanners below as more
    artifact types come online (callable masks, recombination maps, etc.).
    Schema: engines/schemas/set_registry.schema.json."""
    rows = []
    step1 = cfg.get("STEP1_DIR", "")

    # ── set_type='cov' — PCAngsd covariance matrices ─────────────────────
    pcangsd_dir = os.path.join(step1, "05_pcangsd_byLG") if step1 else ""
    if pcangsd_dir and os.path.isdir(pcangsd_dir):
        for cov_path in sorted(glob.glob(os.path.join(pcangsd_dir, "**", "*.cov"), recursive=True)):
            chrom = _detect_chrom_in_path(cov_path)
            rows.append({
                "set_id":      os.path.basename(cov_path).replace(".cov", ""),
                "set_type":    "cov",
                "set_subtype": "pcangsd_default",
                "scope":       "chromosome" if chrom else "global",
                "chrom":       chrom or "genome",
                "path":        cov_path,
                "description": "PCAngsd .cov" + (f" for {chrom}" if chrom else ""),
            })

    # ── set_type='sites' — ANGSD sites lists ─────────────────────────────
    # Single-file SITES_FILE convention, or scan a SITES_DIR.
    sites_file = cfg.get("SITES_FILE", "")
    if sites_file and os.path.isfile(sites_file):
        rows.append({
            "set_id":      os.path.basename(sites_file).replace(".sites", ""),
            "set_type":    "sites",
            "set_subtype": cfg.get("SITES_SUBTYPE", ""),
            "scope":       "global",
            "chrom":       "genome",
            "path":        sites_file,
            "description": "ANGSD sites list from $SITES_FILE",
        })
    sites_dir = cfg.get("SITES_DIR", "")
    if sites_dir and os.path.isdir(sites_dir):
        for sf in sorted(glob.glob(os.path.join(sites_dir, "**", "*.sites"), recursive=True)):
            chrom = _detect_chrom_in_path(sf)
            rows.append({
                "set_id":      os.path.basename(sf).replace(".sites", ""),
                "set_type":    "sites",
                "set_subtype": "",
                "scope":       "chromosome" if chrom else "global",
                "chrom":       chrom or "genome",
                "path":        sf,
                "description": "",
            })

    # ── set_type='callable_bed' — callable-region masks ──────────────────
    callable_bed = cfg.get("CALLABLE_BED", "")
    if callable_bed and os.path.isfile(callable_bed):
        rows.append({
            "set_id":      os.path.basename(callable_bed).replace(".bed", ""),
            "set_type":    "callable_bed",
            "set_subtype": cfg.get("CALLABLE_SUBTYPE", ""),
            "scope":       "global",
            "chrom":       "genome",
            "path":        callable_bed,
            "description": "Callable regions from $CALLABLE_BED",
        })

    # ── set_type='mask' — any other regions-to-mask BED ──────────────────
    mask_bed = cfg.get("MASK_BED", "")
    if mask_bed and os.path.isfile(mask_bed):
        rows.append({
            "set_id":      os.path.basename(mask_bed).replace(".bed", ""),
            "set_type":    "mask",
            "set_subtype": "",
            "scope":       "global",
            "chrom":       "genome",
            "path":        mask_bed,
            "description": "Mask BED from $MASK_BED",
        })

    # ── set_type='recombination_map' — per-genome rate map ───────────────
    rec_map = cfg.get("RECOMBINATION_MAP", "")
    if rec_map and os.path.isfile(rec_map):
        rows.append({
            "set_id":      os.path.basename(rec_map).rsplit(".", 1)[0],
            "set_type":    "recombination_map",
            "set_subtype": "",
            "scope":       "global",
            "chrom":       "genome",
            "path":        rec_map,
            "description": "Recombination rate map from $RECOMBINATION_MAP",
        })

    return rows


# Backward-compat alias — older callers may still import build_cov_registry.
def build_cov_registry(cfg):
    """Deprecated. Use build_set_registry; this returns only cov-typed rows
    in the OLD column shape ({cov_id, scope, chrom, cov_path})."""
    rows = []
    for r in build_set_registry(cfg):
        if r["set_type"] != "cov": continue
        rows.append({"cov_id": r["set_id"], "scope": r["scope"],
                      "chrom": r["chrom"], "cov_path": r["path"]})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default="registries/")
    ap.add_argument("--candidates", default=None)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    cfg = parse_config(args.config)

    # Intervals
    print("[INFO] Building interval registry...")
    intervals = build_intervals(cfg, args.candidates)
    iv_path = os.path.join(args.outdir, "interval_registry.tsv")
    fields = ["interval_id", "scale_type", "chrom", "start", "end",
              "interval_bp", "window_family", "mode", "n_markers", "note"]
    with open(iv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(intervals)

    from collections import Counter
    counts = Counter(i["scale_type"] for i in intervals)
    print(f"[OK] {len(intervals)} intervals: {iv_path}")
    for k, v in sorted(counts.items()):
        print(f"  {k}: {v}")

    # Subsets
    print("\n[INFO] Building sample subset registry...")
    subsets = build_subsets(cfg)
    sub_path = os.path.join(args.outdir, "sample_subsets.tsv")
    with open(sub_path, "w", newline="") as f:
        w = csv.DictWriter(f,
            fieldnames=["subset_id", "subset_type", "description", "n_samples", "sample_list_path"],
            delimiter="\t")
        w.writeheader()
        w.writerows(subsets)
    print(f"[OK] {len(subsets)} subsets: {sub_path}")
    for s in subsets:
        print(f"  {s['subset_id']}: {s['n_samples']} samples")

    # Set registry — generic file sets (cov / sites / mask / ...)
    print("\n[INFO] Building set registry...")
    sets = build_set_registry(cfg)
    set_path = os.path.join(args.outdir, "set_registry.tsv")
    if sets:
        fields = ["set_id", "set_type", "set_subtype", "scope", "chrom", "path", "description"]
        with open(set_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(sets)
        from collections import Counter
        by_type = Counter(r["set_type"] for r in sets)
        print(f"[OK] {len(sets)} sets: {set_path}")
        for k, v in sorted(by_type.items()):
            print(f"  {k}: {v}")
    else:
        print("[INFO] No set files found")


if __name__ == "__main__":
    main()
