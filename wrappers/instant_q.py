#!/usr/bin/env python3
"""
instant_q.py — Python wrapper for the Instant Q Engine.

Provides:
    get_Q(chr, start, end)        — single region query
    get_Q_summary(chr)            — per-window summary for a chromosome
    get_region_stats(chr, start, end, what, groups)  — unified stats dispatcher

Usage:
    from instant_q import configure, get_Q, get_Q_summary, get_region_stats

    configure(config_file="00_ancestry_config.sh")
    q = get_Q("C_gar_LG01", start=35e6, end=41e6)
    summary = get_Q_summary("C_gar_LG01")
    stats = get_region_stats("C_gar_LG01", 35e6, 41e6,
                              what=["Q", "Fst", "Hobs"],
                              groups={"inv": inv_samples, "ref": ref_samples})
"""

import os
import sys
import gzip
import subprocess
import tempfile
import shutil
import glob
import math
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Union

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

EPS_Q = 1e-10

# =============================================================================
# CONFIGURATION
# =============================================================================

_config = {}

def configure(
    config_file: str = None,
    binary: str = None,
    beagle_dir: str = None,
    fopt_file: str = None,
    qinit_file: str = None,
    cache_dir: str = None,
    sample_list: str = None,
    window_size: int = 100,
    window_step: int = 20,
    em_iter: int = 20,
    K: int = None,
    ncores: int = 4
):
    """Configure the Instant Q engine."""
    global _config

    if config_file and os.path.isfile(config_file):
        cfg = _parse_shell_config(config_file)
        beagle_dir  = beagle_dir  or cfg.get("BEAGLE_DIR") or cfg.get("STEP1_BEAGLE_DIR")
        fopt_file   = fopt_file   or cfg.get("BEST_FOPT")
        qinit_file  = qinit_file  or cfg.get("BEST_QOPT")
        cache_dir   = cache_dir   or cfg.get("LOCAL_Q_DIR")
        sample_list = sample_list or cfg.get("SAMPLE_LIST") or cfg.get("STEP1_SAMPLE_LIST")
        binary      = binary      or cfg.get("INSTANT_Q_BIN")

    _config["binary"]      = binary or _find_binary()
    _config["beagle_dir"]  = beagle_dir
    _config["fopt_file"]   = fopt_file
    _config["qinit_file"]  = qinit_file
    _config["cache_dir"]   = cache_dir or "local_Q"
    _config["sample_list"] = sample_list
    _config["window_size"] = window_size
    _config["window_step"] = window_step
    _config["em_iter"]     = em_iter
    _config["ncores"]      = ncores

    # Load sample names
    if sample_list and os.path.isfile(sample_list):
        with open(sample_list) as f:
            _config["samples"] = [l.strip() for l in f if l.strip()]

    # Auto-detect K
    if qinit_file and os.path.isfile(qinit_file):
        with open(qinit_file) as f:
            first = f.readline().strip().split()
            _config["K"] = len(first)
    if K is not None:
        _config["K"] = K

    os.makedirs(_config["cache_dir"], exist_ok=True)

    print(f"[instant_q.py] Configured: K={_config.get('K', 'auto')}, "
          f"cache={_config['cache_dir']}, "
          f"binary={_config.get('binary', 'NOT FOUND')}",
          file=sys.stderr)


def _parse_shell_config(path: str) -> dict:
    cfg = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or "=" not in line:
                continue
            line = line.lstrip("export").strip()
            eq = line.find("=")
            if eq > 0:
                key = line[:eq].strip()
                val = line[eq+1:].strip().strip('"').strip("'")
                cfg[key] = val
    return cfg


def _find_binary() -> Optional[str]:
    candidates = [
        os.environ.get("INSTANT_Q_BIN", ""),
        os.path.join(os.path.dirname(__file__), "src", "instant_q"),
        os.path.join(os.path.dirname(__file__), "instant_q"),
        os.path.expanduser("~/.local/bin/instant_q"),
    ]
    for p in candidates:
        if p and os.path.isfile(p) and os.access(p, os.X_OK):
            return os.path.abspath(p)
    w = shutil.which("instant_q")
    return w


def _resolve_beagle(chrom: str) -> str:
    bd = _config.get("beagle_dir")
    if not bd:
        raise ValueError("[instant_q.py] beagle_dir not configured")
    patterns = [
        f"*{chrom}*.beagle.gz",
        f"catfish.{chrom}.rf.txt.thin_500.beagle.gz",
        f"catfish.{chrom}*.beagle.gz",
    ]
    for pat in patterns:
        matches = glob.glob(os.path.join(bd, pat))
        if matches:
            return matches[0]
    raise FileNotFoundError(f"[instant_q.py] No BEAGLE for {chrom} in {bd}")


def _cache_summary(chrom: str) -> str:
    return os.path.join(_config["cache_dir"], f"{chrom}.local_Q_summary.tsv.gz")


def _cache_exists(chrom: str) -> bool:
    p = _cache_summary(chrom)
    return os.path.isfile(p) or os.path.isfile(p.replace(".gz", ""))


# =============================================================================
# get_Q() — Single region
# =============================================================================

def get_Q(chrom: str, start: int, end: int, force: bool = False):
    """
    Get Q vectors for a single genomic region.

    Returns: pandas DataFrame if available, else list of dicts.
    """
    binary = _config.get("binary")
    if not binary or not os.path.isfile(binary):
        raise RuntimeError("[instant_q.py] C++ binary not found. Compile first.")

    beagle = _resolve_beagle(chrom)

    cmd = [
        binary,
        "--beagle", beagle,
        "--fopt", _config["fopt_file"],
        "--qinit", _config["qinit_file"],
        "--chr", chrom,
        "--start", str(int(start)),
        "--end", str(int(end)),
        "--em_iter", str(_config["em_iter"]),
    ]

    print(f"[instant_q.py] {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd, capture_output=True, text=True)

    if proc.returncode != 0:
        print(proc.stderr, file=sys.stderr)
        raise RuntimeError(f"instant_q returned {proc.returncode}")

    if HAS_PANDAS:
        return pd.read_csv(StringIO(proc.stdout), sep="\t")
    else:
        # Parse TSV manually
        lines = proc.stdout.strip().split("\n")
        if not lines:
            return []
        header = lines[0].split("\t")
        return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


# =============================================================================
# get_Q_summary() — Per-window summary
# =============================================================================

def get_Q_summary(chrom: str, force: bool = False):
    """
    Get per-window Q summary for a chromosome (from cache).
    Triggers precompute if cache is missing.
    """
    cf = _cache_summary(chrom)
    cf_plain = cf.replace(".gz", "")

    for path in [cf, cf_plain]:
        if os.path.isfile(path) and not force:
            if HAS_PANDAS:
                return pd.read_csv(path, sep="\t")
            else:
                return _read_tsv(path)

    # Precompute
    print(f"[instant_q.py] No cache for {chrom} — precomputing", file=sys.stderr)
    precompute(chrom, force=force)

    for path in [cf, cf_plain]:
        if os.path.isfile(path):
            if HAS_PANDAS:
                return pd.read_csv(path, sep="\t")
            else:
                return _read_tsv(path)

    raise RuntimeError(f"Precompute failed for {chrom}")


# =============================================================================
# precompute() — Run precompute for a chromosome
# =============================================================================

def precompute(chrom: str, sample_output: bool = False, force: bool = False):
    """Precompute all windows for a chromosome."""
    binary = _config.get("binary")
    if not binary or not os.path.isfile(binary):
        raise RuntimeError("[instant_q.py] C++ binary not found. Compile first.")

    if _cache_exists(chrom) and not force:
        print(f"[instant_q.py] Cache exists for {chrom}", file=sys.stderr)
        return

    beagle = _resolve_beagle(chrom)
    tmpdir = tempfile.mkdtemp(prefix="iq_")

    cmd = [
        binary,
        "--beagle", beagle,
        "--fopt", _config["fopt_file"],
        "--qinit", _config["qinit_file"],
        "--precompute",
        "--outdir", tmpdir,
        "--chr", chrom,
        "--window_size", str(_config["window_size"]),
        "--window_step", str(_config["window_step"]),
        "--em_iter", str(_config["em_iter"]),
        "--ncores", str(_config["ncores"]),
    ]
    if sample_output:
        cmd.append("--sample_output")

    print(f"[instant_q.py] {' '.join(cmd)}", file=sys.stderr)
    ret = subprocess.run(cmd)

    if ret.returncode != 0:
        raise RuntimeError(f"instant_q returned {ret.returncode}")

    # Move results to cache
    os.makedirs(_config["cache_dir"], exist_ok=True)

    summary_raw = os.path.join(tmpdir, f"{chrom}.local_Q_summary.tsv")
    if os.path.isfile(summary_raw):
        # Gzip and move
        target = _cache_summary(chrom)
        if HAS_PANDAS:
            df = pd.read_csv(summary_raw, sep="\t")
            df.to_csv(target, sep="\t", index=False)
        else:
            shutil.copy2(summary_raw, target.replace(".gz", ""))

    meta_raw = os.path.join(tmpdir, f"{chrom}.local_Q_meta.tsv")
    if os.path.isfile(meta_raw):
        shutil.copy2(meta_raw,
                     os.path.join(_config["cache_dir"], f"{chrom}.local_Q_meta.tsv"))

    shutil.rmtree(tmpdir, ignore_errors=True)
    print(f"[instant_q.py] Cached {chrom}", file=sys.stderr)


# =============================================================================
# get_region_stats() — Unified stats dispatcher
# =============================================================================

def get_region_stats(
    chrom: str,
    start: int,
    end: int,
    what: Union[str, List[str]] = "Q",
    groups: Optional[Dict[str, List[str]]] = None
) -> dict:
    """
    Unified stats dispatcher for any genomic region.

    Args:
        chrom: Chromosome name
        start: Start position (bp)
        end: End position (bp)
        what: Stat(s) to compute. Options: Q, Fst, Hobs, dXY, theta_pi, HWE,
              dosage_matrix
        groups: Named dict of sample lists for pairwise stats.

    Returns:
        Dict with requested stats.
    """
    if not HAS_NUMPY:
        raise ImportError("get_region_stats requires numpy")

    if isinstance(what, str):
        what = [what]

    result = {}

    # Q
    if "Q" in what:
        result["Q"] = get_Q(chrom, start, end)

    # Load dosage for dosage-based stats
    need_dosage = any(s in what for s in
                      ["Fst", "Hobs", "dXY", "theta_pi", "HWE", "dosage_matrix"])
    dos = None
    sample_ids = None

    if need_dosage:
        dos, sample_ids, site_pos = _load_dosage(chrom, start, end)
        if dos is not None and "dosage_matrix" in what:
            result["dosage_matrix"] = dos
            result["snp_positions"] = site_pos

    if dos is None and need_dosage:
        print(f"[instant_q.py] No dosage data for {chrom}:{start}-{end}", file=sys.stderr)
        return result

    n_sites, n_ind = dos.shape if dos is not None else (0, 0)

    # Hobs
    if "Hobs" in what and dos is not None:
        het = (dos >= 0.3) & (dos <= 1.7)
        per_sample = np.nanmean(het, axis=0)
        result["Hobs"] = {
            "per_sample": dict(zip(sample_ids, np.round(per_sample, 4))),
            "mean": round(float(np.nanmean(per_sample)), 4),
            "sd": round(float(np.nanstd(per_sample)), 4),
        }

    # theta_pi
    if "theta_pi" in what and dos is not None:
        p_hat = np.nanmean(dos, axis=1) / 2.0
        pi_site = 2 * p_hat * (1 - p_hat) * n_ind / max(n_ind - 1, 1)
        result["theta_pi"] = round(float(np.nanmean(pi_site)), 6)

    # HWE
    if "HWE" in what and dos is not None:
        p_hat = np.nanmean(dos, axis=1) / 2.0
        obs_het = np.nanmean((dos >= 0.3) & (dos <= 1.7), axis=1)
        exp_het = 2 * p_hat * (1 - p_hat)
        excess = obs_het > exp_het
        result["HWE_excess_het_frac"] = round(float(np.nanmean(excess)), 4)

    # Fst (Hudson)
    if "Fst" in what and dos is not None and groups and len(groups) >= 2:
        grp_names = list(groups.keys())
        fst_out = {}
        for ai in range(len(grp_names) - 1):
            for bi in range(ai + 1, len(grp_names)):
                ga, gb = grp_names[ai], grp_names[bi]
                idx_a = [i for i, s in enumerate(sample_ids) if s in set(groups[ga])]
                idx_b = [i for i, s in enumerate(sample_ids) if s in set(groups[gb])]
                if len(idx_a) < 2 or len(idx_b) < 2:
                    continue
                p1 = np.nanmean(dos[:, idx_a], axis=1) / 2.0
                p2 = np.nanmean(dos[:, idx_b], axis=1) / 2.0
                n1, n2 = len(idx_a), len(idx_b)
                num = (p1 - p2)**2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
                den = p1*(1-p2) + p2*(1-p1)
                valid = den > EPS_Q
                fst = float(np.sum(num[valid]) / np.sum(den[valid])) if np.any(valid) else float('nan')
                fst_out[f"Fst_{ga}_{gb}"] = round(max(0.0, fst), 4)
        result["Fst"] = fst_out

    # dXY
    if "dXY" in what and dos is not None and groups and len(groups) >= 2:
        grp_names = list(groups.keys())
        dxy_out = {}
        for ai in range(len(grp_names) - 1):
            for bi in range(ai + 1, len(grp_names)):
                ga, gb = grp_names[ai], grp_names[bi]
                idx_a = [i for i, s in enumerate(sample_ids) if s in set(groups[ga])]
                idx_b = [i for i, s in enumerate(sample_ids) if s in set(groups[gb])]
                if len(idx_a) < 2 or len(idx_b) < 2:
                    continue
                p1 = np.nanmean(dos[:, idx_a], axis=1) / 2.0
                p2 = np.nanmean(dos[:, idx_b], axis=1) / 2.0
                dxy_site = p1 * (1 - p2) + p2 * (1 - p1)
                dxy_out[f"dXY_{ga}_{gb}"] = round(float(np.nanmean(dxy_site)), 6)
        result["dXY"] = dxy_out

    return result


def _load_dosage(chrom, start, end):
    """Load dosage matrix from BEAGLE file for a region."""
    beagle = _resolve_beagle(chrom)

    opener = gzip.open if beagle.endswith(".gz") else open
    with opener(beagle, "rt") as f:
        header = f.readline().strip().split("\t")
        n_ind = (len(header) - 3) // 3

        # Extract sample names
        sample_ids = []
        for i in range(n_ind):
            nm = header[3 + i * 3]
            us = nm.rfind("_")
            if us > 0:
                nm = nm[:us]
            if not sample_ids or sample_ids[-1] != nm:
                sample_ids.append(nm)

        dos_rows = []
        pos_list = []

        for line in f:
            parts = line.strip().split("\t")
            marker = parts[0]
            us = marker.rfind("_")
            if us < 0:
                continue
            pos = int(marker[us + 1:])
            if pos < start or pos > end:
                continue

            pos_list.append(pos)
            gl = [max(float(x), EPS_Q) for x in parts[3:]]
            # Expected dosage
            row = []
            for i in range(n_ind):
                row.append(gl[i*3 + 1] + 2 * gl[i*3 + 2])
            dos_rows.append(row)

    if not dos_rows:
        return None, None, None

    dos = np.array(dos_rows, dtype=np.float64)
    return dos, sample_ids[:n_ind], np.array(pos_list)


def _read_tsv(path):
    """Read TSV without pandas."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline().strip().split("\t")
        rows = []
        for line in f:
            vals = line.strip().split("\t")
            rows.append(dict(zip(header, vals)))
    return rows


# =============================================================================
# CLI
# =============================================================================

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Instant Q Engine (Python CLI)")
    parser.add_argument("--config", help="Config file path")
    parser.add_argument("--beagle_dir", help="BEAGLE directory")
    parser.add_argument("--fopt", help="F matrix file")
    parser.add_argument("--qinit", help="Q init file")
    parser.add_argument("--outdir", default="local_Q", help="Cache directory")
    parser.add_argument("--chr", help="Chromosome")
    parser.add_argument("--start", type=int, help="Start bp")
    parser.add_argument("--end", type=int, help="End bp")
    parser.add_argument("--precompute", action="store_true")
    parser.add_argument("--all", action="store_true", help="All chromosomes")
    parser.add_argument("--K", type=int)
    parser.add_argument("--window_size", type=int, default=100)
    parser.add_argument("--window_step", type=int, default=20)
    parser.add_argument("--em_iter", type=int, default=20)
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--force", action="store_true")

    args = parser.parse_args()

    configure(
        config_file=args.config,
        beagle_dir=args.beagle_dir,
        fopt_file=args.fopt,
        qinit_file=args.qinit,
        cache_dir=args.outdir,
        window_size=args.window_size,
        window_step=args.window_step,
        em_iter=args.em_iter,
        K=args.K,
        ncores=args.ncores,
    )

    if args.precompute:
        if args.all:
            # Discover chromosomes from beagle dir
            bd = _config["beagle_dir"]
            for bf in sorted(glob.glob(os.path.join(bd, "*.beagle.gz"))):
                bn = os.path.basename(bf)
                import re
                m = re.search(r"(C_gar_LG\d+)", bn)
                if m:
                    precompute(m.group(1), force=args.force)
        elif args.chr:
            precompute(args.chr, force=args.force)
        else:
            print("Specify --chr or --all for precompute mode", file=sys.stderr)
            sys.exit(1)
    elif args.chr and args.start is not None and args.end is not None:
        result = get_Q(args.chr, args.start, args.end)
        if HAS_PANDAS:
            print(result.to_csv(sep="\t", index=False))
        else:
            for row in result:
                print("\t".join(str(v) for v in row.values()))
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
