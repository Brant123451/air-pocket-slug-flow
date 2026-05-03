#!/usr/bin/env python3
"""
analyze_slug_pdf.py
===================

Offline post-processor for the per-slug tracking CSV produced by LASSI
when ``enable_slug_track = T`` is set in the case input file.

Computes:
    * slug-length distribution (mean, std, P10/P50/P90, histogram)
    * inter-slug-arrival distribution at each ``x_obs`` station
    * slug-frequency f_s = 1 / mean(inter-arrival time) per station

This is the offline replacement for the Fortran-level diagnostic
originally listed as TIER-D.1 in the LASSI/thesis-alignment audit.
The Fortran-side ``write_slug_track`` already produces the raw data
(t, slug_id, zL, zR, L, U_m), which contains strictly more information
than any pre-aggregated PDF, so doing the histogramming offline keeps
the simulation kernel small and lets the user re-bin without re-running.

Usage
-----
    python analyze_slug_pdf.py outputs/<tag>_slug_track.csv \
        [--x-obs X1 X2 X3 ...] [--out-dir outputs/<tag>_slug_pdf]

Example
-------
    python analyze_slug_pdf.py outputs/paper52_case4_slug_track.csv \
        --x-obs 9 18 27 34
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_track(csv_path: Path) -> pd.DataFrame:
    """Load the slug_track CSV (whitespace-delimited, hash-comment header)."""
    cols = ["t", "slug_id", "zL", "zR", "L", "Um"]
    df = pd.read_csv(
        csv_path,
        comment="#",
        sep=r"\s+",
        names=cols,
        engine="python",
    )
    df = df.dropna()
    df["slug_id"] = df["slug_id"].astype(int)
    return df


def length_pdf(df: pd.DataFrame, out_dir: Path, tag: str) -> dict:
    """Per-slug length statistics: take peak length over each slug's life."""
    grp = df.groupby("slug_id")["L"]
    peak_L = grp.max()
    if len(peak_L) == 0:
        return {"n_slugs": 0}
    stats = {
        "n_slugs": int(len(peak_L)),
        "L_mean": float(peak_L.mean()),
        "L_std":  float(peak_L.std()),
        "L_p10":  float(np.percentile(peak_L, 10)),
        "L_p50":  float(np.percentile(peak_L, 50)),
        "L_p90":  float(np.percentile(peak_L, 90)),
        "L_max":  float(peak_L.max()),
    }
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    bins = np.linspace(0, max(peak_L.max() * 1.05, 0.1), 30)
    ax.hist(peak_L.values, bins=bins, edgecolor="black", color="#5b8def")
    ax.axvline(stats["L_mean"], ls="--", color="r",
               label=f"mean = {stats['L_mean']:.2f} m")
    ax.axvline(stats["L_p50"], ls=":", color="k",
               label=f"P50 = {stats['L_p50']:.2f} m")
    ax.set_xlabel("Slug peak length L [m]")
    ax.set_ylabel("Count")
    ax.set_title(f"{tag} — slug length PDF (n={stats['n_slugs']})")
    ax.legend()
    fig.savefig(out_dir / f"{tag}_length_pdf.png", dpi=120)
    plt.close(fig)
    peak_L.to_csv(out_dir / f"{tag}_peak_length.csv",
                  header=["L_peak"], index_label="slug_id")
    return stats


def slug_passages_at_x(df: pd.DataFrame, x_obs: float) -> np.ndarray:
    """Return ascending list of times at which a slug crosses x_obs.

    A "crossing" is detected as the first sample of a given slug_id whose
    front (zR) exceeds x_obs while its back (zL) is still upstream (zL <= x_obs).
    Slugs that are born already containing x_obs are also counted on their
    first appearance.  This gives a deterministic time-of-passage even when
    the snapshot interval is coarser than the slug residence time.
    """
    times = []
    for sid, sub in df.groupby("slug_id"):
        sub = sub.sort_values("t")
        contains = (sub["zL"] <= x_obs) & (sub["zR"] >= x_obs)
        if not contains.any():
            continue
        t_first = sub.loc[contains, "t"].iloc[0]
        times.append(float(t_first))
    return np.array(sorted(times))


def frequency_pdf(df: pd.DataFrame, x_obs: float, out_dir: Path,
                  tag: str) -> dict:
    """Inter-arrival statistics + slug frequency at x_obs."""
    t_pass = slug_passages_at_x(df, x_obs)
    n = len(t_pass)
    out: dict = {"x_obs": float(x_obs), "n_passages": int(n)}
    if n < 2:
        out["f_s"] = float("nan")
        return out
    dt = np.diff(t_pass)
    out["dt_mean"] = float(dt.mean())
    out["dt_std"]  = float(dt.std())
    out["dt_p50"]  = float(np.percentile(dt, 50))
    out["f_s"]     = 1.0 / float(dt.mean())

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.hist(dt, bins=20, edgecolor="black", color="#a8d18b")
    ax.axvline(out["dt_mean"], ls="--", color="r",
               label=f"mean Δt = {out['dt_mean']:.2f} s\n"
                     f"f_s = {out['f_s']:.3f} Hz")
    ax.set_xlabel("Inter-slug interval Δt [s]")
    ax.set_ylabel("Count")
    ax.set_title(f"{tag} — slug arrivals at x = {x_obs:.1f} m (n={n})")
    ax.legend()
    fig.savefig(out_dir / f"{tag}_freq_x{x_obs:.1f}.png", dpi=120)
    plt.close(fig)
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("track_csv", type=Path,
                   help="Path to <tag>_slug_track.csv produced by LASSI")
    p.add_argument("--x-obs", type=float, nargs="+", default=[9.0, 18.0, 27.0, 34.0],
                   help="Observation stations [m] for slug frequency PDF")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Output directory (default: alongside track_csv)")
    args = p.parse_args(argv)

    if not args.track_csv.exists():
        print(f"ERROR: {args.track_csv} not found", file=sys.stderr)
        return 1

    tag = args.track_csv.stem.replace("_slug_track", "")
    out_dir = args.out_dir or args.track_csv.parent / f"{tag}_slug_pdf"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = read_track(args.track_csv)
    if df.empty:
        print("WARN: track CSV is empty (no slugs were generated)")
        return 0

    print(f"Read {len(df)} slug-track rows, {df['slug_id'].nunique()} unique slugs")

    L_stats = length_pdf(df, out_dir, tag)
    print("\nSlug length statistics:")
    for k, v in L_stats.items():
        print(f"  {k:10s} = {v}")

    print("\nSlug frequency per station:")
    rows = []
    for x in args.x_obs:
        rec = frequency_pdf(df, x, out_dir, tag)
        rows.append(rec)
        if rec["n_passages"] >= 2:
            print(f"  x = {x:6.2f} m  n = {rec['n_passages']:4d}  "
                  f"<Δt> = {rec['dt_mean']:.3f} s  f_s = {rec['f_s']:.4f} Hz")
        else:
            print(f"  x = {x:6.2f} m  n = {rec['n_passages']:4d}  "
                  f"(too few slugs — f_s undefined)")

    pd.DataFrame(rows).to_csv(out_dir / f"{tag}_freq_table.csv", index=False)
    print(f"\nOutputs written to {out_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
