"""
Plot LASSI K-L 2018 case-4 x-t maps from snapshot files.
Each snapshot lists the linked-list state at one time; we resample β
onto a uniform z grid and stack the rows to form an x-t heatmap, then
also produce a slug-count time series.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_snapshot(path: Path):
    """Return (t, list of (kind, zR, L, beta, Ul))."""
    t = None
    rows = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                m = re.search(r"t\s*=\s*([0-9.+\-Ee]+)", s)
                if m:
                    t = float(m.group(1))
                continue
            cols = s.split()
            if len(cols) < 8:
                continue
            kind = cols[1]
            zR   = float(cols[2])
            L    = float(cols[3])
            beta = float(cols[4])
            Ul   = float(cols[5])
            rows.append((kind, zR, L, beta, Ul))
    return t, rows


def resample_beta(rows, z_grid):
    """Piecewise-constant resample of beta onto the uniform z_grid."""
    beta = np.full_like(z_grid, np.nan)
    # walk left-to-right; sections/slugs cover (zR-L, zR)
    for kind, zR, L, b, _ in rows:
        if kind not in ("SECTION", "SLUG"):
            continue
        if L <= 0.0:
            continue
        zL = zR - L
        mask = (z_grid >= zL) & (z_grid < zR)
        if kind == "SLUG":
            beta[mask] = 1.0
        else:
            beta[mask] = b
    return beta


def count_slugs(rows):
    return sum(1 for r in rows if r[0] == "SLUG")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snap-glob", default="outputs/KL2018_case4_snap_*.dat")
    ap.add_argument("--Lpipe", type=float, default=36.0)
    ap.add_argument("--nz", type=int, default=720)
    ap.add_argument("--out-prefix", default="outputs/KL2018_case4_xt")
    args = ap.parse_args()

    snap_paths = sorted(Path(".").glob(args.snap_glob))
    if not snap_paths:
        raise SystemExit(f"no snapshots match {args.snap_glob}")

    z_grid = np.linspace(0.0, args.Lpipe, args.nz)
    times, beta_rows, slug_counts = [], [], []
    for p in snap_paths:
        t, rows = parse_snapshot(p)
        if t is None:
            continue
        times.append(t)
        beta_rows.append(resample_beta(rows, z_grid))
        slug_counts.append(count_slugs(rows))

    times = np.array(times)
    beta_arr = np.vstack(beta_rows)

    # x-t heatmap of beta
    fig, ax = plt.subplots(figsize=(7.0, 4.5), dpi=150)
    img = ax.pcolormesh(z_grid, times, beta_arr,
                        cmap="viridis", vmin=0.0, vmax=1.0,
                        shading="auto")
    cb = fig.colorbar(img, ax=ax, label=r"$\beta = A_l / A$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("t [s]")
    ax.set_title(r"LASSI K-L 2018 case 4: liquid holdup $\beta(z,t)$")
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_beta.png")
    fig.savefig(f"{args.out_prefix}_beta.pdf")
    plt.close(fig)

    # slug count time series
    fig, ax = plt.subplots(figsize=(7.0, 3.0), dpi=150)
    ax.plot(times, slug_counts, "-o", lw=1.5, ms=3)
    ax.set_xlabel("t [s]")
    ax.set_ylabel("slug count")
    ax.set_title("LASSI K-L 2018 case 4: number of slugs vs time")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_slugcount.png")
    plt.close(fig)

    print(f"wrote {args.out_prefix}_beta.png/pdf and {args.out_prefix}_slugcount.png")
    print(f"final slug count = {slug_counts[-1]}, peak = {max(slug_counts)}")


if __name__ == "__main__":
    main()
