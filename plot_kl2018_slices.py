"""
Side-view pipe slices of K-L 2018 case 4 at selected times.
Each panel shows the pipe drawn as a horizontal channel of diameter D;
the liquid is filled from the bottom up to the geometric depth h_l(beta)
in stratified sections, and the entire cross-section is filled wherever
a SLUG sits.  A second panel underneath each pipe view plots beta(z) for
the same instant.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


def gamma_from_beta(beta: np.ndarray) -> np.ndarray:
    """Invert beta = (gamma - sin gamma cos gamma)/pi by bisection."""
    b = np.clip(beta, 0.0, 1.0)
    lo = np.zeros_like(b)
    hi = np.full_like(b, np.pi)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f = (mid - np.sin(mid) * np.cos(mid)) / np.pi - b
        hi = np.where(f > 0.0, mid, hi)
        lo = np.where(f > 0.0, lo, mid)
    return 0.5 * (lo + hi)


def hl_from_beta(beta: np.ndarray, D: float) -> np.ndarray:
    g = gamma_from_beta(beta)
    return 0.5 * D * (1.0 - np.cos(g))


def parse_snapshot(path: Path):
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
            zR = float(cols[2])
            L = float(cols[3])
            beta = float(cols[4])
            Ul = float(cols[5])
            rows.append((kind, zR, L, beta, Ul))
    return t, rows


def draw_pipe_panel(ax_pipe, ax_beta, rows, Lpipe, D, t,
                    slug_minwidth=0.08):
    """Draw the pipe side view + beta(z) using object-native (zL, zR)
    rectangles, so thin slugs remain visible at 36 m pipe scale."""
    # pipe outline
    ax_pipe.add_patch(Rectangle((0.0, 0.0), Lpipe, D,
                                fill=False, edgecolor="black", lw=1.0))

    # background: ambient gas (white) is the default

    # collect the section poly-line for the liquid surface and a list of
    # slug rectangles (visually widened to slug_minwidth so single-cell
    # slugs do not vanish on a 36 m axis)
    surface_x: list[float] = []
    surface_y: list[float] = []
    slug_rects: list[tuple[float, float]] = []  # (zL_visual, L_visual)
    sections_for_beta: list[tuple[float, float, float]] = []
    slugs_for_beta: list[tuple[float, float]] = []

    for kind, zR, L, b, _ in rows:
        if kind not in ("SECTION", "SLUG"):
            continue
        if L <= 0.0:
            continue
        zL = zR - L
        if kind == "SLUG":
            visual_L = max(L, slug_minwidth)
            visual_zL = zL - 0.5 * (visual_L - L)
            slug_rects.append((visual_zL, visual_L))
            slugs_for_beta.append((zL, zR))
        else:
            h = float(hl_from_beta(np.array([b]), D)[0])
            ax_pipe.add_patch(Rectangle((zL, 0.0), L, h,
                                        color="#6cb4f2", lw=0))
            surface_x.extend([zL, zR])
            surface_y.extend([h, h])
            sections_for_beta.append((zL, zR, b))

    # slug rectangles drawn last so they sit on top
    for zL_v, Lv in slug_rects:
        ax_pipe.add_patch(Rectangle((zL_v, 0.0), Lv, D, color="#1f4e8a", lw=0))

    # liquid-surface guide line on the pipe view
    if surface_x:
        ax_pipe.plot(surface_x, surface_y, color="#0a3060", lw=0.6, alpha=0.6)

    ax_pipe.set_xlim(0.0, Lpipe)
    ax_pipe.set_ylim(-0.005, D + 0.005)
    ax_pipe.set_yticks([0.0, D])
    ax_pipe.set_yticklabels(["0", f"{D*1000:.0f} mm"])
    ax_pipe.set_xticklabels([])
    ax_pipe.set_aspect("auto")
    ax_pipe.set_title(f"t = {t:.2f} s", loc="left", fontsize=10)

    # beta(z) panel: stratified pieces as a step plot, slugs as full bars
    for zL, zR, b in sections_for_beta:
        ax_beta.hlines(b, zL, zR, color="#1f4e8a", lw=1.0)
    for zL, zR in slugs_for_beta:
        Lv = max(zR - zL, slug_minwidth)
        ax_beta.add_patch(Rectangle((zL - 0.5 * (Lv - (zR - zL)), 0.0),
                                    Lv, 1.0, color="#1f4e8a", alpha=0.85, lw=0))
    ax_beta.set_ylim(0.0, 1.05)
    ax_beta.set_xlim(0.0, Lpipe)
    ax_beta.set_yticks([0.0, 0.5, 1.0])
    ax_beta.set_ylabel(r"$\beta$", fontsize=9)
    ax_beta.grid(alpha=0.25)
    # annotate slug count
    ax_beta.text(0.99, 0.85, f"{len(slugs_for_beta)} slugs",
                 transform=ax_beta.transAxes, ha="right", va="top",
                 fontsize=8, color="#1f4e8a")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snap-dir", default="outputs")
    ap.add_argument("--tag", default="KL2018_case4")
    ap.add_argument("--Lpipe", type=float, default=36.0)
    ap.add_argument("--D", type=float, default=0.078)
    ap.add_argument("--nz", type=int, default=720)
    ap.add_argument("--out", default="outputs/KL2018_case4_slices.png")
    ap.add_argument("--snap-indices", type=int, nargs="+",
                    default=[0, 7, 14, 21, 30, 41])
    args = ap.parse_args()

    snap_dir = Path(args.snap_dir)

    n_panels = len(args.snap_indices)
    fig, axes = plt.subplots(2 * n_panels, 1, figsize=(11.0, 1.8 * n_panels),
                             dpi=160,
                             gridspec_kw={"height_ratios": [1.6, 1.0] * n_panels,
                                          "hspace": 0.35})

    for k, idx in enumerate(args.snap_indices):
        path = snap_dir / f"{args.tag}_snap_{idx:05d}.dat"
        if not path.exists():
            print(f"[skip] {path}")
            continue
        t, rows = parse_snapshot(path)
        ax_pipe = axes[2 * k]
        ax_beta = axes[2 * k + 1]
        draw_pipe_panel(ax_pipe, ax_beta, rows, args.Lpipe, args.D, t or 0.0)
        if k == n_panels - 1:
            ax_beta.set_xlabel("z [m]")
        else:
            ax_beta.set_xticklabels([])

    fig.suptitle("LASSI K-L 2018 case 4: pipe-side view at selected times "
                 r"(light blue = stratified $h_l(\beta)$, dark = slug)",
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(args.out)
    fig.savefig(args.out.replace(".png", ".pdf"))
    plt.close(fig)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
