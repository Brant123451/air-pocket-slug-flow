"""Plot LASSI Appendix C case 1 holdup profiles for Fig C.2 comparison."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "outputs_appendixC_case1"
FIG_DIR = ROOT / "figs_appendixC"
FIG_DIR.mkdir(exist_ok=True)
TAG = "appendixC_case1"
L_PIPE = 2.5
TARGETS = [0.65, 0.70, 0.75, 0.80, 0.95]


def load_curve(path: Path):
    with path.open() as f:
        lines = f.readlines()
    t_val = float(lines[1].split("=")[-1])
    rows = []
    for ln in lines[2:]:
        p = ln.split()
        if len(p) < 6:
            continue
        kind = p[1]
        if kind not in ("SECTION", "SLUG"):
            continue
        try:
            zR = float(p[2]); Lsec = float(p[3]); beta = float(p[4])
        except ValueError:
            continue
        if Lsec <= 0:
            continue
        rows.append((zR - Lsec, zR, beta if kind == "SECTION" else 1.0))
    rows.sort(key=lambda r: r[0])
    xs, ys = [], []
    for zL, zR, alpha in rows:
        xs.extend([zL, zR])
        ys.extend([alpha, alpha])
    return t_val, np.array(xs), np.array(ys)


def time_index(snap_dir: Path):
    out = []
    for p in sorted(snap_dir.glob(f"{TAG}_snap_*.dat")):
        with p.open() as f:
            f.readline()
            t = float(f.readline().split("=")[-1])
        out.append((t, p))
    return out


def pick(target, idx):
    return min(idx, key=lambda kv: abs(kv[0] - target))


def main():
    idx = time_index(OUT_DIR)
    sel = [pick(t, idx) for t in TARGETS]
    n = len(sel)
    fig, axes = plt.subplots(n, 1, figsize=(7.5, 1.6 * n), sharex=True)
    for ax, tt, (ts, path) in zip(axes, TARGETS, sel):
        t_val, x, beta = load_curve(path)
        h = np.clip(beta, 0, 1)
        ax.fill_between(x, 0, h, step="pre", color="#1f77b4", lw=0)
        ax.fill_between(x, h, 1, step="pre", color="#d0d0d0", lw=0)
        ax.set_title(f"Simulated liquid height profile (LASSI) after time: {ts:.2f} s",
                     fontsize=9)
        ax.set_xlim(0, L_PIPE)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("adim. liquid height", fontsize=8)
    axes[-1].set_xlabel("position in the pipe (m)")
    fig.tight_layout()
    out_png = FIG_DIR / "appendixC_case1_figC2_style.png"
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
    print(f"saved -> {out_png}")
    print("time stats:")
    for tt, (ts, path) in zip(TARGETS, sel):
        t_val, x, beta = load_curve(path)
        print(f"  t_target={tt:.2f}  t_snap={ts:.3f}  beta_min={beta.min():.4f}  beta_max={beta.max():.4f}")


if __name__ == "__main__":
    main()
