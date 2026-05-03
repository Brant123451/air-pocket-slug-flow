from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs_paper1_fig18_uniform_no_perturb"
FIG = ROOT / "figs_paper1_fig18_uniform_no_perturb"
FIG.mkdir(exist_ok=True)

SNAP_RE = re.compile(r"_snap_(\d+)\.dat$")


def load_snapshot(path: Path):
    t = 0.0
    beta_xs = []
    beta_vals = []
    section_blocks = []
    slug_blocks = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        first = f.readline()
        second = f.readline()
        if "=" in second:
            try:
                t = float(second.split("=")[-1])
            except ValueError:
                t = 0.0
        for line in f:
            parts = line.split()
            if len(parts) < 7:
                continue
            kind = parts[1]
            if kind not in {"SECTION", "SLUG"}:
                continue
            zR = float(parts[2])
            L = float(parts[3])
            beta = float(parts[4])
            if kind == "SLUG":
                beta = 1.0
            zL = zR - L
            if L <= 0.0:
                continue
            beta_xs.extend([zL, zR, np.nan])
            beta_vals.extend([beta, beta, np.nan])
            if kind == "SLUG":
                slug_blocks.append((zL, zR))
            else:
                section_blocks.append((zL, zR, beta))
    return t, np.array(beta_xs), np.array(beta_vals), section_blocks, slug_blocks


def draw_profile(ax, xs, betas, section_blocks, slug_blocks):
    for zL, zR, beta in section_blocks:
        ax.fill_between([zL, zR], [0.0, 0.0], [beta, beta], color="0.72", linewidth=0.0)
    for zL, zR in slug_blocks:
        ax.axvspan(zL, zR, ymin=0.0, ymax=1.0, color="tab:red", alpha=0.35, linewidth=0.0)
    ax.plot(xs, betas, color="black", lw=1.4)


def snapshots():
    files = sorted(OUT.glob("*_snap_*.dat"), key=lambda p: int(SNAP_RE.search(p.name).group(1)))
    return files


def make_slices(files):
    targets = [0.003, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 4.5]
    loaded = [(p, *load_snapshot(p)) for p in files]
    chosen = []
    for target in targets:
        p, t, xs, betas, section_blocks, slug_blocks = min(loaded, key=lambda item: abs(item[1] - target))
        if p not in [c[0] for c in chosen]:
            chosen.append((p, t, xs, betas, section_blocks, slug_blocks))
    fig, axes = plt.subplots(len(chosen), 1, figsize=(11, 1.7 * len(chosen)), sharex=True)
    if len(chosen) == 1:
        axes = [axes]
    for ax, (p, t, xs, betas, section_blocks, slug_blocks) in zip(axes, chosen):
        draw_profile(ax, xs, betas, section_blocks, slug_blocks)
        ax.axhline(0.98, color="tab:red", lw=0.8, ls="--")
        ax.set_ylim(0, 1.05)
        ax.set_xlim(0, 20)
        ax.set_ylabel(r"$\beta$")
        ax.set_title(f"t = {t:.3f} s, slug blocks = {len(slug_blocks)}", loc="left", fontsize=10)
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("x (m)")
    fig.suptitle("LASSI slug initiation slices", y=0.995)
    fig.tight_layout()
    out = FIG / "fig18_slug_init_slices.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def make_gif(files):
    stride = max(1, len(files) // 70)
    use_files = files[::stride]
    loaded = [load_snapshot(p) for p in use_files]
    fig, ax = plt.subplots(figsize=(10, 3.2))
    line, = ax.plot([], [], color="black", lw=1.5)
    artists = [[]]
    title = ax.text(0.02, 0.92, "", transform=ax.transAxes)
    ax.axhline(0.98, color="tab:red", lw=0.9, ls="--")
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("x (m)")
    ax.set_ylabel(r"$\beta$")
    ax.grid(True, alpha=0.25)

    def update(i):
        for artist in artists[0]:
            artist.remove()
        artists[0] = []
        t, xs, betas, section_blocks, slug_blocks = loaded[i]
        for zL, zR, beta in section_blocks:
            artists[0].append(ax.fill_between([zL, zR], [0.0, 0.0], [beta, beta], color="0.72", linewidth=0.0))
        for zL, zR in slug_blocks:
            artists[0].append(ax.axvspan(zL, zR, ymin=0.0, ymax=1.0, color="tab:red", alpha=0.35, linewidth=0.0))
        line.set_data(xs, betas)
        title.set_text(f"t = {t:.3f} s, slug blocks = {len(slug_blocks)}")
        return [line, title] + artists[0]

    anim = FuncAnimation(fig, update, frames=len(loaded), interval=100, blit=False)
    out = FIG / "fig18_slug_init.gif"
    anim.save(out, writer=PillowWriter(fps=10))
    plt.close(fig)
    return out


def main():
    files = snapshots()
    if not files:
        raise SystemExit(f"No snapshots found in {OUT}")
    slice_path = make_slices(files)
    gif_path = make_gif(files)
    print(slice_path)
    print(gif_path)


if __name__ == "__main__":
    main()
