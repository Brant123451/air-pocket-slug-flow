"""Side-by-side comparison: LASSI 5.1 BEFORE vs AFTER dalpha/dh geom fix."""
from pathlib import Path
import matplotlib.pyplot as plt


def load_snap(path):
    text = Path(path).read_text(errors="ignore").splitlines()
    t = None
    xs, betas = [], []
    for line in text:
        if line.startswith("# t ="):
            try:
                t = float(line.split("=")[1])
            except Exception:
                pass
            continue
        s = line.split()
        if len(s) < 9 or s[1] != "SECTION":
            continue
        zR = float(s[2])
        L = float(s[3])
        beta = float(s[4])
        xs.append(zR - 0.5 * L)
        betas.append(beta)
    return t, xs, betas


def panel(ax, dirpath, title):
    snaps = sorted(Path(dirpath).glob("paper51_KH_wave_snap_*.dat"))
    for p in snaps:
        t, xs, betas = load_snap(p)
        if not xs:
            continue
        if t is None:
            lab = p.stem
            lw, alpha = 0.55, 0.25
        elif abs(t - 6.0) < 0.05:
            lab, lw, alpha = "t=6.0s", 1.4, 1.0
        elif abs(t - 4.0) < 0.05:
            lab, lw, alpha = "t=4.0s", 1.1, 0.95
        elif abs(t - 2.0) < 0.05:
            lab, lw, alpha = "t=2.0s", 1.1, 0.95
        elif abs(t - 0.0) < 0.05:
            lab, lw, alpha = "t=0.0s", 1.1, 0.95
        else:
            lab, lw, alpha = None, 0.55, 0.25
        ax.plot(xs, betas, lw=lw, alpha=alpha, label=lab)
    ax.axhline(0.5311, color="k", ls="--", lw=0.9, label=r"base $\alpha_l=0.5311$")
    ax.set_xlabel("x (m)")
    ax.set_ylabel(r"$\alpha_l$")
    ax.set_title(title, fontsize=10)
    ax.set_xlim(0, 30)
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=7, ncol=2)


fig, axs = plt.subplots(2, 1, figsize=(11.5, 8.5), constrained_layout=True)
panel(
    axs[0],
    "outputs_paper_tests_51",
    r"(a) BEFORE geom fix: $d\alpha/dh_l = 2\sin\gamma/(\pi D)$ (factor-2 missing)  $\Rightarrow$ $\kappa$ underestimated by 2$\times$",
)
panel(
    axs[1],
    "outputs_paper_tests_51_fix",
    r"(b) AFTER geom fix: $d\alpha/dh_l = 4\sin\gamma/(\pi D)$ (LASSI thesis App. D, Table D.1)  $\Rightarrow$ $\kappa$ correct",
)
fig.suptitle(
    "LASSI Paper Section 5.1 — effect of dalpha/dh_l geom fix on KH wave evolution\n"
    "(N=300, t<=6s, P-M step still a stub: gas pressure pinned at p_out)",
    fontsize=11,
)
out = Path("outputs_paper_tests_51_fix/compare_before_after_geom_fix.png")
fig.savefig(out, dpi=150)
print(out)
