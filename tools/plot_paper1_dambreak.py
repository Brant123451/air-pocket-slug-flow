"""Compare LASSI numerical dam-break against analytical MSW Riemann solution.
Reproduces Paper 1 §4.2 Figure 12 verification (LASSI vs analytical).
"""
from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "outputs_paper1_dambreak"
FIG_DIR = ROOT / "figs_paper1"
FIG_DIR.mkdir(exist_ok=True)
TAG = "paper1_dambreak"

# ---- Test parameters (must match INPUT_paper1_dambreak.txt) -----------------
D       = 1.5
L_PIPE  = 5.0
X_DAM   = 2.5
BETA_L  = 0.7217
BETA_R  = 0.2973
RHO_L   = 1000.0
RHO_G   = 1.18837   # default rho_g_ref from grid.f90 (and what LASSI uses)
G_ACC   = 9.80665

# Flat-top D-shape geometry
GAMMA_MAX = math.acos(-0.9)
AREA_FAC = GAMMA_MAX - math.sin(GAMMA_MAX) * math.cos(GAMMA_MAX)


def gamma_from_beta(beta: float) -> float:
    lo, hi = 0.0, GAMMA_MAX
    for _ in range(80):
        m = 0.5 * (lo + hi)
        if (m - math.sin(m) * math.cos(m)) / AREA_FAC > beta:
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)


def kappa_lassi(beta: float, U_l: float = 0.0, U_g: float = 0.0) -> float:
    gamma = gamma_from_beta(beta)
    dalpha_dh = 4.0 * math.sin(gamma) / (AREA_FAC * D)
    inv_dalpha = 1.0 / max(dalpha_dh, 1e-12)
    alpha = max(1.0 - beta, 1e-9)
    return (RHO_L - RHO_G) / RHO_L * G_ACC * inv_dalpha \
        - (1.0 / alpha) * (RHO_G / RHO_L) * (U_g - U_l) ** 2


# ---- Analytical Riemann solution for LASSI MSW (κ averaged at the dam) -----
def analytical_riemann(beta_L: float, beta_R: float, U_L: float = 0.0, U_R: float = 0.0):
    kappa = 0.5 * (kappa_lassi(beta_L) + kappa_lassi(beta_R))
    sqk = math.sqrt(kappa)

    def uL_curve(b):
        if b <= beta_L:
            return U_L - 2.0 * sqk * (math.sqrt(b) - math.sqrt(beta_L))
        return U_L - (1.0 / math.sqrt(2.0)) * sqk * (b - beta_L) \
            * math.sqrt(1.0 / b + 1.0 / beta_L)

    def uR_curve(b):
        if b <= beta_R:
            return U_R + 2.0 * sqk * (math.sqrt(b) - math.sqrt(beta_R))
        return U_R + (1.0 / math.sqrt(2.0)) * sqk * (b - beta_R) \
            * math.sqrt(1.0 / b + 1.0 / beta_R)

    lo, hi = 1e-4, 1.0 - 1e-4
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f = uL_curve(mid) - uR_curve(mid)
        if f * (uL_curve(lo) - uR_curve(lo)) < 0:
            hi = mid
        else:
            lo = mid
    beta_M = 0.5 * (lo + hi)
    U_M = 0.5 * (uL_curve(beta_M) + uR_curve(beta_M))

    s_LL = U_L - sqk * math.sqrt(beta_L)        # left rarefaction head speed
    s_LR = U_M - sqk * math.sqrt(beta_M)        # left rarefaction tail speed
    if abs(beta_M - beta_R) < 1e-12:
        s_R = U_M
    else:
        s_R = (beta_M * U_M - beta_R * U_R) / (beta_M - beta_R)
    return dict(kappa=kappa, sqk=sqk, beta_M=beta_M, U_M=U_M,
                s_LL=s_LL, s_LR=s_LR, s_R=s_R)


def analytical_profile(t: float, x_arr: np.ndarray):
    rs = analytical_riemann(BETA_L, BETA_R)
    sqk = rs["sqk"]
    beta = np.empty_like(x_arr)
    for i, x in enumerate(x_arr):
        xi = (x - X_DAM) / max(t, 1e-12)        # similarity variable
        if xi <= rs["s_LL"]:
            beta[i] = BETA_L
        elif xi <= rs["s_LR"]:
            # left rarefaction: U + 2*sqrt(kappa*β) = U_L + 2*sqrt(kappa*β_L),
            # and xi = U - sqrt(kappa*β), therefore
            # sqrt(β) = (U_L + 2*sqrt(kappa*β_L) - xi) / (3*sqrt(kappa)).
            sb = (2.0 * math.sqrt(BETA_L) - xi / sqk) / 3.0
            beta[i] = max(min(sb * sb, BETA_L), 1e-6)
        elif xi <= rs["s_R"]:
            beta[i] = rs["beta_M"]
        else:
            beta[i] = BETA_R
    return beta, rs


# ---- Load LASSI snapshots ---------------------------------------------------
def load_snapshot(path: Path):
    with path.open() as f:
        f.readline()
        t_val = float(f.readline().split("=")[-1])
        rows = []
        for ln in f:
            v = ln.split()
            if len(v) < 6:
                continue
            kind = v[1]
            if kind not in ("SECTION", "SLUG"):
                continue
            zR = float(v[2]); Lsec = float(v[3]); beta = float(v[4])
            if Lsec <= 0:
                continue
            rows.append((zR - Lsec, zR, beta if kind == "SECTION" else 1.0))
    rows.sort()
    xs, ys = [], []
    for zL, zR, b in rows:
        xs.extend([zL, zR])
        ys.extend([b, b])
    return t_val, np.asarray(xs), np.asarray(ys)


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
    rs = analytical_riemann(BETA_L, BETA_R)
    print("Analytical Riemann solution for the LASSI MSW dam-break:")
    print(f"  κ_avg     = {rs['kappa']:.4f}")
    print(f"  √(κ)      = {rs['sqk']:.4f}")
    print(f"  β_M       = {rs['beta_M']:.4f}")
    print(f"  U_M       = {rs['U_M']:+.4f} m/s")
    print(f"  s_LL  (head)  = {rs['s_LL']:+.4f} m/s   (left rarefaction head)")
    print(f"  s_LR  (tail)  = {rs['s_LR']:+.4f} m/s   (left rarefaction tail)")
    print(f"  s_R   (shock) = {rs['s_R']:+.4f} m/s   (right shock)")

    idx = time_index(OUT_DIR)
    targets = [0.06, 0.10, 0.20, 0.30, 0.40]
    selected = [pick(t, idx) for t in targets]
    n = len(selected)
    fig, axes = plt.subplots(n, 1, figsize=(9, 1.7 * n), sharex=True)
    if n == 1:
        axes = [axes]
    for ax, tt, (ts, path) in zip(axes, targets, selected):
        t_val, x_lassi, beta_lassi = load_snapshot(path)
        x_anal = np.linspace(0, L_PIPE, 1001)
        beta_anal, _ = analytical_profile(t_val, x_anal)
        ax.plot(x_lassi, beta_lassi, "k-", lw=1.6, label="LASSI", drawstyle="steps-pre")
        ax.plot(x_anal, beta_anal, "r--", lw=1.4, label="Analytical")
        ax.axhline(BETA_L, color="grey", lw=0.6, ls=":")
        ax.axhline(BETA_R, color="grey", lw=0.6, ls=":")
        ax.axhline(rs["beta_M"], color="orange", lw=0.6, ls=":")
        ax.set_ylabel(r"$\beta$", fontsize=10)
        ax.set_title(f"t = {ts:.2f} s   (target {tt:.2f} s)",
                     fontsize=10, loc="left")
        ax.set_xlim(1.0, 4.0)
        ax.set_ylim(0.0, 1.0)
        ax.grid(True, lw=0.3, alpha=0.4)
    axes[0].legend(loc="upper right", fontsize=9)
    axes[-1].set_xlabel("x (m)")
    fig.suptitle("Paper 1 §4.2 Moving Dam Break — LASSI vs analytical Riemann",
                 fontsize=11, y=0.995)
    fig.tight_layout()
    out_png = FIG_DIR / "paper1_dambreak_compare.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"saved -> {out_png}")


if __name__ == "__main__":
    main()
