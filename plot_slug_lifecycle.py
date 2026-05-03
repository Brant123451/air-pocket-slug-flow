"""
Plot the slug life cycle as a SIDE VIEW of the pipe.

For every snapshot we compute the local liquid level h_l(z)/D inside
the pipe.  In stratified sections we invert the geometric relation
    beta = (gamma - sin gamma cos gamma)/pi
    h_l  = D/2 * (1 - cos gamma)
and inside slugs we set h_l/D = 1 (slug = liquid contact with pipe top).

The result is one continuous water-level curve along the pipe.  Slugs
appear naturally as plateaus at h_l/D = 1.  We highlight the tracked
slug #48 with a thick black outline and annotation.

Output: figs/slug_lifecycle_kl2018_case4.png  (slice plot)
        figs/slug_lifecycle_track_kl2018_case4.png  (z(t) tracking)
"""
import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
SNAP_DIR = os.path.join(ROOT, 'outputs')
OUT_DIR  = os.path.join(ROOT, 'figs')
os.makedirs(OUT_DIR, exist_ok=True)

LPIPE = 36.0
TARGET_ID = 48
TARGET_TIMES = [3.9, 4.5, 6.0, 8.0, 10.0, 12.0, 14.0, 15.5, 16.5, 17.4]


# --------------------------------------------------------------------
# Geometry: invert beta -> gamma -> h_l/D  (matches lassi_geom.f90)
# --------------------------------------------------------------------
def gamma_from_beta(beta):
    """Bisect  beta = (gamma - sin(gamma) cos(gamma))/pi  for gamma in [0, pi]."""
    beta = np.clip(beta, 0.0, 1.0)
    lo = 0.0
    hi = np.pi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f = (mid - np.sin(mid) * np.cos(mid)) / np.pi - beta
        # f is monotonic increasing in mid; if f<0 we need larger gamma
        lo = np.where(f < 0, mid, lo)
        hi = np.where(f < 0, hi,  mid)
    return 0.5 * (lo + hi)


def h_over_D_from_beta(beta):
    g = gamma_from_beta(np.asarray(beta, dtype=float))
    return 0.5 * (1.0 - np.cos(g))


# --------------------------------------------------------------------
# Snapshot parsing
# --------------------------------------------------------------------
def parse_snapshot(path):
    objs = []
    t = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                m = re.search(r't\s*=\s*([0-9.+\-Ee]+)', s)
                if m:
                    t = float(m.group(1))
                continue
            cols = s.split()
            if len(cols) < 10:
                continue
            try:
                objs.append({
                    'id'  : int(cols[0]),
                    'kind': cols[1],
                    'zR'  : float(cols[2]),
                    'L'   : float(cols[3]),
                    'beta': float(cols[4]),
                    'Ul'  : float(cols[5]),
                    'Um'  : float(cols[6]),
                    'p'   : float(cols[7]),
                    'rhog': float(cols[8]),
                })
            except ValueError:
                continue
    return t, objs


def liquid_level_curve(objs):
    """Build a piecewise water-level curve (z, h/D) from the object list."""
    zs = []
    hs = []
    ids = []
    kinds = []
    for o in objs:
        zL = o['zR'] - o['L']
        zR = o['zR']
        if o['kind'] == 'SECTION':
            h_over_D = h_over_D_from_beta(o['beta'])
        elif o['kind'] == 'SLUG':
            h_over_D = 1.0
        elif o['kind'] in ('INLET', 'OUTLET'):
            continue
        elif o['kind'] == 'BUBBLE':
            # bubbles are point-like (zL == zR); skip in level curve
            continue
        else:
            continue
        # piecewise constant within object: emit two points (zL, zR)
        zs.append(zL); hs.append(h_over_D); ids.append(o['id']); kinds.append(o['kind'])
        zs.append(zR); hs.append(h_over_D); ids.append(o['id']); kinds.append(o['kind'])
    return np.array(zs), np.array(hs), np.array(ids), np.array(kinds)


def slug_rects(objs):
    out = []
    for o in objs:
        if o['kind'] == 'SLUG':
            zL = o['zR'] - o['L']
            out.append((zL, o['L'], o['Um'], o['id']))
    return out


def find_snapshot_at(snaps, t_target):
    return min(snaps, key=lambda x: abs(x[0] - t_target))


# --------------------------------------------------------------------
# Main figure
# --------------------------------------------------------------------
def main():
    files = sorted(glob.glob(os.path.join(SNAP_DIR, 'KL2018_case4_snap_*.dat')))
    snaps = []
    for f in files:
        t, objs = parse_snapshot(f)
        if t is None:
            continue
        snaps.append((t, objs))
    print(f'Loaded {len(snaps)} snapshots, t={snaps[0][0]:.2f} -> {snaps[-1][0]:.2f}s')

    n = len(TARGET_TIMES)
    fig, axes = plt.subplots(
        n, 2, figsize=(15.5, 1.6 * n),
        gridspec_kw={'width_ratios': [4.5, 1.0], 'wspace': 0.06,
                     'hspace': 0.4},
        sharey='row',
    )
    if n == 1:
        axes = np.array([axes])

    def draw_panel(ax, objs, xL, xR, draw_pipewalls=True, zoom=False):
        if draw_pipewalls:
            ax.axhline(1.0, color='dimgray', lw=1.2, zorder=1)
            ax.axhline(0.0, color='dimgray', lw=1.2, zorder=1)
        zs, hs, ids, kinds = liquid_level_curve(objs)
        if len(zs) > 0:
            ax.fill_between(zs, 0.0, hs, color='#3b8bd9', alpha=0.55,
                            step='pre', zorder=2, linewidth=0)
            ax.plot(zs, hs, color='#003e7a', lw=0.9 if not zoom else 1.4,
                    drawstyle='steps-pre', zorder=3)
        # Outline the tracked slug
        for zL, L, Um, sid in slug_rects(objs):
            if sid == TARGET_ID:
                ax.add_patch(plt.Rectangle(
                    (zL, 0.0), L, 1.0,
                    facecolor='none', edgecolor='black', linewidth=2.5,
                    zorder=5,
                ))
        ax.set_xlim(xL, xR)

    for row, t_t in enumerate(TARGET_TIMES):
        ax_full = axes[row, 0]
        ax_zoom = axes[row, 1]
        t, objs = find_snapshot_at(snaps, t_t)

        # ----- full-pipe panel -----
        draw_panel(ax_full, objs, 0.0, LPIPE)

        # tracked slug annotation on full panel
        target_z = None; target_um = None; target_L = None
        for zL, L, Um, sid in slug_rects(objs):
            if sid == TARGET_ID:
                target_z = zL + 0.5 * L
                target_um = Um
                target_L = L
        if target_z is not None:
            ax_full.annotate(
                f'#{TARGET_ID}: $z$={target_z:.1f}m  $L$={target_L:.3f}m  $U_m$={target_um:.2f} m/s',
                xy=(target_z, 1.0), xytext=(target_z, 1.55),
                ha='center', va='bottom', fontsize=8, color='black',
                arrowprops=dict(arrowstyle='->', color='black', lw=0.8),
                zorder=6,
            )
            # zoom window of ±1 m around tracked slug
            zoom_xL = max(0.0,   target_z - 1.0)
            zoom_xR = min(LPIPE, target_z + 1.0)
        else:
            zoom_xL, zoom_xR = 0.0, 2.0  # arbitrary fallback

        ax_full.set_ylim(-0.05, 1.85)
        ax_full.set_ylabel(f't = {t:.2f} s', rotation=0, ha='right', va='center', fontsize=10)
        ax_full.tick_params(axis='y', labelsize=8)
        ax_full.set_yticks([0.0, 0.5, 1.0])
        ax_full.set_yticklabels(['bottom', '0.5 D', 'top'])
        ax_full.set_xlim(0.0, LPIPE)
        if row != n - 1:
            ax_full.set_xticklabels([])

        # ----- zoom panel -----
        draw_panel(ax_zoom, objs, zoom_xL, zoom_xR, zoom=True)
        ax_zoom.set_ylim(-0.05, 1.85)
        ax_zoom.tick_params(axis='y', labelsize=8)
        ax_zoom.set_yticks([0.0, 0.5, 1.0])
        ax_zoom.set_yticklabels([])
        # zoom range annotation inside panel
        ax_zoom.text(0.02, 0.94, f'zoom: z={zoom_xL:.1f}-{zoom_xR:.1f}m',
                     transform=ax_zoom.transAxes, fontsize=7, color='dimgray',
                     ha='left', va='top',
                     bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
        ax_zoom.tick_params(axis='x', labelsize=7)

        # Mark zoom window on the full panel
        ax_full.add_patch(plt.Rectangle(
            (zoom_xL, -0.05), zoom_xR - zoom_xL, 1.85 + 0.05,
            facecolor='none', edgecolor='orange', lw=0.7, ls='--', zorder=4,
        ))

    axes[-1, 0].set_xlabel('z (m) — full pipe', fontsize=11)
    axes[-1, 1].set_xlabel('z (m) — zoom', fontsize=10)

    # Title via figtext (no tight_layout conflict)
    fig.text(0.5, 0.985,
             'LASSI side-view of liquid level $h_l/D$ — K-L 2018 case 4 '
             '(USL=1, USG=2 m/s, D=0.078 m, L=36 m)',
             ha='center', va='top', fontsize=12, fontweight='bold')
    fig.text(0.5, 0.965,
             f'tracked slug #{TARGET_ID}: born at z=8.8 m (t=3.9 s), '
             'exits at outlet (t=17.4 s) — zoom panel (right) shows '
             'a ±1 m window around it',
             ha='center', va='top', fontsize=10, color='dimgray')

    # Figure-level legend
    blue_patch = plt.Rectangle((0, 0), 1, 1, facecolor='#3b8bd9', alpha=0.55,
                               edgecolor='#003e7a', linewidth=0.9)
    track_patch = plt.Rectangle((0, 0), 1, 1, facecolor='none',
                                edgecolor='black', linewidth=2.5)
    zoom_patch = plt.Rectangle((0, 0), 1, 1, facecolor='none',
                               edgecolor='orange', linewidth=1.0, ls='--')
    fig.legend(
        [blue_patch, track_patch, zoom_patch],
        ['liquid below $h_l/D$ (slug $\\Rightarrow h_l/D=1$)',
         f'tracked slug #{TARGET_ID}',
         'zoom window'],
        loc='upper center', bbox_to_anchor=(0.5, 0.945),
        ncol=3, fontsize=9, frameon=False,
    )

    fig.subplots_adjust(left=0.06, right=0.99, top=0.93, bottom=0.04)
    out_path = os.path.join(OUT_DIR, 'slug_lifecycle_kl2018_case4.png')
    fig.savefig(out_path, dpi=160, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {out_path}')

    # ------------------------------------------------------------------
    # Tracking plot (z(t), L(t)) for slug #48
    # ------------------------------------------------------------------
    trace = []
    for t, objs in snaps:
        for o in objs:
            if o['kind'] == 'SLUG' and o['id'] == TARGET_ID:
                trace.append((t, o['zR'] - 0.5 * o['L'], o['L'], o['Um']))
                break

    if trace:
        fig2, axL = plt.subplots(1, 1, figsize=(8.5, 4.2))
        ts  = np.array([p[0] for p in trace])
        zs  = np.array([p[1] for p in trace])
        Ls  = np.array([p[2] for p in trace])
        ums = np.array([p[3] for p in trace])

        axL.plot(ts, zs, '-o', color='red', label='slug centroid z(t)', markersize=3)
        axL.set_xlabel('t (s)')
        axL.set_ylabel('z (m)', color='red')
        axL.tick_params(axis='y', labelcolor='red')
        axL.set_ylim(0, LPIPE)

        axR = axL.twinx()
        axR.plot(ts, Ls, '-s', color='blue', label='slug length L(t)', markersize=3)
        axR.set_ylabel('L (m)', color='blue')
        axR.tick_params(axis='y', labelcolor='blue')

        coef = np.polyfit(ts, zs, 1)
        axL.plot(ts, np.polyval(coef, ts), '--', color='gray',
                 label=f'fit U_b = {coef[0]:.2f} m/s')

        axL.set_title(f'Slug #{TARGET_ID} life cycle  (mean U_m = {ums.mean():.2f} m/s, '
                      f'mean L = {Ls.mean():.3f} m)')
        axL.legend(loc='upper left', fontsize=9)
        axR.legend(loc='lower right', fontsize=9)
        axL.grid(True, alpha=0.3)

        plt.tight_layout()
        out2 = os.path.join(OUT_DIR, 'slug_lifecycle_track_kl2018_case4.png')
        fig2.savefig(out2, dpi=160)
        plt.close(fig2)
        print(f'wrote {out2}')


if __name__ == '__main__':
    main()
