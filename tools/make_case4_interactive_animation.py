from __future__ import annotations

from pathlib import Path
import html
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs_current_paper52_case4_N200_noinlet_thesis_shedfix"
FIG = ROOT / "figs_current_paper52_case4_N200_noinlet_thesis_shedfix"
TAG = "paper52_case4"
SNAP_RE = re.compile(r"_snap_(\d+)\.dat$")


def load_snapshot(path: Path):
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        f.readline()
        t = float(f.readline().split("=")[-1])
        for line in f:
            p = line.split()
            if len(p) < 7:
                continue
            kind = p[1]
            if kind not in {"SECTION", "SLUG"}:
                continue
            zR = float(p[2])
            length = float(p[3])
            beta = float(p[4])
            if length <= 0.0:
                continue
            zL = zR - length
            rows.append((zL, zR, 1.0 if kind == "SLUG" else beta, kind))
    rows.sort(key=lambda item: item[0])
    xs = []
    ys = []
    slugs = []
    for zL, zR, beta, kind in rows:
        xs.extend([zL, zR, np.nan])
        ys.extend([beta, beta, np.nan])
        if kind == "SLUG":
            slugs.append((zL, zR))
    return t, np.asarray(xs), np.asarray(ys), slugs


def snapshots():
    files = []
    for path in OUT.glob(f"{TAG}_snap_*.dat"):
        match = SNAP_RE.search(path.name)
        if match:
            files.append((int(match.group(1)), path))
    files.sort(key=lambda item: item[0])
    return [path for _, path in files]


def draw_frame(ax, t, xs, ys, slugs):
    for zL, zR in slugs:
        ax.axvspan(zL, zR, color="tab:red", alpha=0.28, lw=0)
    ax.plot(xs, ys, color="black", lw=1.25)
    ax.axhline(0.98, color="tab:red", lw=0.8, ls="--", alpha=0.65)
    ax.set_xlim(0.0, 36.0)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("x (m)")
    ax.set_ylabel(r"$\alpha_l$")
    ax.grid(True, alpha=0.25, lw=0.4)
    ax.set_title(f"CASE 4 N=200 no inlet perturbation | t = {t:.3f} s | slug blocks = {len(slugs)}", loc="left", fontsize=10)


def make_frames(loaded):
    frame_dir = FIG / "paper52_case4_N200_noinlet_thesis_shedfix_frames"
    frame_dir.mkdir(parents=True, exist_ok=True)
    frame_paths = []
    for i, (t, xs, ys, slugs) in enumerate(loaded):
        fig, ax = plt.subplots(figsize=(11, 3.4))
        draw_frame(ax, t, xs, ys, slugs)
        fig.tight_layout()
        path = frame_dir / f"frame_{i:04d}.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        frame_paths.append(path)
    return frame_dir, frame_paths


def make_gif(loaded):
    fig, ax = plt.subplots(figsize=(11, 3.4))

    def update(i):
        ax.clear()
        draw_frame(ax, *loaded[i])
        fig.tight_layout()

    anim = FuncAnimation(fig, update, frames=len(loaded), interval=180, repeat=True)
    gif_path = FIG / "paper52_case4_N200_noinlet_thesis_shedfix_animation.gif"
    anim.save(gif_path, writer=PillowWriter(fps=6))
    plt.close(fig)
    return gif_path


def make_html(loaded, frame_dir, frame_paths):
    html_path = FIG / "paper52_case4_N200_noinlet_thesis_shedfix_frame_viewer.html"
    rel_paths = [path.relative_to(html_path.parent).as_posix() for path in frame_paths]
    labels = [f"frame {i + 1}/{len(loaded)} | t={loaded[i][0]:.3f} s | slug blocks={len(loaded[i][3])}" for i in range(len(loaded))]
    image_items = ",\n".join(f"        {path!r}" for path in rel_paths)
    label_items = ",\n".join(f"        {label!r}" for label in labels)
    content = f"""<!doctype html>
<html lang=\"zh-CN\">
<head>
<meta charset=\"utf-8\">
<title>LASSI CASE4 frame viewer</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 22px; background: #f6f7f9; color: #111827; }}
.container {{ max-width: 1180px; margin: 0 auto; background: white; padding: 18px 22px; border-radius: 14px; box-shadow: 0 8px 28px rgba(15,23,42,0.12); }}
h1 {{ font-size: 20px; margin: 0 0 8px; }}
.meta {{ color: #4b5563; margin-bottom: 14px; }}
.controls {{ display: flex; gap: 10px; align-items: center; flex-wrap: wrap; margin: 12px 0 16px; }}
button {{ border: 0; background: #2563eb; color: white; padding: 8px 13px; border-radius: 8px; cursor: pointer; font-weight: 600; }}
button.secondary {{ background: #4b5563; }}
input[type=range] {{ flex: 1; min-width: 360px; }}
#label {{ font-weight: 700; color: #111827; min-width: 360px; }}
#frame {{ width: 100%; border: 1px solid #d1d5db; border-radius: 10px; background: white; }}
.hint {{ color: #6b7280; font-size: 13px; margin-top: 10px; }}
</style>
</head>
<body>
<div class=\"container\">
<h1>LASSI CASE4 N=200 无入口扰动：逐帧查看器</h1>
<div class=\"meta\">结果目录：{html.escape(str(OUT))}<br>帧目录：{html.escape(str(frame_dir))}</div>
<div class=\"controls\">
<button onclick=\"prevFrame()\">上一帧</button>
<button onclick=\"nextFrame()\">下一帧</button>
<button class=\"secondary\" onclick=\"togglePlay()\" id=\"playButton\">播放</button>
<input id=\"slider\" type=\"range\" min=\"0\" max=\"{len(rel_paths) - 1}\" value=\"0\" step=\"1\" oninput=\"setFrame(Number(this.value))\">
<div id=\"label\"></div>
</div>
<img id=\"frame\" alt=\"LASSI frame\">
<div class=\"hint\">提示：可以用滑块逐帧拖动，也可以用键盘左右方向键切换帧。</div>
</div>
<script>
const images = [
{image_items}
];
const labels = [
{label_items}
];
let idx = 0;
let timer = null;
function setFrame(i) {{
    idx = Math.max(0, Math.min(images.length - 1, i));
    document.getElementById('frame').src = images[idx];
    document.getElementById('label').textContent = labels[idx];
    document.getElementById('slider').value = idx;
}}
function prevFrame() {{ setFrame(idx - 1); }}
function nextFrame() {{ setFrame(idx + 1); }}
function togglePlay() {{
    const button = document.getElementById('playButton');
    if (timer) {{
        clearInterval(timer);
        timer = null;
        button.textContent = '播放';
    }} else {{
        timer = setInterval(() => setFrame((idx + 1) % images.length), 180);
        button.textContent = '暂停';
    }}
}}
document.addEventListener('keydown', event => {{
    if (event.key === 'ArrowLeft') prevFrame();
    if (event.key === 'ArrowRight') nextFrame();
}});
setFrame(0);
</script>
</body>
</html>
"""
    html_path.write_text(content, encoding="utf-8")
    return html_path


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    files = snapshots()
    if not files:
        raise SystemExit(f"No snapshots found in {OUT}")
    loaded = [load_snapshot(path) for path in files]
    frame_dir, frame_paths = make_frames(loaded)
    gif_path = make_gif(loaded)
    html_path = make_html(loaded, frame_dir, frame_paths)
    print(f"frames={len(frame_paths)}")
    print(f"frame_dir={frame_dir}")
    print(f"gif={gif_path}")
    print(f"html={html_path}")
    for i, (t, _, _, slugs) in enumerate(loaded):
        if i == 0 or i == len(loaded) - 1 or abs((t * 2) - round(t * 2)) < 0.02:
            lengths = [zR - zL for zL, zR in slugs]
            print(f"frame={i:03d} t={t:.3f} nslug={len(slugs)} max_slug_L={(max(lengths) if lengths else 0.0):.3f}")


if __name__ == "__main__":
    main()
