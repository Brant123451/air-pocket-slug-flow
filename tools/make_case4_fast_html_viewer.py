from __future__ import annotations

from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs_current_paper52_case4_N200_noinlet_thesis_shedfix"
FIG = ROOT / "figs_current_paper52_case4_N200_noinlet_thesis_shedfix"
TAG = "paper52_case4"
SNAP_RE = re.compile(r"_snap_(\d+)\.dat$")


def read_snapshot(path: Path):
    segments = []
    slugs = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        f.readline()
        t = float(f.readline().split("=")[-1])
        for line in f:
            parts = line.split()
            if len(parts) < 7:
                continue
            kind = parts[1]
            if kind not in {"SECTION", "SLUG"}:
                continue
            zR = float(parts[2])
            length = float(parts[3])
            beta = float(parts[4])
            if length <= 0.0:
                continue
            zL = zR - length
            value = 1.0 if kind == "SLUG" else beta
            segments.append([zL, zR, value, kind])
            if kind == "SLUG":
                slugs.append([zL, zR])
    segments.sort(key=lambda item: item[0])
    return {"t": t, "segments": segments, "slugs": slugs}


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    files = []
    for path in OUT.glob(f"{TAG}_snap_*.dat"):
        match = SNAP_RE.search(path.name)
        if match:
            files.append((int(match.group(1)), path))
    files.sort(key=lambda item: item[0])
    if not files:
        raise SystemExit(f"No snapshots found in {OUT}")

    frames = [read_snapshot(path) for _, path in files]
    payload = json.dumps(frames, ensure_ascii=False, separators=(",", ":"))
    html = f"""<!doctype html>
<html lang=\"zh-CN\">
<head>
<meta charset=\"utf-8\">
<title>LASSI CASE4 interactive viewer</title>
<style>
body {{ margin: 0; padding: 22px; font-family: Arial, sans-serif; background: #f3f4f6; color: #111827; }}
main {{ max-width: 1180px; margin: auto; background: white; border-radius: 14px; padding: 18px 22px; box-shadow: 0 10px 30px rgba(15,23,42,.12); }}
h1 {{ font-size: 20px; margin: 0 0 8px; }}
.meta {{ color: #4b5563; font-size: 14px; line-height: 1.6; }}
.controls {{ margin: 16px 0; display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }}
button {{ background: #2563eb; color: white; border: 0; border-radius: 8px; padding: 8px 13px; font-weight: 700; cursor: pointer; }}
button.gray {{ background: #4b5563; }}
input[type=range] {{ flex: 1; min-width: 380px; }}
#label {{ min-width: 420px; font-weight: 700; }}
canvas {{ width: 100%; height: 360px; border: 1px solid #d1d5db; border-radius: 10px; background: white; }}
.hint {{ color: #6b7280; font-size: 13px; margin-top: 10px; }}
</style>
</head>
<body>
<main>
<h1>LASSI CASE4 N=200 无入口扰动：交互逐帧动图</h1>
<div class=\"meta\">输出目录：{OUT}<br>帧数：{len(frames)}。红色半透明区域是代码中的 SLUG 对象；黑线是液相体积分数/持液率。</div>
<div class=\"controls\">
<button onclick=\"step(-1)\">上一帧</button>
<button onclick=\"step(1)\">下一帧</button>
<button class=\"gray\" id=\"play\" onclick=\"togglePlay()\">播放</button>
<input id=\"slider\" type=\"range\" min=\"0\" max=\"{len(frames)-1}\" value=\"0\" step=\"1\" oninput=\"setFrame(+this.value)\">
<span id=\"label\"></span>
</div>
<canvas id=\"plot\" width=\"1100\" height=\"360\"></canvas>
<div class=\"hint\">可拖动滑块逐帧查看；键盘左右方向键也可以切换。播放速度约 6 fps。</div>
</main>
<script>
const frames = {payload};
const canvas = document.getElementById('plot');
const ctx = canvas.getContext('2d');
const slider = document.getElementById('slider');
const label = document.getElementById('label');
let idx = 0;
let timer = null;
const margin = {{left: 56, right: 18, top: 24, bottom: 46}};
function sx(x) {{ return margin.left + x / 36.0 * (canvas.width - margin.left - margin.right); }}
function sy(y) {{ return margin.top + (1.05 - y) / 1.05 * (canvas.height - margin.top - margin.bottom); }}
function drawAxes() {{
  ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(sx(0), sy(0)); ctx.lineTo(sx(36), sy(0)); ctx.moveTo(sx(0), sy(0)); ctx.lineTo(sx(0), sy(1.05)); ctx.stroke();
  ctx.fillStyle = '#374151'; ctx.font = '13px Arial';
  for (let x=0; x<=36; x+=6) {{ ctx.fillText(String(x), sx(x)-6, canvas.height-20); ctx.beginPath(); ctx.moveTo(sx(x), sy(0)); ctx.lineTo(sx(x), sy(1.05)); ctx.strokeStyle='#eef2f7'; ctx.stroke(); }}
  for (let y=0; y<=1.0; y+=0.2) {{ ctx.fillText(y.toFixed(1), 18, sy(y)+4); ctx.beginPath(); ctx.moveTo(sx(0), sy(y)); ctx.lineTo(sx(36), sy(y)); ctx.strokeStyle='#eef2f7'; ctx.stroke(); }}
  ctx.fillStyle = '#111827'; ctx.fillText('x (m)', canvas.width/2-18, canvas.height-8); ctx.fillText('alpha_l', 8, 18);
}}
function drawFrame(i) {{
  idx = Math.max(0, Math.min(frames.length-1, i));
  const f = frames[idx];
  ctx.clearRect(0,0,canvas.width,canvas.height);
  drawAxes();
  ctx.fillStyle = 'rgba(220,38,38,0.25)';
  for (const s of f.slugs) {{ ctx.fillRect(sx(s[0]), sy(1.0), sx(s[1])-sx(s[0]), sy(0)-sy(1.0)); }}
  ctx.strokeStyle = 'rgba(220,38,38,0.7)'; ctx.setLineDash([6,4]); ctx.beginPath(); ctx.moveTo(sx(0), sy(0.98)); ctx.lineTo(sx(36), sy(0.98)); ctx.stroke(); ctx.setLineDash([]);
  ctx.strokeStyle = '#111827'; ctx.lineWidth = 2; ctx.beginPath();
  let first = true;
  for (const seg of f.segments) {{
    const x0 = sx(seg[0]), x1 = sx(seg[1]), y = sy(seg[2]);
    if (first) {{ ctx.moveTo(x0, y); first = false; }} else {{ ctx.lineTo(x0, y); }}
    ctx.lineTo(x1, y);
  }}
  ctx.stroke();
  slider.value = idx;
  let maxL = 0, totalL = 0;
  for (const s of f.slugs) {{ const L = s[1]-s[0]; maxL = Math.max(maxL, L); totalL += L; }}
  label.textContent = `frame ${{idx+1}}/${{frames.length}} | t=${{f.t.toFixed(3)}} s | slug blocks=${{f.slugs.length}} | maxL=${{maxL.toFixed(2)}} m`;
}}
function setFrame(i) {{ drawFrame(i); }}
function step(di) {{ drawFrame(idx + di); }}
function togglePlay() {{
  const b = document.getElementById('play');
  if (timer) {{ clearInterval(timer); timer = null; b.textContent='播放'; }}
  else {{ timer = setInterval(() => drawFrame((idx+1)%frames.length), 170); b.textContent='暂停'; }}
}}
document.addEventListener('keydown', e => {{ if (e.key==='ArrowLeft') step(-1); if (e.key==='ArrowRight') step(1); }});
drawFrame(0);
</script>
</body>
</html>
"""
    path = FIG / "paper52_case4_N200_noinlet_thesis_shedfix_fast_viewer.html"
    path.write_text(html, encoding="utf-8")
    print(path)
    print(f"frames={len(frames)}")
    for target in (6.4, 10.7, 15.0, 19.3):
        j = min(range(len(frames)), key=lambda k: abs(frames[k]["t"] - target))
        f = frames[j]
        lengths = [s[1] - s[0] for s in f["slugs"]]
        print(f"target={target:.1f}s snap={f['t']:.3f}s nslug={len(f['slugs'])} maxL={(max(lengths) if lengths else 0.0):.3f}")


if __name__ == "__main__":
    main()
