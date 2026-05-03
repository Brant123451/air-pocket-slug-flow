# K–L 2018 Case 4 — 论文对标后的调参清单

**背景**：2026-04-30 把 LASSI Fortran 代码全部对标论文公式
（`τ` 系数 1/2 → 1/8、`λ_i` 默认 `λ_g`、Bendiksen 双段+max+cos/sin 外提、
SR/SS Riemann `U_b` 修正、Moissis–Griffith W_eff 函数）。

**已知风险**：旧版（带 1/2·τ + Andritsos–Hanratty `λ_i = (1+75β)·λ_g`）
通过过度摩擦耗散，碰巧把 K-L 2018 Case 4 的 slug 频率调到
**0.453 Hz**（论文实验 0.46 Hz，1.5 % 匹配）。新代码下这种隐式
抵消被打破，需重新校准。

---

## 1. 第一步：先跑一次新代码看现状

```powershell
cd "e:\PHD Thesis\lassi"
.\build.ps1 -Compiler gfortran
.\build\lassi.exe cases\INPUT_paper52_case4.txt outputs_thesis_align
python plot_kl2018_xt.py outputs_thesis_align/KL2018_case4_*.dat
```

记录：

- `slug_count_at_outlet`：通过 z=30 m 监测点的 slug 数 / 时长 = 频率
- `outputs_thesis_align/KL2018_case4_xt_beta.png`：xt 等值线（slug 是否覆盖 10–35 m）
- `outputs_thesis_align/KL2018_case4_xt_slugcount.png`：稳态 slug 内含数

**判定**：
- 若频率 ∈ [0.4, 0.55] Hz ⇒ 自然匹配，**无需调参**。
- 若频率 < 0.3 Hz ⇒ 摩擦/扰动太弱，slug 形成不足。⇒ 调表 A。
- 若频率 > 0.6 Hz ⇒ 过度不稳定，slug 过频生成。⇒ 调表 B。

---

## 2. 调参表 A：频率偏低（slug 太少 / 过晚生成）

| 参数 | INPUT 字段 | 当前默认 | 建议值 | 物理解释 | 影响量级 |
|---|---|---|---|---|---|
| 入口 β 扰动 | `perturb_beta` | 0.02 | **0.04 ~ 0.08** | 增大入口噪声幅度 ⇒ KH 不稳定更易触发 | 🔴 大：直接控制 slug 起始位置 |
| 入口 U_l 扰动 | `perturb_Ul` | 0.10 | **0.20 ~ 0.30** | 同上 | 🔴 大 |
| 管壁粗糙度 | `rough` (代码) | 0 | **5e-5（50 μm）** | 论文 page 19 默认；增加湍流摩擦 ⇒ stratified 平衡 holdup ↑ ⇒ 更易跨越 β_init | 🟡 中 |
| Slug 起始阈值 | `beta_init` | 0.98 | **0.95 ~ 0.97** | 降低 slug 形成门槛 | 🟡 中：会提高假阳性 |
| 启用 Andritsos-Hanratty | `use_AH` (代码) | `.false.` | **改回 `.true.`** | 强界面摩擦 ⇒ 气液耦合 ⇒ 干扰放大 | 🔴🔴 大：等于回退一半改动 |
| AH 强度乘子 | `ai_factor` | 1.0 | 同上前置 | 进一步缩放 AH 公式 | 视 use_AH 而定 |

---

## 3. 调参表 B：频率偏高（slug 过频 / 数值噪声引入伪 slug）

| 参数 | INPUT 字段 | 当前默认 | 建议值 | 物理解释 |
|---|---|---|---|---|
| 入口扰动 | `perturb_beta`, `perturb_Ul` | 0.02 / 0.10 | **0.01 / 0.05** | 减小入口噪声 |
| Slug 起始阈值 | `beta_init` | 0.98 | **0.99** | 提高生成门槛，过滤数值伪 slug |
| TargetLength | `TargetLength` | 0.18 | **0.30 ~ 0.50** | 较粗网格 ⇒ 高频 KH 模式被截断 |
| Wake-effect 上限 | `W_eff` | 1.0 | **保持 1.0** | 不要开 Moissis–Griffith；否则短 slug 加速碰撞会增频 |

---

## 4. 调参表 C：形态/速度不对（频率正确但 slug 长度/速度偏差）

| 现象 | 调整项 | 方向 |
|---|---|---|
| Slug 平均长度过短 | `TargetLength` ↑ 0.18 → 0.25 | 减小数值耗散 |
| Slug 平均长度过长 | `TargetLength` ↓ 0.18 → 0.10 | 增加分裂细度 |
| Slug 跑得太快（U_nose 过大） | 改 `bendiksen_nose` 强制用 low-Fr 段 | 在 `lassi_friction.f90:107` 改 `Un = U1` 而非 `max(U1, U2)` |
| Slug 跑得太慢 | 启用 wake effect | INPUT 中 `W_eff = 1.5 ~ 2.0` |
| 短 slug 越追越快导致兼并失稳 | 关闭 wake effect | `W_eff = 1.0` |

---

## 5. 验证流程模板

```text
1. 跑 baseline (不改 INPUT)
   ├─ frequency_obs = ?
   ├─ qualitative behaviour (xt etc.)
   └─ 与 paper case4 对比 (target 0.46 Hz, 7-9 slugs steady-state)

2. 如不匹配，按表 A 或 B 单参数扫
   ├─ 一次只动 1 个参数
   ├─ 对比频率变化
   └─ 收敛到 |f_obs - 0.46| < 0.05 Hz

3. 验证泛化
   ├─ 复跑 case1, case2, case3 (USL, USG 不同)
   └─ 同一组参数应给类似实验质量的结果
```

---

## 6. 关键源码定位

| 待调代码 | 文件 | 行 | 当前值 |
|---|---|---|---|
| use_AH 默认 | `src/lassi_friction.f90` | 19 | `.false.` |
| τ 系数 | `src/lassi_friction.f90` | 76-78 | `0.125` (Darcy 1/8) |
| λ_i = λ_g 分支 | `src/lassi_friction.f90` | 72 | `lam_i = lam_g` |
| Bendiksen max(low, high Fr) | `src/lassi_friction.f90` | 107 | `Un = max(U1, U2)` |
| wake_effect 函数 | `src/lassi_friction.f90` | 117 | 新增 |
| SR `U_b` | `src/lassi_riemann.f90` | 240 | `0.5*(U_LR+U_RL)` |
| SS `U_b` | `src/lassi_riemann.f90` | 271 | `sR` |
| SSS `β_MR/U_MR` | `src/lassi_riemann.f90` | 303 | `(β_R, U_R)` |

---

## 7. 历史性能记录（修改前 vs 修改后对照基线）

| 量 | 修改前 (旧) | 修改后 (新, 待测) | K-L 2018 case 4 论文 |
|---|---|---|---|
| Slug 频率 (z=30m) | 0.453 Hz ✓ | **TBD** | 0.46 Hz |
| 稳态 slug 数 (在管内) | 7–9 | **TBD** | 7–9 |
| Slug 起始位置 | z ≈ 10 m | **TBD** | z ≈ 10 m |
| 远端最大瞬时 U | 3.0 m/s | **TBD** | – |
| 总耗时 (20 s 仿真) | ~10 min | **TBD** | – |

跑完 §1 后填进 TBD。
