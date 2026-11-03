# model/ — 协议分辨的微观动力学模型

最小可行模型，用来回答一件事：**在什么条件下，加料顺序能够改变产物分布？**

## 文件

| 文件 | 作用 |
|---|---|
| `microkinetic.py` | 反应网络 + 三种加料协议 + 求解器；核心库 |
| `analysis.py` | 参数扫描、相图、作图 |
| `build_species.py`, `build_solvated.py` | 手搓起始几何（CH₃MgCl 体系及其 THF 溶剂化物） |
| `run_species.py`, `run_solvated.py` | GFN2-xTB 优化 + 频率，做方法标定 |
| `figures/` | 三张图 |
| `baseline.json`, `sweep_results.json`, `schlenk_*.json` | 结果 |

## 运行

```bash
cd model
MPLCONFIGDIR=/tmp/mpl python3 analysis.py     # 全部重跑，约 2–3 分钟
python3 microkinetic.py                        # 只看基线
```

xTB 部分需要 `/opt/anaconda3/envs/qc/bin/xtb`（GFN2-xTB 6.7.1）：

```bash
cd model/xtb
python3 build_solvated.py && python3 run_solvated.py
```

## 网络

```
F   RX + Mg(site)  -> RMgX          k_form
R   S  + Mg(site)  -> Red           k_red
W   RMgX + RX      -> W             k_wurtz
A1  RMgX + S       -> P             k_G
A2  R2Mg + 2 S     -> 2 P           k_D
S1  2 RMgX <-> R2Mg + MgX2          k_sf, k_sr
Pn  2 S + Mg(site) -> Pin           k_pin
```

表面用竞争 Langmuir：`θ₀ = 1/(1 + K_RX[RX] + K_S[S])`。底物因此会毒化"制造试剂的那块表面"——这是模型里唯一把进料顺序和化学连起来的机制。

## 协议（唯一差别是进料轨迹）

```
preformed : RX 0→t_feed；静置 t_age；S 于 t_feed+t_age 开始，按 t_add_S 加完
barbier   : RX 与 S 同时 0→t_feed
inverse   : S 于 t=0 一次性投入；RX 0→t_feed
```

三者的总计量、S 的加料速率、后处理时间（t_hold）完全相同。**这是"控制总计量、只变顺序"的最小对照设计。**

## 基线结果（`Params()` 默认值）

| 协议 | P | 还原 | ∫[S]θ dt |
|---|---|---|---|
| preformed | 0.939 | 0.040 | 134.5 |
| barbier | 0.910 | 0.073 | 243.6 |
| inverse | 0.885 | 0.096 | 319.9 |

还原产物严格正比于 `∫[S]·θ·sites dt`（比值 3.000×10⁻⁴）。**这是零模型的定量形式。**

## 这一层不做什么

- 速率常数是**明确的假设值**，不是拟合值。模型的价值在结构（哪些无量纲量控制效应），不在数值。
- Schlenk 平衡被当作瞬时平衡（τ = 0.6 s），因此 **H1（物种分布假说）被显式关掉了**。这正是设计意图：先把平凡解释算到尽头，再去看还剩什么。
- 不含 Mg 表面电子结构、不含聚集体的显式结构、不含局部浓度梯度。

这些留给 `研究方案.md` 的 L1–L3 三层。
