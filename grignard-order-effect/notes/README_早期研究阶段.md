# grignard-order-effect

**加料顺序 → 有机镁微观状态分布 → 选择性**：文献、判据、计算方案与可行性验证。

一句话：在平衡态下加料顺序不该留下痕迹，所以顺序效应只能来自**动力学记忆**或**第二相（Mg 表面）**；这两条可以分别证伪。用已发表的 Schlenk 势垒做数量级估算，Me/Et/i-Pr 的记忆在毫秒内被抹平，只有 t-Bu 落在"加料时间尺度"的窗口里——这就是这个问题真正可以被检验的地方。

---

## 一张图看懂

![logic](figures/fig0_logic.png)

① 三个方案，同一锅料、同样加料速度，只差"底物什么时候进来"；淡红是底物与金属共处的时间。
② 顺序只可能通过两条路径留下痕迹：池子来不及平衡（A），或底物占住金属表面（B）。
③ 判据是"池子自愈时间" vs "加料时间"：Me/Et/i-Pr 快 5–8 个数量级，A 不可能；t-Bu 卡在判据线上，误差棒正好跨过去。
④ 判定顺序：先扣接触时间积分 → 换 t-Bu → 加竞争吸附剂。

---

## 先读这一个

**[研究方案.md](研究方案.md)** — 完整方案：文献判读、三层问题重述、三个假设及其判别指纹、L1–L3 计算路线、实验设计、风险。

核心结论三条：

1. **空白在建模，不在“没人想到”。** 工程界（React. Chem. Eng. 2023）已证明加料速率/停留时间分布改变 Wurtz 选择性；Ashby 1980 已指出镁的制备史影响机理；有机锌领域（Acc. Chem. Res. 2026）已明确提出"动态物种分布是隐藏驱动力"。缺的是**把加料协议当作自变量写进动力学模型**。
2. **τ_Schlenk 是判据。** 用 J. Catal. 2025 的势垒（45–84 kJ/mol）做 Eyring 估算：Me 0.3 ms / Et 12 µs / i-Pr 18 ms / t-Bu **84 s**。前三个体系里物种分布记忆不可能存在；t-Bu 是唯一的窗口——但 ±12 kJ/mol 的误差把 τ 摊到 0.7 s–10⁴ s，**四个数量级的不确定性本身就是该做的计算**。
3. **零模型会吃掉大部分效应。** pilot 模型里，三种协议的还原副产物差异（4.0 / 7.3 / 9.6%）**严格正比于** `∫[S]·θ·sites dt`。所以"看到顺序效应"不构成任何物种分布机制的证据；必须先把这个平凡项回归掉。

---

## 目录

```
grignard-order-effect/
├── 研究方案.md              ← 主文档
├── lit/                     文献（7 篇开放获取 PDF + 抽取文本）
│   ├── 文献清单.md           下载状态、每篇的作用、4 篇付费墙文献说明
│   ├── extraction_report.md  对 Peltzer2017 / JACS2020 的逐条事实抽取（带行号）
│   └── *.pdf / *.txt
├── model/                   协议分辨的微观动力学模型（Python + scipy）
│   ├── microkinetic.py       反应网络、三种协议、求解
│   ├── analysis.py           参数扫描、相图
│   ├── README.md             模型说明与运行方法
│   ├── xtb/                  GFN2-xTB 方法标定（结论：半经验不可用）
│   └── figures/
└── notes/                   （预留）
```

## 三张图

| 图 | 说明 |
|---|---|
| ![fig1](model/figures/fig1_timecourses.png) | 三协议浓度-时间曲线。相同总计量、相同加料速率，差别只在"自由底物 × 自由表面"的重叠窗口。 |
| ![fig2](model/figures/fig2_phase_map.png) | 协议效应相图：`ΔYield(Barbier − preformed)` 随底物还原速率与表面吸附常数的变化。效应只在"强吸附 + 易还原"的角落显著。 |
| ![fig3](model/figures/fig3_protocol_curves.png) | 三个协议的产物产率随还原速率的分叉。 |

---

## 复现

```bash
cd model
MPLCONFIGDIR=/tmp/mpl python3 analysis.py      # ~2–3 min，重生成全部图与 sweep_results.json
python3 microkinetic.py                        # 基线三协议
```

xTB 标定（GFN2-xTB 6.7.1，`/opt/anaconda3/envs/qc/bin/xtb`）：

```bash
cd model/xtb && python3 build_solvated.py && python3 run_solvated.py
```

## 需要人工处理的一件事

`lit/文献清单.md` 里列了 4 篇付费墙文献（Chen 2012 OBC、Yamazaki 2002 JOC、Blomberg 1977 Synthesis、俄语学派 CTC 2013）。它们**不是本方案的判据性文献**——判据是 J. Catal. 2025（τ）和 JACS 2020（k_i），两篇都拿到了开放版本。如需全文，走学校图书馆代理。
