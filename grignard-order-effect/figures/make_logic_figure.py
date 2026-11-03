#!/usr/bin/env python3
"""A single explanatory figure for the whole research plan."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))

plt.rcParams["font.family"] = ["Hiragino Sans GB", "PingFang HK", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"

RX_C = "#2c6fb5"
S_C = "#c0392b"
H1_C = "#2e7d32"
H2_C = "#d9822b"
GREY = "#555555"
LIGHT = "#f2f2f2"


def box(ax, x, y, w, h, text, fc="white", ec=GREY, fs=11, lw=1.2, bold=False, tc="black"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.03",
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, zorder=3, fontweight="bold" if bold else "normal", linespacing=1.6)


def arrow(ax, p0, p1, color=GREY, lw=1.6, style="-|>"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=14,
                                 color=color, lw=lw, zorder=1))


fig = plt.figure(figsize=(16.8, 11.6))
fig.patch.set_facecolor("white")

# =============================================================== Panel A
axA = fig.add_axes([0.075, 0.575, 0.395, 0.285])
axA.set_xlim(-3, 170)
axA.set_ylim(-0.05, 3.05)
axA.axis("off")
axA.set_title("① 同一锅料，只改“底物什么时候进来”", fontsize=13.5, loc="left", pad=10)

rows = [("预形成法\npreformed", 2.45, (0, 60), (75, 135)),
        ("一锅法\nbarbier", 1.45, (0, 60), (0, 60)),
        ("反加法\ninverse", 0.45, (0, 60), (0, 90))]

for label, y, (r0, r1), (s0, s1) in rows:
    axA.add_patch(Rectangle((s0, y - 0.36), s1 - s0, 0.74, fc=S_C, alpha=.10, ec="none", zorder=0))
    axA.add_patch(Rectangle((r0, y + 0.10), r1 - r0, 0.20, fc=RX_C, ec="none", zorder=2))
    axA.add_patch(Rectangle((s0, y - 0.30), s1 - s0, 0.20, fc=S_C, ec="none", zorder=2))
    axA.text(-4, y - 0.02, label, ha="right", va="center", fontsize=10, linespacing=1.5)
    axA.text((r0 + r1) / 2, y + 0.20, "RX", ha="center", va="center", color="white", fontsize=8.5)
    axA.text((s0 + s1) / 2, y - 0.20, "底物", ha="center", va="center", color="white", fontsize=8.5)

axA.annotate("", xy=(0, -0.02), xytext=(162, -0.02),
             arrowprops=dict(arrowstyle="->", color="black", lw=1))
axA.text(162, -0.12, "时间 / min", ha="right", va="top", fontsize=9.5)
axA.text(88, 2.98, "淡红 = 底物与金属共处的时间窗", fontsize=9, color=S_C, ha="center")

# =============================================================== Panel B
axB = fig.add_axes([0.525, 0.555, 0.45, 0.325])
axB.set_xlim(0, 10); axB.set_ylim(0, 10); axB.axis("off")
axB.set_title("② 顺序只可能通过两条路径留下痕迹", fontsize=13.5, loc="left", pad=10)

box(axB, 0.15, 4.4, 1.75, 1.3, "加料\n顺序", fc=LIGHT, bold=True, fs=11)
arrow(axB, (1.95, 5.55), (2.75, 6.95), color=H1_C)
arrow(axB, (1.95, 4.55), (2.75, 3.15), color=H2_C)

box(axB, 2.8, 6.25, 3.05, 1.75,
    "路径 A：池子来不及平衡\n2 RMgX $\\rightleftharpoons$ R$_2$Mg + MgX$_2$\n"
    "→ 组成 $x_i$ 与平衡值不同", fc="#eaf3ea", ec=H1_C, fs=9.8)
box(axB, 2.8, 1.75, 3.05, 1.75,
    "路径 B：底物占住金属表面\n底物吸附 → 自由位点 $\\theta$ 下降\n"
    "→ 生成变慢、副反应变多", fc="#fdf0e3", ec=H2_C, fs=9.8)

box(axB, 6.55, 6.45, 2.75, 1.35, "速率 $= \\sum_i x_i\\,k_i$\n少数组分也可能主导",
    fc="white", ec=H1_C, fs=10)
box(axB, 6.55, 1.95, 2.75, 1.35, "还原 / 频哪醇 / Wurtz\n与加成竞争",
    fc="white", ec=H2_C, fs=10)
arrow(axB, (5.9, 7.12), (6.52, 7.12), color=H1_C)
arrow(axB, (5.9, 2.62), (6.52, 2.62), color=H2_C)

box(axB, 2.55, 4.35, 5.35, 0.95, "产物分布 / 产率", fc="#333333", tc="white", bold=True, fs=12)
arrow(axB, (7.9, 6.4), (6.9, 5.35), color=H1_C)
arrow(axB, (7.9, 3.4), (6.9, 4.35), color=H2_C)

axB.text(2.8, 8.35, "只有寿命 ≳ 加料时间才留得下", color=H1_C, fontsize=9.5)
axB.text(2.8, 0.85, "任何条件下都可能发生", color=H2_C, fontsize=9.5)

# =============================================================== Panel C
axC = fig.add_axes([0.085, 0.105, 0.375, 0.275])
axC.set_xlim(-5.7, 4.7)
axC.set_ylim(-1.0, 1.35)
axC.set_yticks([])
axC.spines[["left", "right", "top"]].set_visible(False)
axC.set_title("③ 判据：池子自愈得快，还是慢？", fontsize=13.5, loc="left", pad=10)

lo, hi = np.log10(600), np.log10(7200)
axC.add_patch(Rectangle((lo, -0.62), hi - lo, 1.24, fc="#d8d8ef", ec="none", zorder=0))
axC.text((lo + hi) / 2, 0.40, "典型加料\n10 min – 2 h", ha="center", va="bottom",
         fontsize=8.8, color="#3b3b8f")
axC.axvline(np.log10(3600), color="#3b3b8f", ls="--", lw=1.3, zorder=1)
axC.text(np.log10(3600), 1.12, "判据线：τ = 加料时间", ha="center", va="bottom",
         fontsize=8.8, color="#3b3b8f")

axC.axhline(0, color="black", lw=1.2, xmin=0.02, xmax=0.99)

data = [("Me", 3.1e-4, "0.3 ms", 0.30), ("Et", 1.2e-5, "12 µs", -0.34),
        ("i-Pr", 1.8e-2, "18 ms", 0.30)]
for name, t, txt, dy in data:
    x = np.log10(t)
    axC.plot([x], [0], "o", ms=9, color=H1_C, zorder=4, mec="white", mew=1.2)
    axC.annotate(f"{name}\n{txt}", xy=(x, 0), xytext=(x, dy), ha="center",
                 va="center", fontsize=8.8, color=H1_C,
                 arrowprops=dict(arrowstyle="-", color=H1_C, lw=.7))

# t-Bu with its error bar
xt = np.log10(84)
xlo, xhi = np.log10(0.66), np.log10(1.06e4)
axC.plot([xlo, xhi], [0, 0], color="#b03a2e", lw=5, alpha=.35,
         solid_capstyle="butt", zorder=3)
axC.plot([xt], [0], "o", ms=11, color="#b03a2e", zorder=5, mec="white", mew=1.5)
axC.annotate("t-Bu\n84 s", xy=(xt, 0.05), xytext=(xt, 0.62), ha="center", va="center",
             fontsize=9.5, color="#b03a2e", fontweight="bold",
             arrowprops=dict(arrowstyle="-", color="#b03a2e", lw=.9))
axC.text(np.mean([xlo, xhi]), -0.30, "±12 kJ/mol 的误差棒\n正好跨过判据线",
         ha="center", va="center", fontsize=8.5, color="#b03a2e")

axC.text(-5.6, 0.92, "快 5–8 个数量级 → 记忆不可能留下", fontsize=9, color=H1_C)
axC.set_xlabel("log₁₀ ( 池子自愈时间 τ / s )", fontsize=9.5, labelpad=2)
axC.set_xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4])
axC.set_xticklabels(["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "1", "10", "10²", "10³", "10⁴"],
                    fontsize=8.5)

# =============================================================== Panel D
axD = fig.add_axes([0.545, 0.075, 0.42, 0.315])
axD.set_xlim(0, 10); axD.set_ylim(0, 10); axD.axis("off")
axD.set_title("④ 怎么判定：三个问题，依次排除", fontsize=13.5, loc="left", pad=10)

steps = [
    (7.9, "Q1  扣掉“底物 × 金属”的接触时间积分了吗？", "没扣 → 看到的可能只是平凡的接触时间差", "#7f8c8d"),
    (5.1, "Q2  换成 t-Bu（慢池子），效应变大了吗？", "变大 → 路径 A（物种分布）成立", H1_C),
    (2.3, "Q3  加竞争吸附剂（只占表面不反应），效应消失了吗？", "消失 → 路径 B（表面）成立", H2_C),
]
for y, head, tail, col in steps:
    axD.add_patch(FancyBboxPatch((0.1, y), 9.6, 1.5, boxstyle="round,pad=0.02,rounding_size=0.06",
                                 fc="white", ec=col, lw=1.6))
    axD.text(0.4, y + 1.06, head, fontsize=10, va="center", fontweight="bold", color="#222222")
    axD.text(0.4, y + 0.44, tail, fontsize=9, va="center", color=col)
for y0 in (7.7, 4.9):
    arrow(axD, (5.0, y0), (5.0, y0 - 1.05), color="#aaaaaa", lw=1.4)

axD.text(0.1, 1.55, "三项都排除 → 效应属于“混合 / 传质”，与化学无关",
         fontsize=9, color=GREY, style="italic")

fig.text(0.035, 0.975, "加料顺序为什么会影响产率：一张图",
         fontsize=19, fontweight="bold", ha="left", va="top")
fig.text(0.035, 0.938,
         "平衡态下 Schlenk 分布只由组成决定，与“怎么倒进去”无关。所以顺序效应只能来自：池子没来得及平衡，或底物改变了金属表面。",
         fontsize=10.5, color="#444444", ha="left", va="top")

fig.savefig(os.path.join(HERE, "fig0_logic.png"), dpi=190, facecolor="white", bbox_inches="tight")
print("ok")
