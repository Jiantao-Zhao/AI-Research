#!/usr/bin/env python3
"""Plot the zero-shot MACE-POLAR barrier profiles against the JACS 2020 values."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams["font.family"] = ["Hiragino Sans GB", "PingFang HK", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

d = json.load(open(os.path.join(HERE, "scan_results.json")))
meta = {
    "Agem_like": dict(label="Agem-like  CH$_3$MgCl(THF)(ACA)", color="#c0392b",
                      lit=8.7, lit_err=0.8, lit_txt="8.7 ± 0.8"),
    "Bgem_like": dict(label="Bgem-like  Mg(CH$_3$)$_2$(THF)(ACA)", color="#2c6fb5",
                      lit=6.5, lit_err=0.7, lit_txt="6.5 ± 0.7"),
}

fig, ax = plt.subplots(figsize=(6.4, 4.4))
for key, m in meta.items():
    prof = d[key]["profile"]
    r = [p["r"] for p in prof]
    e = [p["dE_kcal"] for p in prof]
    ax.plot(r, e, "o-", ms=4, lw=1.8, color=m["color"], label=m["label"])
    peak = d[key]["max_dE_kcal"]
    ax.plot([d[key]["r_at_max"]], [peak], "*", ms=16, color=m["color"], zorder=5)
    dx, dy = (10, -4) if key == "Agem_like" else (-26, -6)
    ax.annotate(f"{peak:.1f}", xy=(d[key]["r_at_max"], peak), xytext=(dx, dy),
                textcoords="offset points", fontsize=9, color=m["color"],
                fontweight="bold")
    ax.axhline(m["lit"], color=m["color"], ls="--", lw=1.2, alpha=.6)
    ax.text(4.32, m["lit"] + 0.4, f"JACS 2020  ΔA‡ = {m['lit_txt']}",
            fontsize=8.5, color=m["color"], ha="left")

ax.axhline(0, color="0.6", lw=.8)
ax.set_xlabel(r"C(CH$_3$)$\cdots$C(C=O) 距离 / $\mathrm{\AA}$")
ax.set_ylabel("相对能量 / kcal mol$^{-1}$")
ax.set_title("MACE-POLAR-1 零样本：格氏试剂 + 乙醛的加成势垒", fontsize=11)
ax.invert_xaxis()
ax.grid(alpha=.3)
ax.legend(fontsize=8.5, frameon=False, loc="lower left")
ax.set_ylim(-25, 12)



fig.tight_layout()
out = os.path.join(HERE, "barrier_bench.png")
fig.savefig(out, dpi=200, facecolor="white")
print("saved", out)
