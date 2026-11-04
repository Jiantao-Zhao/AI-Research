"""读 WHAM 结果，画 PMF，取势垒。"""
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["Hiragino Sans GB", "PingFang HK", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "run/pmf.json"))
r = np.array(d["r"]); pmf = np.array(d["pmf"]); nb = np.array(d["counts"])
m = nb > 0
r, pmf = r[m], pmf[m]

# 参考态：r >= 3.4 的最小值
far = r >= 3.4
refmask = far if far.sum() >= 3 else np.ones_like(r, bool)
ref = pmf[refmask].min()
p = pmf - ref
i_max = int(np.argmax(np.where(r < 3.2, p, -1e9)))
barrier = float(p[i_max])
print(f"参考态 (r={r[refmask][np.argmin(pmf[refmask])]:.2f}) 取 0")
print(f"势垒 {barrier:.1f} kcal/mol at r = {r[i_max]:.2f} Å")
print(f"文献 8.7 ± 0.8")

fig, ax = plt.subplots(figsize=(6.4, 4.0))
ax.plot(r, p, "o-", color="#c0392b", ms=3.5, lw=1.4)
ax.axhline(8.7, color="#2c6fb5", ls="--", lw=1.2)
ax.text(r.min()+0.05, 9.0, "文献 ΔA‡ = 8.7", fontsize=8.5, color="#2c6fb5")
ax.annotate(f"{barrier:.1f}", xy=(r[i_max], barrier), xytext=(6, 10),
            textcoords="offset points", fontsize=10, color="#c0392b", fontweight="bold")
ax.set_xlabel("C1–C20 距离 / Å"); ax.set_ylabel("PMF / kcal mol⁻¹")
ax.set_title("显式 THF 溶剂中的加成自由能面（伞形采样 + WHAM）", fontsize=10.5)
ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig("pmf_solvent.png", dpi=170, facecolor="white")
print("saved pmf_solvent.png")
json.dump({"barrier": barrier, "r_at_barrier": float(r[i_max]),
           "n_windows": int(len(r))}, open("pmf_summary.json", "w"), indent=1)
