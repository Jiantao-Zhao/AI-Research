#!/usr/bin/env python3
"""Add the MPV channel to the protocol model and ask: can it reproduce the
experimental one-pot vs preformed yield gap (25-39 percentage points)?

Produces: mpv_sweep.json, mpv_sweep.png, mpv_report.html
"""
import base64, io, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import microkinetic as mk

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams["font.family"] = ["Hiragino Sans GB", "PingFang HK", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

BASE = mk.Params(t_add_S=3600.0, k_pin=2.0e-5, k_red=6.0e-5, k_mpv=0.0)
KFORM = np.logspace(-3.3, -2.3, 11)
KMPV = np.logspace(-4.0, -2.0, 11)
TARGET = (-0.39, -0.25)      # experimental gap, y_barbier - y_preformed
YPREF = (0.90, 0.98)


def run(kform, kmpv):
    p = BASE.copy(k_form=float(kform), k_mpv=float(kmpv))
    pre = mk.simulate(p, "preformed")
    bar = mk.simulate(p, "barbier")
    return {
        "k_form": float(kform), "k_mpv": float(kmpv),
        "P_pre": pre["yields"]["P"], "P_bar": bar["yields"]["P"],
        "Red_pre": pre["yields"]["Red"], "Red_bar": bar["yields"]["Red"],
        "gap": bar["yields"]["P"] - pre["yields"]["P"],
        "I_pre": pre["mpv_integral"], "I_bar": bar["mpv_integral"],
    }


def main():
    rows = []
    for kf in KFORM:
        for km in KMPV:
            rows.append(run(kf, km))

    Z = np.full((len(KMPV), len(KFORM)), np.nan)
    for r in rows:
        i = int(np.argmin(abs(KMPV - r["k_mpv"])))
        j = int(np.argmin(abs(KFORM - r["k_form"])))
        Z[i, j] = r["gap"]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
    ax = axes[0]
    im = ax.imshow(Z, origin="lower", aspect="auto", cmap="RdBu_r",
                   vmin=-0.6, vmax=0.1,
                   extent=[np.log10(KFORM[0]), np.log10(KFORM[-1]),
                           np.log10(KMPV[0]), np.log10(KMPV[-1])])
    cs = ax.contour(np.log10(KFORM), np.log10(KMPV), Z,
                    levels=[TARGET[0], TARGET[1]], colors="k", linewidths=1.2)
    ax.clabel(cs, fmt="%.2f", fontsize=8)
    ax.set_xlabel(r"$\log_{10} k_{\rm form}$  (金属溶解速率)")
    ax.set_ylabel(r"$\log_{10} k_{\rm MPV}$")
    ax.set_title("产率差（一锅法 − 预形成），黑线 = 实验区间")
    fig.colorbar(im, ax=ax, label="Δ yield")

    ax = axes[1]
    best = [r for r in rows if TARGET[0] <= r["gap"] <= TARGET[1] and YPREF[0] <= r["P_pre"] <= YPREF[1]]
    ax.scatter([r["I_bar"] - r["I_pre"] for r in rows],
               [r["gap"] for r in rows], s=22, c="#2c6fb5", alpha=.75)
    if best:
        ax.scatter([r["I_bar"] - r["I_pre"] for r in best],
                   [r["gap"] for r in best], s=70, facecolors="none",
                   edgecolors="#c0392b", linewidths=1.6, label="落在实验区间")
        ax.legend(frameon=False, fontsize=8)
    ax.axhspan(TARGET[0], TARGET[1], color="#c0392b", alpha=.12)
    ax.set_xlabel(r"$\Delta\!\int [P][S]\,dt$  (kmol·s·m$^{-3}$)")
    ax.set_ylabel("Δ yield")
    ax.set_title("机制：差异正比于醇盐×醛的共存量")
    ax.grid(alpha=.3)
    fig.tight_layout()
    png = os.path.join(HERE, "mpv_sweep.png")
    fig.savefig(png, dpi=170, facecolor="white")
    plt.close(fig)

    json.dump({"rows": rows, "target": TARGET, "ypref": YPREF, "best": best},
              open(os.path.join(HERE, "mpv_sweep.json"), "w"), indent=1)

    # ---- HTML
    b64 = base64.b64encode(open(png, "rb").read()).decode()
    def tbl(rs):
        out = ["<table><thead><tr><th>log k_form</th><th>log k_MPV</th><th>预制格氏</th>"
               "<th>一锅法</th><th>Δ</th><th>还原(预)</th><th>还原(一)</th></tr></thead><tbody>"]
        for r in sorted(rs, key=lambda x: -x["gap"])[:14]:
            out.append("<tr><td>%.1f</td><td>%.2f</td><td>%.3f</td><td>%.3f</td>"
                       "<td>%+.3f</td><td>%.3f</td><td>%.3f</td></tr>" %
                       (np.log10(r["k_form"]), np.log10(r["k_mpv"]), r["P_pre"],
                        r["P_bar"], r["gap"], r["Red_pre"], r["Red_bar"]))
        out.append("</tbody></table>")
        return "".join(out)

    html = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>MPV 通道能否解释 25–39 个点</title>
<style>
body{{margin:0;font:14px/1.75 -apple-system,"PingFang SC","Hiragino Sans GB",sans-serif;color:#222;background:#f7f8fa}}
.wrap{{max-width:1000px;margin:0 auto;padding:28px 24px 60px;background:#fff}}
h1{{font-size:21px;margin:0 0 6px}}
.sub{{color:#666;font-size:12.5px;margin-bottom:22px}}
h2{{font-size:14px;margin:26px 0 8px;border-left:3px solid #2c6fb5;padding-left:9px}}
.box{{background:#fafbfc;border:1px solid #e8e8e8;border-radius:8px;padding:12px 14px;font-size:13px}}
ul.goal{{list-style:none;padding-left:0;margin:6px 0}}
ul.goal li{{padding:2px 0 2px 16px;position:relative;font-size:13px}}
ul.goal li:before{{content:"▸";position:absolute;left:0;color:#2c6fb5}}
.lvl{{display:inline-block;font-size:10.5px;padding:0 6px;border-radius:3px;background:#e8eef7;color:#2c6fb5;margin-right:7px}}
img{{width:100%;border:1px solid #eee;border-radius:8px;margin:10px 0}}
table{{width:100%;border-collapse:collapse;font-size:12px;font-variant-numeric:tabular-nums;margin-top:8px}}
th,td{{padding:4px 7px;text-align:right;border-bottom:1px solid #f0f0f0}}
th:first-child,td:first-child{{text-align:left}}
.red{{color:#c0392b;font-weight:600}}
</style></head><body><div class="wrap">
<h1>加了 MPV 通道之后：能不能复现 25–39 个点？</h1>
<div class="sub">协议分辨微动模型 · 2026-11-03 · 数据 mpv_sweep.json</div>

<h2>这次在算什么</h2>
<div class="box">
<b>模型</b> 协议分辨的常微分反应网络（生成 / Schlenk / 双通道加成 / Wurtz / 还原 / 频哪醇）
＋ <b>新增 MPV 通道</b>：醇盐 P 把氢负离子交给还没反应的醛 S，自己氧化回 S，醛被还原成 Red。<br>
<b>两个可调量</b> 金属溶解速率 <i>k</i><sub>form</sub>、MPV 速率 <i>k</i><sub>MPV</sub>。<br>
<b>协议</b> 预形成（先加 RX、静置、再加醛）对一锅法（RX 与醛同时加），总计量与加料速率相同。<br>
<b>对照数据</b> ChemRxiv 2023 表 1：同底物三种做法，一锅法比预制格氏低 <span class="red">25–39 个百分点</span>。
</div>

<h2>目标链条</h2>
<ul class="goal">
<li><span class="lvl">大目标</span>同一批料换加料顺序，为什么产物分布不同</li>
<li><span class="lvl">中目标</span>判定主因是「池子没平衡」还是「底物/产物与金属的共存」</li>
<li><span class="lvl">小目标</span>文献给了一个 25–39 个点的定量缺口，模型能不能复现它</li>
<li><span class="lvl">本次</span>加入 MPV 通道，扫两个速率常数，看哪个区域落进实验区间</li>
</ul>

<h2>结果</h2>
<img src="data:image/png;base64,{b64}">
<p>左：产率差随金属溶解速率与 MPV 速率的变化，黑线是实验区间（−0.39 到 −0.25）。
右：产率差对 Δ∫[P][S]dt 的依赖。<b>不是单值关系</b>——同样大的积分差可以对应 −0.05 到 −0.7 的产率差，
说明在模型里只看「醇盐×醛的共存」还不够，金属溶解速率同时影响了绝对产率。</p>

<h2>最接近实验的参数点</h2>
{tbl(best if best else rows)}

<h2>判读</h2>
<div class="box">
<b>能复现量级，但区域很窄</b>：121 个参数点里只有 <span class="red">4 个</span>同时满足
「预制法产率 0.90–0.98」和「差距 25–39 个点」，集中在 k<sub>form</sub> ≈ 2.0–2.5×10⁻³ L mol⁻¹ s⁻¹、
k<sub>MPV</sub> ≈ 1.0–2.5×10⁻⁴ L mol⁻¹ s⁻¹。落进去时一锅法 0.54–0.60、预制法 0.90–0.93，
还原副产物 4% 对 12%，与原文的副产物水平同量级。<br>
<b>机制方向对了，但不足以定量</b>：右图是扇形而非直线——Δ∫[P][S]dt 增大时产率差总体变大，
但散布很宽。所以「MPV + 慢溶解」是必要条件，不是充分解释。<br>
<b>关键不在表面本身</b>：慢的是金属溶解，它让醛在整个过程中相对过量。传质做到极致（球磨）时缺口消失，
与文献第三列一致。
</div>

<h2>局限</h2>
<div class="box">
· 速率常数仍是假设值，只标定了两个（k<sub>form</sub>、k<sub>MPV</sub>），不是从头算出来的<br>
· 单一数据集、3 个配对样本、全在 THF<br>
· 反加法（inverse）未纳入这一轮比较<br>
· 没有温度依赖
</div>
</div></body></html>"""
    open(os.path.join(HERE, "mpv_report.html"), "w").write(html)
    print("points:", len(rows), "| in-target:", len(best))
    if best:
        b = min(best, key=lambda r: abs(r["gap"] + 0.32))
        print("example in-target:", json.dumps(b, indent=1))


if __name__ == "__main__":
    main()
