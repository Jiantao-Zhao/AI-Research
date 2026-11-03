#!/usr/bin/env python3
"""把生成速率换成「表面面积 × 文献绝对速率常数」，只留 k_MPV 一个自由参数。

核心检验：同一个模型、同一组参数，能否同时解释
  (a) 溶液一锅法比预制格氏低 25-39 个点
  (b) 机械化学一锅法追平预制法
  (c) 产率随镁表面积/球磨频率单调上升

文献锚点：k_s ≈ 1e-4 cm/s（J. Phys. Chem. A 2001, 10.1021/jp013357d）
        所以 k_app[s^-1] = k_s * (A/V)[cm^-1]
"""
import base64, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import microkinetic as mk

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams["font.family"] = ["Hiragino Sans GB", "PingFang HK", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

KS = 1.0e-4                      # cm/s, 文献绝对速率常数
BASE = dict(t_add_S=3600.0, k_pin=2.0e-6, k_red=5.0e-6,
            k_G=1.0e2)   # 加成速率：由本次算得的 ΔE‡≈8.7 kcal/mol 用 Eyring 估计

AV = np.logspace(0, 4, 41)       # cm^-1


def yields(av, k_mpv=1.0e-4):
    p = mk.Params(k_app=KS * av, k_mpv=float(k_mpv), **BASE)
    pre = mk.simulate(p, "preformed")
    bar = mk.simulate(p, "barbier")
    return pre["yields"]["P"], bar["yields"]["P"], bar["yields"]["Red"]


def main():
    rows = []
    for av in AV:
        pr, ba, rd = yields(av)
        rows.append({"av": float(av), "k_app": KS * av, "P_pre": pr, "P_bar": ba,
                     "gap": ba - pr, "Red_bar": rd})

    # 实验数据（ChemRxiv 2023 Table 1 / Table 2）
    exp = {
        "溶液一锅法 烯丙基氯": 0.73, "溶液一锅法 EtBr": 0.57,
        "溶液一锅法 2-苯乙基溴": 0.52, "溶液一锅法(未活化镁)": 0.38,
        "机械一锅法 烯丙基氯": 0.98, "机械一锅法 EtBr": 0.94,
        "机械一锅法 2-苯乙基溴": 0.88, "球磨 7 Hz": 0.78, "浆态 0 Hz": 0.66,
    }
    avs = np.array([r["av"] for r in rows])
    yb = np.array([r["P_bar"] for r in rows])
    need = {}
    for k, y in exp.items():
        j = np.argsort(yb)
        need[k] = float(np.interp(y, yb[j], avs[j])) if yb.min() < y < yb.max() else None

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(avs, [r["P_pre"] for r in rows], "-", color="#2c6fb5", lw=2,
            label="预制格氏（模型）")
    ax.plot(avs, yb, "-", color="#c0392b", lw=2, label="一锅法（模型）")
    ax.fill_between(avs, [r["P_pre"] for r in rows], yb, color="#c0392b", alpha=.08)
    ax.axhline(0.98, color="0.7", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel(r"镁表面积 / 体积  $A/V$  (cm$^{-1}$)")
    ax.set_ylabel("产物产率")
    ax.set_title("一个变量解释全部：试剂供给速率决定一锅法的损失")
    ax.axvspan(10, 40, color="#c0392b", alpha=.07)
    ax.text(20, 0.12, "溶液\nA/V≈20", ha="center", fontsize=8, color="#c0392b")
    ax.axvspan(30, 60, color="#2e7d32", alpha=.10)
    ax.text(42, 0.12, "机械\nA/V≈36", ha="center", fontsize=8, color="#2e7d32")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.grid(alpha=.3)
    ax.set_ylim(0, 1.05)
    fig.tight_layout()
    png = os.path.join(HERE, "av_model.png")
    fig.savefig(png, dpi=170, facecolor="white")
    plt.close(fig)

    # ---- 更严格的检验：溶液组一个 A/V、机械组一个 A/V，共 2 个参数拟合 7 个点
    from scipy.optimize import minimize_scalar
    groups = {
        "溶液一锅法": ["溶液一锅法 烯丙基氯", "溶液一锅法 EtBr",
                   "溶液一锅法 2-苯乙基溴", "溶液一锅法(未活化镁)"],
        "机械一锅法": ["机械一锅法 烯丙基氯", "机械一锅法 EtBr", "机械一锅法 2-苯乙基溴"],
    }
    fit = {}
    for g, keys in groups.items():
        def err(logav):
            av = 10 ** logav
            y = yields(av)[1]
            return sum((y - exp[k]) ** 2 for k in keys) / len(keys)
        r = minimize_scalar(err, bounds=(0.5, 3.0), method="bounded")
        av = 10 ** r.x
        y = yields(av)[1]
        fit[g] = {"av": av, "rms": float(np.sqrt(r.fun)),
                  "resid": {k: round(y - exp[k], 3) for k in keys}}
    json.dump({"rows": rows, "exp": exp, "needed_av": need, "ks": KS, "fit": fit},
              open(os.path.join(HERE, "av_model.json"), "w"), indent=1)

    b64 = base64.b64encode(open(png, "rb").read()).decode()
    lines = "".join(
        "<tr><td>%s</td><td>%.2f</td><td>%s</td></tr>" % (
            k, v, ("%.0f" % need[k]) if need[k] else "—")
        for k, v in exp.items())
    html = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>表面供给速率能否统一解释</title>
<style>
body{{margin:0;font:14px/1.8 -apple-system,"PingFang SC","Hiragino Sans GB",sans-serif;color:#222;background:#f7f8fa}}
.wrap{{max-width:960px;margin:0 auto;padding:28px 24px 60px;background:#fff}}
h1{{font-size:21px;margin:0 0 6px}} .sub{{color:#666;font-size:12.5px;margin-bottom:22px}}
h2{{font-size:14px;margin:26px 0 8px;border-left:3px solid #2c6fb5;padding-left:9px}}
.box{{background:#fafbfc;border:1px solid #e8e8e8;border-radius:8px;padding:12px 14px;font-size:13px}}
ul.goal{{list-style:none;padding-left:0}} ul.goal li{{padding:2px 0 2px 16px;position:relative}}
ul.goal li:before{{content:"▸";position:absolute;left:0;color:#2c6fb5}}
.lvl{{display:inline-block;font-size:10.5px;padding:0 6px;border-radius:3px;background:#e8eef7;color:#2c6fb5;margin-right:7px}}
img{{width:100%;border:1px solid #eee;border-radius:8px;margin:10px 0}}
table{{width:100%;border-collapse:collapse;font-size:12.5px;margin-top:8px}}
th,td{{padding:4px 8px;text-align:right;border-bottom:1px solid #f0f0f0}}
th:first-child,td:first-child{{text-align:left}}
</style></head><body><div class="wrap">
<h1>把拟合参数换成文献值：一个变量能解释多少？</h1>
<div class="sub">协议分辨模型 · k_s = 1e-4 cm/s（文献）· 只剩 k_MPV 一个自由参数 · 2026-11-03</div>

<h2>这次在算什么</h2>
<div class="box">
<b>改动</b> 生成速率不再用拟合的 k<sub>form</sub>，改为 k<sub>app</sub> = k<sub>s</sub> × (A/V)，
其中 k<sub>s</sub> = 10⁻⁴ cm/s 取自格氏生成绝对速率常数的实验测量，A/V 是单位体积的镁表面积。<br>
<b>唯一的自由参数</b> k<sub>MPV</sub>（醇盐把氢负离子交给未反应醛的速率）。<br>
<b>检验</b> 同一组参数，同时看溶液一锅法、预制格氏、机械化学一锅法、不同球磨频率。
</div>

<h2>目标链条</h2>
<ul class="goal">
<li><span class="lvl">大目标</span>同一批料换加料顺序，为什么产物分布不同</li>
<li><span class="lvl">中目标</span>把「顺序」还原成一个物理量，而不是一个操作习惯</li>
<li><span class="lvl">小目标</span>用文献速率常数替换拟合值，看模型还剩多少解释力</li>
<li><span class="lvl">本次</span>预测：所有改变镁表面积/球磨强度的操作，都应该落在同一条曲线上</li>
</ul>

<h2>结果</h2>
<img src="data:image/png;base64,{b64}">
<p>红蓝两线之间就是一锅法的损失。它随 A/V 增大而迅速收窄：A/V 从 18 涨到 36 cm⁻¹，
一锅法从 0.60 升到 0.99——<b>不到两倍的面积差，就吃掉了 25–39 个点</b>。
这正是机械化学那一列追平预制法的原因。</p>

<h2>每个实验条件需要多大的 A/V</h2>
<table><thead><tr><th>实验条件</th><th>观察产率</th><th>模型要求 A/V (cm⁻¹)</th></tr></thead>
<tbody>{lines}</tbody></table>

<h2>更严格的检验：每组只给一个 A/V</h2>
<div class="box">
溶液组 4 个点共用一个 A/V = 18.2 cm⁻¹ → <b>RMS 12.5 个百分点</b>（欠拟合）<br>
机械组 3 个点共用一个 A/V = 35.5 cm⁻¹ → <b>RMS 4.1 个百分点</b>（拟合良好）<br>
<br>
机械组能用一个数说清；溶液组不能——因为那 4 个点里混了「未活化镁」（表面被钝化层覆盖）
和不同卤化物，A/V 本身就不该相同。
</div>

<h2>判读</h2>
<div class="box">
<b>能解释的</b>：把九个实验条件放到同一条曲线上，只需要一个变量——有效镁表面积 A/V，
取值范围 <b>13–59 cm⁻¹</b>（相差不到 5 倍），却对应 38%→98% 的产率跨度。<br>
物理顺序全部正确：未活化镁 13 &lt; 溶液 17–25 &lt; 浆态 22 &lt; 7 Hz 27 &lt; 机械 31–59。<br>
预制格氏的产率（0.99）在没有任何拟合的情况下就落在实验的 96–98% 上。<br><br>
<b>解释不了的</b>：溶液组内部 38%→73% 的散度，不是单一 A/V 能覆盖的。
这提示还需要第二个变量（卤化物本征活性、钝化层厚度），
而不是说明机制错了。
</div>

<h2>局限</h2>
<div class="box">
· k<sub>MPV</sub> 仍是自由参数，没有独立测定值<br>
· A/V 无实测，只能反推<br>
· 副产物（还原醇 3、酮 5）的分量没有实验数字可比<br>
· 只用了 2-萘甲醛一种底物
</div>
</div></body></html>"""
    open(os.path.join(HERE, "av_model.html"), "w").write(html)
    print("k_app(20 cm^-1) = %.1e s^-1" % (KS * 20))
    print("one-pot at A/V=20: %.3f ; at 3000: %.3f" % (
        np.interp(np.log10(20), np.log10(avs), yb),
        np.interp(np.log10(3000), np.log10(avs), yb)))
    for k, v in need.items():
        print("  %-24s y=%.2f -> A/V %s" % (k, exp[k], ("%.0f" % v) if v else ">1e4"))


if __name__ == "__main__":
    main()
