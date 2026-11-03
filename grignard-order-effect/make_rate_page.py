#!/usr/bin/env python3
"""生成《速率是怎么算的》HTML：实验 → 量子化学 → 速率 → 产率 的完整链条。"""
import base64, os

HERE = os.path.dirname(os.path.abspath(__file__))


def img(path):
    with open(os.path.join(HERE, path), "rb") as f:
        return base64.b64encode(f.read()).decode()


FIG_MPV = img("bench/mpv_profile.png")
FIG_FRAMES = img("bench/frames_Agem.png")
FIG_RMS = img("model/summary_rms.png")
FIG_REGIME = img("model/regime.png")

html = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>速率是怎么算的：实验、量子化学与产率模型如何接上</title>
<style>
body{{margin:0;font:15px/1.9 -apple-system,"PingFang SC","Hiragino Sans GB",sans-serif;color:#1c1c1c;background:#f6f7f9}}
.wrap{{max-width:940px;margin:0 auto;padding:36px 28px 80px;background:#fff}}
h1{{font-size:24px;margin:0 0 6px;line-height:1.35}}
.sub{{color:#666;font-size:13px;margin-bottom:28px}}
h2{{font-size:16px;margin:34px 0 10px;border-left:3px solid #2c6fb5;padding-left:10px}}
h3{{font-size:14px;margin:20px 0 6px;color:#333}}
.box{{background:#fafbfc;border:1px solid #e8e8e8;border-radius:8px;padding:13px 16px;font-size:13.5px}}
.key{{background:#f2f6fb;border:1px solid #dde6f2}}
.warn{{background:#fff8e6;border:1px solid #f0d9a8}}
table{{width:100%;border-collapse:collapse;font-size:13px;margin-top:8px}}
th,td{{padding:5px 9px;text-align:right;border-bottom:1px solid #f0f0f0}}
th:first-child,td:first-child{{text-align:left}}
img{{width:100%;border:1px solid #eee;border-radius:8px;margin:12px 0}}
.g{{color:#2e7d32;font-weight:600}} .r{{color:#c0392b;font-weight:600}} .o{{color:#d9822b;font-weight:600}}
.eq{{background:#f7f8fa;border:1px solid #e8e8e8;border-radius:6px;padding:9px 13px;font-family:"SF Mono",Menlo,monospace;font-size:13px;margin:8px 0;overflow-x:auto}}
.flow{{display:flex;align-items:stretch;gap:8px;flex-wrap:wrap;margin:16px 0}}
.node{{flex:1;min-width:170px;background:#fff;border:1px solid #dde3ea;border-radius:8px;padding:11px 12px;font-size:12.5px}}
.node b{{display:block;font-size:13px;margin-bottom:4px;color:#2c6fb5}}
.node .v{{color:#666;font-size:11.5px}}
.arr{{align-self:center;color:#9aa7b5;font-size:18px}}
code{{background:#f0f2f5;border-radius:4px;padding:1px 5px;font-size:12.5px}}
ul{{padding-left:22px}} li{{margin:4px 0}}
.small{{font-size:12.5px;color:#666}}
</style></head><body><div class="wrap">

<h1>速率是怎么算的：实验、量子化学与产率模型如何接上</h1>
<div class="sub">Barbier 一锅法 vs 预制格氏 · 2026-11-03</div>

<h2>一句话</h2>
<div class="box key">
速率常数来自三个互不相同的来源：<b>格氏生成速率</b>取自文献实测；
<b>加成与 MPV 的速率</b>用量子化学（MACE-POLAR 机器学习势）算出势垒再换算；
<b>镁表面积</b>由几何估算。三者放进一个常微分方程网络，算出产率，再和实验对照。
</div>

<h2>一、整条链</h2>
<div class="flow">
<div class="node"><b>实验</b>产率 52–73%（一锅法）<br>90–98%（预制法）<span class="v">ChemRxiv 2023</span></div>
<div class="arr">→</div>
<div class="node"><b>几何</b>镁粉粒径 → A/V<span class="v">400 µm → 17 cm⁻¹<br>200 µm → 35 cm⁻¹</span></div>
<div class="arr">→</div>
<div class="node"><b>量子化学</b>MACE-POLAR 算势垒<span class="v">加成 8.7 kcal/mol<br>MPV 15.5 kcal/mol</span></div>
<div class="arr">→</div>
<div class="node"><b>Eyring</b>势垒 → 速率常数<span class="v">k_G ≈ 10⁶<br>k_MPV ≈ 10⁻⁴</span></div>
<div class="arr">→</div>
<div class="node"><b>微动模型</b>常微分方程网络<span class="v">产率 vs 条件</span></div>
</div>

<h2>二、量子化学到底算了什么</h2>
<div class="box">
用的是 <b>MACE-POLAR-1-M</b>——一个在 OMol25 数据集（1 亿个 ωB97M-V 级别 DFT 结构）上训练的
极化机器学习势。它不是从头算 DFT，而是 DFT 的代理；在本次工作中它扮演"能算显式溶剂、
算得起几十个原子"的角色。
</div>

<h3>算了两个反应</h3>
<div class="box">
<b>（1）羰基加成</b>：CH₃MgCl(THF) + 乙醛 → 醇盐。固定 C···C 距离，放开其余坐标做弛豫扫描，
取能量最高点。
</div>
<img src="data:image/png;base64,{FIG_FRAMES}">
<p class="small">加成路径上的四个瞬间。可以看到 C–C 在缩短，而 Mg–C 在最高点几乎没动——
这正是后来发现"鞍点没验证成功"的地方。</p>

<div class="box">
<b>（2）MPV 氢转移</b>：CH₃O–MgCl + CH₃CHO → HCHO + CH₃CH₂O–MgCl。
固定 H···C(羰基) 距离做弛豫扫描，势垒 15.5 kcal/mol。
</div>
<img src="data:image/png;base64,{FIG_MPV}">

<h2>三、从势垒到速率常数</h2>
<div class="eq">k = (k<sub>B</sub>T / h) × (1/C°) × exp( −ΔG‡ / RT )</div>
<div class="box">
指前因子 k<sub>B</sub>T/h ÷ C°（C° = 1 M）在 298 K 是 6.2×10¹² L mol⁻¹ s⁻¹。
难点在于 <b>ΔG‡ 不等于 ΔE‡</b>：量子化学算的是电子能，而速率需要自由能。
MPV 是双分子成环反应，要加一个熵校正 ΔS‡ ≈ −25 cal mol⁻¹ K⁻¹。
</div>

<table><thead><tr><th>反应</th><th>ΔE‡（算）</th><th>ΔG‡（加熵）</th><th>速率常数 k</th></tr></thead><tbody>
<tr><td>羰基加成</td><td>8.7 kcal/mol</td><td>8.7（单分子参考）</td><td class="g">≈ 2.6×10⁶ L mol⁻¹ s⁻¹</td></tr>
<tr><td>MPV 氢转移</td><td>15.5 kcal/mol</td><td>23.0</td><td class="g">≈ 9.3×10⁻⁵ L mol⁻¹ s⁻¹</td></tr>
<tr><td>格氏生成</td><td colspan="2">文献实测 k<sub>s</sub> = 1×10⁻⁴ cm/s</td><td class="g">k<sub>app</sub> = k<sub>s</sub>·A/V</td></tr>
</tbody></table>

<div class="box warn">
<b>这一步是整个模型最脆弱的地方。</b>熵校正是取的双分子成环的教科书值（−20 到 −30），
不是测出来的。模型对它的取值很敏感：ΔS‡ 取 −20 时缺口是 −0.73，取 −25 时 −0.385，
取 −30 时 −0.33。
</div>

<h2>四、从速率到产率</h2>
<div class="box">
反应网络（常微分方程）：
<div class="eq">
Mg + RX → RMgX 　　　　　　　　　 k<sub>s</sub>·(A/V)·θ
2 RMgX ⇌ R₂Mg + MgX₂ 　　　　　　 Schlenk 平衡
RMgX + S → P 　 / 　 R₂Mg + 2S → 2P 　 k<sub>G</sub>
RMgX + RX → R–R 　　　　　　　　　 Wurtz 偶联
S + Mg → 还原产物 　　　　　　　　 表面还原
2 S + Mg → 频哪醇 　　　　　　　　 表面偶联
P + S → 还原产物 　　　　　　　　　 MPV（k<sub>MPV</sub>）
</div>
加料顺序不写成参数，而是写成<b>进料时间表</b>：预制法先加 RX、静置、再加 S；
一锅法两者同时在场。θ 是镁表面的自由位点分数，由 RX 和 S 竞争吸附决定。
</div>

<h2>五、对上了多少</h2>
<img src="data:image/png;base64,{FIG_RMS}">
<img src="data:image/png;base64,{FIG_REGIME}">
<table><thead><tr><th>检验</th><th>点数</th><th>RMS</th><th>结论</th></tr></thead><tbody>
<tr><td>条件依赖（面积 / 球磨 / 镁形态）</td><td>9</td><td class="g">1.1 点</td><td>定量有效</td></tr>
<tr><td>同系一级溴 / 氯</td><td>4</td><td class="g">3.7 点</td><td>定量有效</td></tr>
<tr><td>酮类（NH₄Cl 对照）</td><td>5</td><td class="o">16 点</td><td>失效</td></tr>
<tr><td>碘化物 / 仲卤代物</td><td>11</td><td class="r">25 点</td><td>失效</td></tr>
</tbody></table>

<h2>六、量子化学在这里起到的作用，以及它的边界</h2>
<div class="box">
<b>起到的作用</b>：它提供了两个<b>实验上很难直接测量</b>的速率常数——
羰基加成的 k<sub>G</sub> 和 MPV 的 k<sub>MPV</sub>。没有这两个数，模型里就只有可调参数，
没有预测能力。它也验证了 Mg–C / Mg–Cl 这类化学在它的训练分布之内。
</div>
<div class="box warn">
<b>边界（必须说清楚）</b>：
<ul>
<li>势垒是<b>单构象扫描</b>，不是系综平均。没有多构象 Boltzmann 加权，没有分子动力学采样。</li>
<li>加成路径的"最高点"经过 Dimer + 频率检验，结果是 <b>0 个虚频</b>，最低振动模只有 0.0–0.5 cm⁻¹
——势垒顶部是一片近乎平坦的高原，数值上无法确认它是鞍点。</li>
<li>同一体系在 r = 4.40 Å 处，两个构象相差 <b>3.84 kcal/mol</b>，比文献误差棒（±0.8）还大。
所以这些势垒值本身有未收敛的成分。</li>
<li>MPV 的熵校正是假设值，不是算出来的。</li>
</ul>
</div>

<div class="box">
<b>所以正确的说法是</b>：量子化学在这里提供的是"量级正确、机制方向正确"的速率估计，
不是一个经过严格收敛的定量自由能。模型能解释条件依赖和大部分产率差，
但把它的数字当成精确值是不合适的。
</div>

<h2>七、要让它更硬，需要什么</h2>
<ul>
<li>势垒做构象系综（多构象 + Boltzmann 平均，或机器学习势的分子动力学采样）</li>
<li>MPV 的熵校正用显式溶剂自由能方法复核</li>
<li>实验上补两组数据：未反应醛的回收率、MPV 副产物 3/5 的定量</li>
</ul>

<h2>附：文件索引</h2>
<table><thead><tr><th>文件</th><th>内容</th></tr></thead><tbody>
<tr><td>model/FINAL_summary.html</td><td>总结报告</td></tr>
<tr><td>model/decomposition_report.html</td><td>缺口的 84/16 拆解</td></tr>
<tr><td>bench/mpv_barrier_report.html</td><td>MPV 势垒计算</td></tr>
<tr><td>bench/live_2d.html</td><td>加成反应二维势能面</td></tr>
<tr><td>model/scorecard.html</td><td>全部实验点记分卡</td></tr>
<tr><td>model/geometry_report.html</td><td>A/V 的几何来源</td></tr>
<tr><td>protocol-lit/README.md</td><td>文献定量数据集</td></tr>
</tbody></table>

</div></body></html>"""

with open(os.path.join(HERE, "速率怎么算的.html"), "w") as f:
    f.write(html)
print("written", os.path.join(HERE, "速率怎么算的.html"))
