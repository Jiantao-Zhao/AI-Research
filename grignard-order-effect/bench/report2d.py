#!/usr/bin/env python3
"""Render the 2D scan (live or finished) as a self-contained HTML report.

    python report2d.py            # one shot
    python report2d.py --watch    # regenerate whenever the partial json changes
"""
import argparse, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
JSON = os.path.join(HERE, "scan2d_partial.json")
GEO = os.path.join(HERE, "geoms_2d")
OUT = os.path.join(HERE, "live_2d.html")


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0].split()[0])
    sym, pos = [], []
    for ln in L[2:2 + n]:
        f = ln.split()
        sym.append(f[0])
        pos.append([round(float(f[1]), 4), round(float(f[2]), 4), round(float(f[3]), 4)])
    return {"s": sym, "p": pos}


def build(state):
    pts = state["points"]
    for p in pts:
        fp = os.path.join(GEO, p["xyz"])
        if os.path.exists(fp):
            p["geom"] = read_xyz(fp)
        else:
            p["geom"] = None
    payload = {
        "tag": state["tag"], "cc": state["cc"], "mgc": state["mgc"],
        "done": state["done"], "total": state["total"],
        "finished": state["finished"], "running": state.get("running"),
        "points": pts, "i": 1, "j": 20, "mg": 0,
    }
    return payload


TEMPLATE = r"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
__REFRESH__
<title>__TITLE__</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
 body{margin:0;font:13px/1.6 -apple-system,"PingFang SC","Hiragino Sans GB",sans-serif;color:#222;background:#f7f8fa}
 header{padding:14px 20px;background:#fff;border-bottom:1px solid #e6e6e6}
 h1{font-size:17px;margin:0 0 4px}
 .meta{color:#666;font-size:12px}
 .wrap{display:flex;height:calc(100vh - 78px)}
 .left{width:520px;overflow:auto;padding:16px 20px;background:#fff;border-right:1px solid #e6e6e6}
 .right{flex:1;position:relative}
 #viewer{position:absolute;inset:0}
 h2{font-size:13px;margin:18px 0 6px;color:#111;border-left:3px solid #2c6fb5;padding-left:8px}
 .box{background:#fafbfc;border:1px solid #e8e8e8;border-radius:8px;padding:10px 12px;font-size:12px;color:#333}
 .box b{color:#111}
 ul.goal{list-style:none;padding-left:0;margin:6px 0}
 ul.goal li{padding:3px 0 3px 14px;position:relative;font-size:12px}
 ul.goal li:before{content:"▸";position:absolute;left:0;color:#2c6fb5}
 .lvl{display:inline-block;font-size:10px;padding:0 5px;border-radius:3px;background:#e8eef7;color:#2c6fb5;margin-right:6px}
 .prog{height:6px;background:#eee;border-radius:3px;overflow:hidden;margin:8px 0}
 .prog>div{height:100%;background:#2c6fb5;width:0%}
 table{width:100%;border-collapse:collapse;font-size:11.5px;font-variant-numeric:tabular-nums}
 th,td{padding:2px 5px;text-align:right;border-bottom:1px solid #f0f0f0}
 th:nth-child(-n+2),td:nth-child(-n+2){text-align:left}
 tr.cur{background:#fff3cd}
 tbody tr{cursor:pointer} tbody tr:hover{background:#f2f2f2}
 .readout{font-variant-numeric:tabular-nums;background:#f2f6fb;border:1px solid #dde6f2;border-radius:8px;padding:8px 10px;margin:8px 0}
 .readout b{font-size:15px}
 .note{font-size:11.5px;color:#777}
</style></head><body>
<header>
 <h1>__TITLE__</h1>
 <div class="meta">__SUBTITLE__</div>
 <div class="prog"><div id="bar"></div></div>
</header>
<div class="wrap">
 <div class="left">
   <h2>这次在算什么</h2>
   <div class="box" id="method"></div>

   <h2>目标链条</h2>
   <ul class="goal" id="goals"></ul>

   <h2>二维能量图（C···C × Mg–C）</h2>
   <canvas id="map" width="470" height="300" style="width:100%;height:auto"></canvas>
   <div class="note">横轴 C···C（右＝成键），纵轴 Mg–C（上＝断镁）。灰格＝还没算。</div>

   <h2>当前结构</h2>
   <div class="readout" id="readout">–</div>

   <h2>已完成的点</h2>
   <table><thead><tr><th>#</th><th>C···C</th><th>Mg–C</th><th>ΔE</th></tr></thead>
   <tbody id="tbody"></tbody></table>
 </div>
 <div class="right"><div id="viewer"></div></div>
</div>
<script>
const D = __DATA__;
let cur = -1;

document.getElementById("method").innerHTML = `
 <b>势函数</b> MACE-POLAR-1-M（OMol25 预训练，ωB97M-V 级）零样本<br>
 <b>体系</b> CH₃MgCl(THF)(乙醛)，26 原子，气相、单构象<br>
 <b>约束</b> 同时固定 |C1–C20| 与 |Mg0–C1|，其余全部放开，BFGS fmax 0.10 eV/Å<br>
 <b>能量</b> 电子能，相对弛豫前体；<b>无</b> ZPE / 熵 / 显式溶剂<br>
 <b>耗时</b> ${D.total} 点，约 ${Math.max(1,Math.round(D.total*35/60))} 分钟`;

document.getElementById("goals").innerHTML = `
 <li><span class="lvl">大目标</span>同一批料换加料顺序，为什么产物分布不同</li>
 <li><span class="lvl">中目标</span>判定它来自「池子没来得及平衡」还是「底物占了金属表面」</li>
 <li><span class="lvl">小目标</span>加成速率常数 k 能不能用——决定 k_eff = Σxᵢkᵢ 是否可信</li>
 <li><span class="lvl">本次</span>一维扫描的最大值不是鞍点（0 虚频）。把它展成二维，看是否存在「边成键边断镁」的低能量通道</li>`;

const viewer = $3Dmol.createViewer("viewer", {backgroundColor:"white"});

function fit(){
  const el = document.getElementById("viewer");
  try{ viewer.setWidth(el.clientWidth); viewer.setHeight(el.clientHeight); }catch(e){}
  viewer.zoomTo();
  const cv = el.querySelector("canvas");
  const dpr = (cv && el.clientWidth) ? (cv.width/el.clientWidth) : 1;
  viewer.zoom(1.15/(dpr||1)); viewer.render();
}
function draw(k){
  const p = D.points[k]; if(!p || !p.geom) return;
  cur = k; location.hash = "p"+k;
  viewer.removeAllModels();
  let lines=[p.geom.s.length,"frame"];
  for(let a=0;a<p.geom.s.length;a++) lines.push(p.geom.s[a]+" "+p.geom.p[a].join(" "));
  viewer.addModel(lines.join("\n"),"xyz");
  viewer.setStyle({},{stick:{radius:0.14,colorscheme:"Jmol"},sphere:{scale:0.26}});
  viewer.setStyle({elem:"Mg"},{sphere:{scale:0.34,color:"#8e44ad"},stick:{radius:0.15,color:"#8e44ad"}});
  viewer.setStyle({elem:"Cl"},{sphere:{scale:0.34,color:"#27ae60"},stick:{radius:0.15,color:"#27ae60"}});
  viewer.setStyle({elem:"O"},{sphere:{scale:0.30,color:"#e74c3c"},stick:{radius:0.14,color:"#e74c3c"}});
  viewer.addStyle({index:[D.i,D.j]},{sphere:{scale:0.42,color:"#f39c12"}});
  viewer.zoomTo(); fit();
  document.getElementById("readout").innerHTML =
    `C···C = <b>${p.cc.toFixed(2)}</b> Å &nbsp;|&nbsp; Mg–C = <b>${p.mgc.toFixed(2)}</b> Å<br>`+
    `ΔE = <b>${p.dE_kcal>=0?"+":""}${p.dE_kcal.toFixed(2)}</b> kcal/mol`;
  [...document.querySelectorAll("#tbody tr")].forEach((tr,n)=>{
    tr.className = (+tr.dataset.k === k) ? "cur" : "";
  });
  drawMap();
}
function drawMap(){
  const c=document.getElementById("map"), ctx=c.getContext("2d");
  const W=c.width,H=c.height,pad=40;
  ctx.clearRect(0,0,W,H);
  const xs=D.cc, ys=D.mgc;
  const cw=(W-2*pad)/xs.length, ch=(H-2*pad)/ys.length;
  const idx={}; D.points.forEach(p=>idx[p.cc+","+p.mgc]=p);
  let lo=Infinity,hi=-Infinity;
  D.points.forEach(p=>{lo=Math.min(lo,p.dE_kcal);hi=Math.max(hi,p.dE_kcal);});
  const col=v=>{
    if(!isFinite(v)) return "#eee";
    const t=(v-lo)/((hi-lo)||1);
    const r=Math.round(255*t), b=Math.round(255*(1-t));
    return `rgb(${Math.round(80+175*t)},${Math.round(140-80*Math.abs(t-0.5)*2)},${Math.round(80+175*(1-t))})`;
  };
  ys.forEach((y,j)=>xs.forEach((x,i)=>{
    const px=pad+i*cw, py=H-pad-(j+1)*ch;
    const p=idx[x+","+y];
    ctx.fillStyle = p? col(p.dE_kcal) : "#eee";
    ctx.fillRect(px+1,py+1,cw-2,ch-2);
    ctx.strokeStyle="#fff"; ctx.strokeRect(px+1,py+1,cw-2,ch-2);
    if(p && p===D.points[cur]){
      ctx.strokeStyle="#111"; ctx.lineWidth=2; ctx.strokeRect(px+1,py+1,cw-2,ch-2); ctx.lineWidth=1;
    }
    if(p){
      ctx.fillStyle="#222"; ctx.font="9px sans-serif"; ctx.textAlign="center";
      ctx.fillText(p.dE_kcal.toFixed(1), px+cw/2, py+ch/2+3);
    }
  }));
  ctx.fillStyle="#555"; ctx.font="11px sans-serif"; ctx.textAlign="center";
  xs.forEach((x,i)=>ctx.fillText(x.toFixed(1), pad+(i+0.5)*cw, H-pad+14));
  ctx.textAlign="right";
  ys.forEach((y,j)=>ctx.fillText(y.toFixed(1), pad-5, H-pad-(j+0.5)*ch+4));
  ctx.save(); ctx.translate(12,H/2); ctx.rotate(-Math.PI/2);
  ctx.textAlign="center"; ctx.fillText("Mg–C / Å",0,0); ctx.restore();
  ctx.fillText("C···C / Å", W/2, H-6);
}
const tb=document.getElementById("tbody");
D.points.forEach((p,k)=>{
  const tr=document.createElement("tr"); tr.dataset.k=k;
  tr.innerHTML=`<td>${k+1}</td><td>${p.cc.toFixed(1)}</td><td>${p.mgc.toFixed(1)}</td>
                <td>${p.dE_kcal>=0?"+":""}${p.dE_kcal.toFixed(2)}</td>`;
  tr.onclick=()=>draw(k); tb.appendChild(tr);
});
document.getElementById("bar").style.width = (100*D.done/D.total).toFixed(1)+"%";
let start = D.points.length-1;
const m=location.hash.match(/^#p(\d+)$/); if(m) start=Math.min(+m[1], D.points.length-1);
if(D.points.length) draw(start); else drawMap();
</script></body></html>"""


def write(state):
    payload = build(state)
    title = f"二维势能面：{payload['tag']}（{payload['done']}/{payload['total']} 点）"
    sub = (f"C···C × Mg–C 双约束扫描 · 每点放开其余坐标 · "
           f"{'已完成' if payload['finished'] else ('正在计算 ' + str(payload['running']))}")
    refresh = "" if payload["finished"] else '<meta http-equiv="refresh" content="5">'
    html = (TEMPLATE.replace("__TITLE__", title)
                    .replace("__SUBTITLE__", sub)
                    .replace("__REFRESH__", refresh)
                    .replace("__DATA__", json.dumps(payload, ensure_ascii=False)))
    with open(OUT, "w") as f:
        f.write(html)
    return OUT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--watch", action="store_true")
    a = ap.parse_args()
    last = None
    while True:
        try:
            st = json.load(open(JSON))
        except Exception:
            time.sleep(1); continue
        key = (st["done"], st.get("finished"), (st.get("running") or {}).get("cc"))
        if key != last:
            write(st)
            print(f"{time.strftime('%H:%M:%S')} rendered {st['done']}/{st['total']}", flush=True)
            last = key
        if not a.watch:
            break
        if st.get("finished"):
            write(st)
            print("finished", flush=True)
            break
        time.sleep(1)


if __name__ == "__main__":
    main()
