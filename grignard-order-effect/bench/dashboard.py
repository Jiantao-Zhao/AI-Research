#!/usr/bin/env python3
"""Build a self-contained interactive HTML page for one structure or a scan series.

Usage:
    python dashboard.py --json Agem_like_fine.json --glob 'Agem_like_fine_*.xyz' --out Agem.html
    python dashboard.py --xyz some.xyz --out one.html
"""
import argparse, glob, json, os, sys
import numpy as np


def read(path):
    """Minimal extxyz reader: returns (symbols, Nx3 positions)."""
    L = open(path).read().splitlines()
    n = int(L[0].split()[0])
    sym, pos = [], []
    for ln in L[2:2 + n]:
        f = ln.split()
        sym.append(f[0])
        pos.append([float(f[1]), float(f[2]), float(f[3])])
    return sym, np.array(pos)

HERE = os.path.dirname(os.path.abspath(__file__))


def xyz_block(atoms):
    sym, pos = atoms
    return {"s": sym, "p": [[round(float(v), 4) for v in row] for row in pos]}


def build_payload_scan(prof_json, xyz_files, title, lit=None, i=None, j=None):
    d = json.load(open(prof_json))
    prof = d["prof"] if isinstance(d, dict) and "prof" in d else d
    files = sorted(xyz_files)
    pts = []
    for k, rec in enumerate(prof):
        if k >= len(files):
            break
        at = read(files[k])
        pts.append({"r": rec.get("r"), "dE": rec.get("dE_kcal"),
                    "rMgC": rec.get("r_MgC"), "xyz": xyz_block(at),
                    "file": os.path.basename(files[k])})
    n = len(pts[0]["xyz"]["s"]) if pts else 0
    if i is None or j is None:
        # default: last two carbons that are closest at the highest-energy frame
        k = int(np.argmax([p["dE"] or 0 for p in pts])) if pts else 0
        P = np.array(pts[k]["xyz"]["p"])
        C = [a for a, s in enumerate(pts[k]["xyz"]["s"]) if s == "C"]
        best = None
        for a in range(len(C)):
            for b in range(a + 1, len(C)):
                dd = np.linalg.norm(P[C[a]] - P[C[b]])
                if best is None or dd < best[2]:
                    best = (C[a], C[b], dd)
        i, j = best[0], best[1]
    return {"title": title, "lit": lit, "points": pts, "i": int(i), "j": int(j),
            "n_atoms": n}


def build_payload_single(xyz, title):
    at = read(xyz)
    return {"title": title, "lit": None,
            "points": [{"r": None, "dE": None, "rMgC": None,
                        "xyz": xyz_block(at), "file": os.path.basename(xyz)}],
            "i": 0, "j": 1, "n_atoms": len(at[0])}


TEMPLATE = r"""<!doctype html>
<html lang="zh"><head><meta charset="utf-8">
<title>__TITLE__</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
 body{margin:0;font:13px/1.5 -apple-system,"PingFang SC","Hiragino Sans GB",sans-serif;color:#222;background:#fafafa}
 header{padding:10px 16px;background:#fff;border-bottom:1px solid #e6e6e6}
 h1{font-size:15px;margin:0;font-weight:600}
 .sub{color:#777;font-size:12px;margin-top:2px}
 .wrap{display:flex;gap:0;height:calc(100vh - 62px)}
 .left{width:430px;padding:14px 16px;overflow:auto;background:#fff;border-right:1px solid #e6e6e6}
 .right{flex:1;position:relative}
 #viewer{width:100%;height:100%}
 .ctrl{display:flex;align-items:center;gap:8px;margin:10px 0}
 button{padding:5px 12px;border:1px solid #ccc;background:#fff;border-radius:6px;cursor:pointer;font-size:12px}
 button:hover{background:#f0f0f0}
 input[type=range]{flex:1}
 .readout{background:#f6f7f9;border:1px solid #e6e6e6;border-radius:8px;padding:10px;margin:10px 0;font-variant-numeric:tabular-nums}
 .readout b{font-size:15px}
 table{width:100%;border-collapse:collapse;font-size:12px;font-variant-numeric:tabular-nums}
 th,td{padding:3px 6px;text-align:right;border-bottom:1px solid #f0f0f0}
 th:first-child,td:first-child{text-align:left}
 tr.cur{background:#fff3cd}
 tbody tr{cursor:pointer}
 tbody tr:hover{background:#f5f5f5}
 .tag{display:inline-block;padding:1px 6px;border-radius:4px;font-size:11px;background:#eef;color:#335}
 .lit{color:#c0392b}
</style></head><body>
<header>
 <h1>__TITLE__</h1>
 <div class="sub">拖动滑块或用 ←/→ 逐点查看；点表格行或曲线上的点跳转；红色高亮的是形成/断裂的键</div>
</header>
<div class="wrap">
 <div class="left">
   <canvas id="chart" width="400" height="220" style="width:100%;height:220px"></canvas>
   <div class="ctrl">
     <button id="prev">←</button>
     <input id="slider" type="range" min="0" max="0" value="0" step="1">
     <button id="next">→</button>
   </div>
   <div class="readout">
     <div>第 <b id="idx">1</b> / <span id="tot">1</span> 点 &nbsp;<span class="tag" id="fname"></span></div>
     <div style="margin-top:6px">C···C = <b id="rc">–</b> Å &nbsp;|&nbsp; Mg–C = <b id="rmg">–</b> Å</div>
     <div style="margin-top:4px">ΔE = <b id="de">–</b> kcal/mol <span id="litbox"></span></div>
   </div>
   <table>
     <thead><tr><th>#</th><th>r / Å</th><th>ΔE</th><th>Mg–C</th></tr></thead>
     <tbody id="tbody"></tbody>
   </table>
 </div>
 <div class="right"><div id="viewer"></div></div>
</div>
<script>
const D = __DATA__;
let cur = 0;
const viewer = $3Dmol.createViewer("viewer", {backgroundColor:"white"});
const slider = document.getElementById("slider");

function xyzString(k){
  const p = D.points[k];
  let lines = [p.xyz.s.length, "frame " + k];
  for (let a=0; a<p.xyz.s.length; a++){
    lines.push(p.xyz.s[a] + " " + p.xyz.p[a][0].toFixed(4) + " " + p.xyz.p[a][1].toFixed(4) + " " + p.xyz.p[a][2].toFixed(4));
  }
  return lines.join("\n");
}
function show(k){
  if (k < 0 || k >= D.points.length) return;
  cur = k;
  const p = D.points[k];
  viewer.removeAllModels();
  viewer.addModel(xyzString(k), "xyz");
  viewer.setStyle({}, {stick:{radius:0.14, colorscheme:"Jmol"}, sphere:{scale:0.26}});
  // 让 Mg / Cl 一眼可辨：Jmol 里两者都是绿色，容易混
  viewer.setStyle({elem:"Mg"}, {sphere:{scale:0.34, color:"#8e44ad"}, stick:{radius:0.15, color:"#8e44ad"}});
  viewer.setStyle({elem:"Cl"}, {sphere:{scale:0.34, color:"#27ae60"}, stick:{radius:0.15, color:"#27ae60"}});
  viewer.setStyle({elem:"O"},  {sphere:{scale:0.30, color:"#e74c3c"}, stick:{radius:0.14, color:"#e74c3c"}});
  viewer.addStyle({index:[D.i, D.j]}, {sphere:{scale:0.40, color:"#f39c12"}});
  viewer.addStyle({index:[D.i, D.j]}, {stick:{radius:0.20, color:"#f39c12"}});
  viewer.zoomTo(); viewer.render();
  document.getElementById("idx").textContent = (k+1);
  document.getElementById("fname").textContent = p.file || "";
  document.getElementById("rc").textContent = (p.r==null? "–" : p.r.toFixed(2));
  document.getElementById("rmg").textContent = (p.rMgC==null? "–" : p.rMgC.toFixed(2));
  document.getElementById("de").textContent = (p.dE==null? "–" : p.dE.toFixed(2));
  slider.value = k;
  [...document.querySelectorAll("#tbody tr")].forEach((tr,n)=>tr.className = (n===k?"cur":""));
  drawChart();
}
function drawChart(){
  const c = document.getElementById("chart"), ctx = c.getContext("2d");
  const W=c.width, H=c.height, pad=38;
  ctx.clearRect(0,0,W,H);
  const pts = D.points.filter(p=>p.r!=null && p.dE!=null);
  if (pts.length < 2){ ctx.fillStyle="#888"; ctx.fillText("单结构：无能量曲线", pad, H/2); return; }
  const rs = pts.map(p=>p.r), es = pts.map(p=>p.dE);
  const xmin=Math.min(...rs), xmax=Math.max(...rs);
  let ymin=Math.min(...es), ymax=Math.max(...es);
  const span=Math.max(ymax-ymin, 1e-6); ymin-=span*0.12; ymax+=span*0.12;
  const X = r => pad + (xmax - r)/(xmax-xmin||1) * (W-2*pad);
  const Y = e => H-pad - (e-ymin)/(ymax-ymin||1) * (H-2*pad);
  if (D.lit != null){
    ctx.strokeStyle="#e0a0a0"; ctx.setLineDash([4,3]); ctx.beginPath();
    ctx.moveTo(pad, Y(D.lit)); ctx.lineTo(W-pad, Y(D.lit)); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle="#c0392b"; ctx.font="11px sans-serif";
    ctx.fillText("文献 ΔA‡ " + D.lit.toFixed(1), pad+2, Y(D.lit)-4);
  }
  ctx.strokeStyle="#333"; ctx.lineWidth=1.4; ctx.beginPath();
  pts.forEach((p,n)=>{ const x=X(p.r), y=Y(p.dE); n? ctx.lineTo(x,y): ctx.moveTo(x,y); });
  ctx.stroke();
  pts.forEach((p,n)=>{
    const x=X(p.r), y=Y(p.dE);
    ctx.beginPath(); ctx.arc(x,y,4,0,7);
    ctx.fillStyle = (n===cur)? "#c0392b" : "#2c6fb5"; ctx.fill();
  });
  ctx.fillStyle="#666"; ctx.font="11px sans-serif";
  ctx.fillText("ΔE / kcal mol⁻¹", 4, 14);
  ctx.fillText("C···C / Å", W-58, H-10);
  const xc = X(D.points[cur] && D.points[cur].r != null ? D.points[cur].r : rs[0]);
  ctx.strokeStyle="#c0392b55"; ctx.beginPath(); ctx.moveTo(xc,pad-6); ctx.lineTo(xc,H-pad); ctx.stroke();
}
const tb = document.getElementById("tbody");
D.points.forEach((p,k)=>{
  const tr = document.createElement("tr");
  tr.innerHTML = "<td>"+(k+1)+"</td><td>"+(p.r==null?"–":p.r.toFixed(2))+"</td><td>"+
                 (p.dE==null?"–":p.dE.toFixed(2))+"</td><td>"+
                 (p.rMgC==null?"–":p.rMgC.toFixed(2))+"</td>";
  tr.onclick = ()=>show(k);
  tb.appendChild(tr);
});
document.getElementById("litbox").innerHTML = (D.lit!=null)
  ? '<span class="lit">文献 ΔA‡ = '+D.lit.toFixed(1)+'</span>' : "";
slider.max = D.points.length-1;
slider.oninput = e=>show(+e.target.value);
document.getElementById("prev").onclick = ()=>show(cur-1);
document.getElementById("next").onclick = ()=>show(cur+1);
document.addEventListener("keydown", e=>{
  if(e.key==="ArrowLeft") show(cur-1);
  if(e.key==="ArrowRight") show(cur+1);
});
document.getElementById("tot").textContent = D.points.length;
let start = 0;
if (D.points.some(p=>p.dE!=null)){
  let best=-1e9; D.points.forEach((p,k)=>{ if(p.dE!=null && p.dE>best){best=p.dE; start=k;} });
}
show(start);
function fit(){
  const el = document.getElementById("viewer");
  const w = el.clientWidth || 900, h = el.clientHeight || 800;
  try { if (viewer.setWidth)  viewer.setWidth(w);  } catch(e){}
  try { if (viewer.setHeight) viewer.setHeight(h); } catch(e){}
  viewer.zoomTo();
  // 3Dmol frames using the canvas backing store; on retina (dpr=2) that magnifies
  // the model by the device pixel ratio, so compensate explicitly.
  const cv = el.querySelector("canvas");
  const dpr = (cv && el.clientWidth) ? (cv.width / el.clientWidth) : 1;
  viewer.zoom(1.15 / (dpr || 1));
  viewer.render();
}
setTimeout(fit, 120); setTimeout(fit, 600); setTimeout(fit, 1500);
window.addEventListener("resize", fit);
</script>
</body></html>
"""


def write_html(payload, out):
    html = (TEMPLATE
            .replace("__TITLE__", payload["title"])
            .replace("__DATA__", json.dumps(payload, ensure_ascii=False)))
    with open(out, "w") as f:
        f.write(html)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--json"); ap.add_argument("--glob"); ap.add_argument("--xyz")
    ap.add_argument("--title", default="structure")
    ap.add_argument("--lit", type=float, default=None)
    ap.add_argument("--i", type=int); ap.add_argument("--j", type=int)
    ap.add_argument("--out", default="dashboard.html")
    a = ap.parse_args()
    if a.xyz:
        pl = build_payload_single(a.xyz, a.title)
    else:
        files = sorted(glob.glob(a.glob)) if a.glob else []
        pl = build_payload_scan(a.json, files, a.title, a.lit, a.i, a.j)
    print(write_html(pl, a.out), "|", len(pl["points"]), "points")
