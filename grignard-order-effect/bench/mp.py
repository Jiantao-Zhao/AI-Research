#!/usr/bin/env python3
"""mp — MACE-POLAR 操作台。

每个子命令都是「跑一步、立刻打印人能读的结果、落盘一个文件」，方便手动控制节奏。

  mp sp    <xyz>                              单点：能量 / 力 / 电荷 / 偶极
  mp opt   <xyz> [--fmax 0.05]                几何优化，逐步打印
  mp scan  <xyz> --i 1 --j 20 --r 4.4:1.8:0.2 约束扫描，逐点存几何并自动出看板
  mp view  <xyz> | --json X.json --glob 'X_*.xyz'   生成可交互 HTML

所有命令都接受 --charge / --spin（spin 是多重度 2S+1，单重态填 1）。
"""
import argparse, glob, json, os, sys, time
import numpy as np

from ase.io import read, write
from ase.constraints import FixBondLength
from ase.optimize import BFGS

KCAL = 23.0605
HERE = os.path.dirname(os.path.abspath(__file__))
_CALC = None


def get_calc(model="polar-1-m", dtype="float64"):
    global _CALC
    if _CALC is None:
        import torch
        torch.serialization.add_safe_globals([slice])
        from mace.calculators import mace_polar
        _CALC = mace_polar(model=model, device="cpu", default_dtype=dtype,
                           pbc_handling="realspace")
    return _CALC


def load(path, charge=0, spin=1, calc=None):
    at = read(path)
    at.info["charge"] = charge
    at.info["spin"] = spin
    at.calc = calc or get_calc()
    return at


def hr(title):
    print("\n" + "─" * 68)
    print(title)
    print("─" * 68, flush=True)


# ------------------------------------------------------------------ sp
def cmd_sp(a):
    at = load(a.xyz, a.charge, a.spin)
    t0 = time.time()
    e = at.get_potential_energy()
    f = at.get_forces()
    dt = time.time() - t0
    hr(f"单点  {os.path.basename(a.xyz)}")
    print(f"原子数          {len(at)}")
    print(f"能量            {e:.6f} eV   ({e/len(at):.4f} eV/atom)")
    print(f"最大力          {np.abs(f).max():.4f} eV/Å")
    print(f"耗时            {dt:.2f} s")
    try:
        p = at.calc.results["density_coefficients"]
        ch = p[:, 0]
        print("\n部分电荷（按元素平均）")
        for s in sorted(set(at.get_chemical_symbols())):
            m = [i for i, x in enumerate(at.get_chemical_symbols()) if x == s]
            print(f"  {s:3s}  n={len(m):2d}   mean={ch[m].mean():+.3f} e")
        print(f"  合计 {ch.sum():+.4f} e")
    except Exception:
        pass
    if a.save:
        write(a.save, at)
        print(f"\n写出 -> {a.save}")


# ------------------------------------------------------------------ opt
def cmd_opt(a):
    at = load(a.xyz, a.charge, a.spin)
    e0 = at.get_potential_energy()
    traj = []
    hr(f"几何优化  {os.path.basename(a.xyz)}   fmax={a.fmax}")

    def log():
        e = at.get_potential_energy()
        fmax = np.abs(at.get_forces()).max()
        traj.append((len(traj), e, fmax))
        if len(traj) % a.every == 0:
            print(f"  step {len(traj):4d}   E={e:14.6f} eV   fmax={fmax:.3f}", flush=True)

    opt = BFGS(at, logfile=None, maxstep=a.maxstep)
    opt.attach(log)
    t0 = time.time()
    opt.run(fmax=a.fmax, steps=a.steps)
    e1 = at.get_potential_energy()
    hr("完成")
    print(f"步数            {opt.get_number_of_steps()}")
    print(f"能量            {e0:.6f} -> {e1:.6f} eV   (Δ = {(e1-e0)*KCAL:+.2f} kcal/mol)")
    print(f"最大力          {np.abs(at.get_forces()).max():.4f} eV/Å")
    print(f"耗时            {time.time()-t0:.1f} s")
    out = a.save or os.path.splitext(a.xyz)[0] + "_opt.xyz"
    write(out, at)
    print(f"写出            {out}")


# ------------------------------------------------------------------ scan
def parse_range(s):
    f, t, st = [float(x) for x in s.split(":")]
    n = int(round(abs(f - t) / st)) + 1
    return list(np.round(np.linspace(f, t, n), 3))


def cmd_scan(a):
    at = load(a.xyz, a.charge, a.spin)
    i, j = a.i, a.j
    os.makedirs(a.outdir, exist_ok=True)
    geo_dir = os.path.join(a.outdir, "geoms")
    os.makedirs(geo_dir, exist_ok=True)
    dists = parse_range(a.r)
    e0 = at.get_potential_energy()
    tag = a.tag or os.path.splitext(os.path.basename(a.xyz))[0]

    hr(f"约束扫描  {tag}   键 {i}-{j}   {dists[0]} -> {dists[-1]} Å  ({len(dists)} 点)")
    prof, work = [], at.copy()
    work.calc = at.calc
    t0 = time.time()
    for r in dists:
        v = work.positions[j] - work.positions[i]
        n = np.linalg.norm(v)
        u = v / n
        work.positions[i] -= u * (r - n) / 2
        work.positions[j] += u * (r - n) / 2
        work.set_constraint(FixBondLength(i, j))
        opt = BFGS(work, logfile=None, maxstep=a.maxstep)
        try:
            opt.run(fmax=a.fmax, steps=a.steps)
        except Exception as exc:
            print(f"  !! r={r}: {exc}", flush=True)
        e = work.get_potential_energy()
        rmg = float(work.get_distance(0, i))
        work.set_constraint()
        fn = os.path.join(geo_dir, f"{tag}_r{r:.2f}.xyz")
        write(fn, work)
        prof.append({"r": float(r), "E_eV": float(e), "dE_kcal": float((e - e0) * KCAL),
                     "r_MgC": round(rmg, 3), "file": os.path.basename(fn)})
        print(f"  r={r:5.2f}  ΔE={prof[-1]['dE_kcal']:7.2f} kcal/mol   "
              f"Mg–C={rmg:.2f} Å   ({opt.get_number_of_steps()} 步)", flush=True)

    emax = max(prof, key=lambda p: p["dE_kcal"])
    hr("结果")
    print(f"最高点          r = {emax['r']:.2f} Å,  ΔE = {emax['dE_kcal']:.2f} kcal/mol")
    print(f"总耗时          {time.time()-t0:.1f} s")
    out = os.path.join(a.outdir, f"{tag}_scan.json")
    json.dump({"tag": tag, "i": i, "j": j, "e0": e0, "prof": prof},
              open(out, "w"), indent=1)
    print(f"剖面            {out}")
    print(f"几何            {geo_dir}/  ({len(prof)} 个 extxyz)")

    if not a.no_view:
        sys.path.insert(0, HERE)
        import dashboard
        pl = dashboard.build_payload_scan(out, sorted(glob.glob(os.path.join(geo_dir, f"{tag}_r*.xyz"))),
                                          f"{tag}：{i}–{j} 键扫描", a.lit, i, j)
        html = os.path.join(a.outdir, f"{tag}_dashboard.html")
        dashboard.write_html(pl, html)
        print(f"看板            {html}")


# ------------------------------------------------------------------ view
def cmd_view(a):
    sys.path.insert(0, HERE)
    import dashboard
    if a.xyz:
        pl = dashboard.build_payload_single(a.xyz, a.title or os.path.basename(a.xyz))
    else:
        files = sorted(glob.glob(a.glob))
        pl = dashboard.build_payload_scan(a.json, files, a.title or "scan", a.lit, a.i, a.j)
    out = a.out or "dashboard.html"
    dashboard.write_html(pl, out)
    print(f"写出 {out}  ({len(pl['points'])} 点)")


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(prog="mp", description="MACE-POLAR 操作台")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("xyz")
        p.add_argument("--charge", type=int, default=0)
        p.add_argument("--spin", type=int, default=1, help="多重度 2S+1（单重态=1）")
        p.add_argument("--model", default="polar-1-m")

    p = sub.add_parser("sp", help="单点")
    common(p); p.add_argument("--save"); p.set_defaults(f=cmd_sp)

    p = sub.add_parser("opt", help="几何优化")
    common(p)
    p.add_argument("--fmax", type=float, default=0.05)
    p.add_argument("--steps", type=int, default=400)
    p.add_argument("--maxstep", type=float, default=0.10)
    p.add_argument("--every", type=int, default=5)
    p.add_argument("--save")
    p.set_defaults(f=cmd_opt)

    p = sub.add_parser("scan", help="约束扫描")
    common(p)
    p.add_argument("--i", type=int, required=True)
    p.add_argument("--j", type=int, required=True)
    p.add_argument("--r", default="4.4:1.8:0.2", help="from:to:step")
    p.add_argument("--fmax", type=float, default=0.08)
    p.add_argument("--steps", type=int, default=250)
    p.add_argument("--maxstep", type=float, default=0.12)
    p.add_argument("--outdir", default="scan_out")
    p.add_argument("--tag"); p.add_argument("--lit", type=float)
    p.add_argument("--no-view", action="store_true")
    p.set_defaults(f=cmd_scan)

    p = sub.add_parser("view", help="生成交互看板")
    p.add_argument("xyz", nargs="?")
    p.add_argument("--json"); p.add_argument("--glob")
    p.add_argument("--title"); p.add_argument("--lit", type=float)
    p.add_argument("--i", type=int); p.add_argument("--j", type=int)
    p.add_argument("--out")
    p.set_defaults(f=cmd_view)

    a = ap.parse_args()
    if a.cmd in ("sp", "opt", "scan"):
        get_calc(a.model)
    a.f(a)


if __name__ == "__main__":
    main()
