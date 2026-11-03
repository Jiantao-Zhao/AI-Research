#!/usr/bin/env python3
"""构象系综采样：看势垒能不能收敛到文献值。

对 Agem-like 络合物（CH3MgCl(THF)(乙醛)）：
  1. 随机转动配位 THF 与乙醛（绕各自的 Mg–O 轴，加分面倾斜）和甲基，生成 N 个起点
  2. 各自弛豫，按能量去重
  3. 每个保留的构象做 C···C 约束扫描，取势垒 ΔE‡
  4. 按 298 K Boltzmann 权重平均速率常数（不是平均势垒），再换算成等效 ΔG‡

终点对照：文献 ΔA‡ = 8.7 ± 0.8 kcal/mol（注意那是自由能，这里是电子能）。
"""
import argparse, json, os, time
import numpy as np
import torch
torch.serialization.add_safe_globals([slice])
from ase.io import read, write
from ase.constraints import FixBondLength
from ase.optimize import BFGS
from mace.calculators import mace_polar

HERE = os.path.dirname(os.path.abspath(__file__))
KCAL = 23.0605
R_KCAL = 1.987e-3
T = 298.15
KBH = 6.212e12  # k_B T / h at 298 K

MG, CMET, CCAR = 0, 1, 20
THF = list(range(6, 19))     # 13 原子
ACA = list(range(19, 26))    # O, C, H, C, H, H, H


def rot(axis, angle):
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    c, s = np.cos(angle), np.sin(angle)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) * c + s * K + (1 - c) * np.outer(axis, axis)


def spin(atoms, idx, anchor, axis, angle):
    """把一组原子绕通过 anchor 的 axis 转 angle。"""
    Rm = rot(axis, angle)
    P = atoms.positions
    A = np.asarray(anchor, float)
    for i in idx:
        P[i] = A + Rm @ (P[i] - A)


def make_conformer(base, rng, k):
    at = base.copy()
    mg = at.positions[MG].copy()
    # 配位 THF 的氧
    o_thf = [i for i in THF if at[i].symbol == "O"][0]
    # THF：绕 Mg–O 轴转，再加一个垂直方向的倾斜
    v = at.positions[o_thf] - mg
    spin(at, THF, mg, v, rng.uniform(0, 2 * np.pi))
    perp = np.cross(v, [0, 0, 1.0])
    if np.linalg.norm(perp) < 1e-3:
        perp = np.cross(v, [1.0, 0, 0])
    spin(at, THF, mg, perp, rng.uniform(-0.45, 0.45))
    # 乙醛：绕 Mg–O 轴转 + 倾斜
    v2 = at.positions[19] - mg
    spin(at, ACA, mg, v2, rng.uniform(0, 2 * np.pi))
    perp2 = np.cross(v2, [0, 0, 1.0])
    if np.linalg.norm(perp2) < 1e-3:
        perp2 = np.cross(v2, [1.0, 0, 0])
    spin(at, ACA, mg, perp2, rng.uniform(-0.35, 0.35))
    # 甲基自转
    spin(at, [CMET, 2, 3, 4], mg, at.positions[CMET] - mg, rng.uniform(0, 2 * np.pi / 3))
    at.info["charge"] = 0
    at.info["spin"] = 1
    return at


def relax(at, calc, fmax=0.10, steps=400):
    at.calc = calc
    opt = BFGS(at, logfile=None, maxstep=0.10)
    try:
        opt.run(fmax=fmax, steps=steps)
    except Exception:
        pass
    return float(at.get_potential_energy())


def barrier(at, calc, r_values, step=0.15, fmax=0.12, maxsteps=120):
    """固定 C1–C20，小步逼近每个目标距离，每步后弛豫。返回势垒与剖面。"""
    e0 = float(at.get_potential_energy())
    work = at.copy()
    work.calc = calc
    prof = []
    for r in r_values:
        # 小步走到目标
        for _ in range(80):
            cur = float(work.get_distance(CMET, CCAR))
            gap = float(r) - cur
            if abs(gap) < 3e-3:
                break
            d = float(np.clip(gap, -step, step))
            u = (work.positions[CCAR] - work.positions[CMET]) / cur
            work.positions[CMET] -= u * d / 2
            work.positions[CCAR] += u * d / 2
            work.set_constraint(FixBondLength(CMET, CCAR))
            o = BFGS(work, logfile=None, maxstep=0.05)
            try:
                o.run(fmax=0.30, steps=maxsteps)
            except Exception:
                pass
        work.set_constraint([FixBondLength(CMET, CCAR)])
        o = BFGS(work, logfile=None, maxstep=0.05)
        try:
            o.run(fmax=fmax, steps=150)
        except Exception:
            pass
        e = float(work.get_potential_energy())
        work.set_constraint()
        D = work.get_all_distances(); np.fill_diagonal(D, 9.0)
        prof.append({"r": round(float(r), 2), "dE": float((e - e0) * KCAL),
                     "min_d": round(float(D.min()), 2)})
    best = max(prof, key=lambda p: p["dE"])
    return best["dE"], best["r"], prof


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=14, help="生成的起点数")
    ap.add_argument("--keep", type=int, default=10, help="保留的构象数")
    ap.add_argument("--seed", type=int, default=20261103)
    ap.add_argument("--outdir", default=os.path.join(HERE, "ensemble"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    calc = mace_polar(model="polar-1-m", device="cpu", default_dtype="float64",
                      pbc_handling="realspace")
    base = read(os.path.join(HERE, "Agem_like_min.xyz"))
    base.info["charge"] = 0
    base.info["spin"] = 1
    base.calc = calc
    e_ref = float(base.get_potential_energy())

    rng = np.random.default_rng(a.seed)
    R_SCAN = [3.0, 2.6, 2.35, 2.15, 1.95]

    cands = []
    t0 = time.time()
    for k in range(a.n):
        at = make_conformer(base, rng, k)
        e = relax(at, calc)
        de = (e - e_ref) * KCAL
        print(f"[gen {k+1:2d}/{a.n}] E_rel = {de:+7.2f} kcal/mol", flush=True)
        cands.append({"k": k, "atoms": at, "E": e, "dE": de})

    # 按能量排序，去掉能量过高的（> 6 kcal/mol）与近乎重复的
    cands.sort(key=lambda c: c["E"])
    kept = []
    for c in cands:
        if c["dE"] > 6.0:
            continue
        dup = False
        for q in kept:
            rmsd = np.sqrt(np.mean(np.sum((c["atoms"].positions - q["atoms"].positions) ** 2, axis=1)))
            if rmsd < 0.15:
                dup = True
                break
        if not dup:
            kept.append(c)
        if len(kept) >= a.keep:
            break

    print(f"\n保留 {len(kept)} 个构象，开始扫描势垒", flush=True)
    results = []
    for j, c in enumerate(kept):
        at = c["atoms"]
        at.calc = calc
        b, rmax, prof = barrier(at, calc, R_SCAN)
        write(os.path.join(a.outdir, f"conf{j:02d}.xyz"), at)
        results.append({"idx": j, "dE_conf": c["dE"], "barrier": b, "r": rmax, "profile": prof})
        bad = any(p["min_d"] < 0.85 for p in prof)
        print(f"[scan {j+1:2d}/{len(kept)}] E_conf={c['dE']:+6.2f}  "
              f"ΔE‡={b:6.2f} kcal/mol at r={rmax:.2f}" + ("  [clash]" if bad else ""), flush=True)

    # Boltzmann 权重（按构象能量），再按速率加权
    dEs = np.array([r["dE_conf"] for r in results])
    bars = np.array([r["barrier"] for r in results])
    w = np.exp(-(dEs - dEs.min()) / (R_KCAL * T))
    w /= w.sum()
    k_i = KBH * np.exp(-bars / (R_KCAL * T))
    k_eff = float((w * k_i).sum())
    bar_eff = float(-R_KCAL * T * np.log(k_eff / KBH))

    summary = {
        "n_conformers": len(results),
        "barriers": [round(float(x), 2) for x in bars],
        "barrier_mean": round(float(bars.mean()), 2),
        "barrier_std": round(float(bars.std()), 2),
        "barrier_min": round(float(bars.min()), 2),
        "barrier_max": round(float(bars.max()), 2),
        "weights": [round(float(x), 3) for x in w],
        "barrier_boltzmann": round(bar_eff, 2),
        "literature": 8.7,
        "wall_s": round(time.time() - t0, 1),
        "results": results,
    }
    with open(os.path.join(a.outdir, "ensemble.json"), "w") as f:
        json.dump(summary, f, indent=1)

    print("\n=== 汇总 ===")
    print(f"构象数        {len(results)}")
    print(f"势垒范围      {bars.min():.2f} – {bars.max():.2f} kcal/mol")
    print(f"算术平均      {bars.mean():.2f} ± {bars.std():.2f}")
    print(f"Boltzmann 平均 {bar_eff:.2f}  (按速率加权)")
    print(f"文献 ΔA‡      8.7 ± 0.8")
    print(f"耗时          {summary['wall_s']} s")


if __name__ == "__main__":
    main()
