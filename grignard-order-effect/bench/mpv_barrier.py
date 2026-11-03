#!/usr/bin/env python3
"""MPV 势垒：甲氧基镁把氢负离子交给乙醛。

反应（简化模型，气相）：
    CH3-O-MgCl + CH3CHO  ->  HCHO + CH3CH2-O-MgCl

过渡态是一个六元环：Mg-O-C(methoxy)-H···C(carbonyl)-O(carbonyl)。
用 MACE-POLAR 扫 H 从甲氧基碳迁移到羰基碳，得到 ΔE‡，再用 Eyring 换成 k。
"""
import json, os
import numpy as np
import torch
torch.serialization.add_safe_globals([slice])
from ase import Atoms
from ase.io import write
from ase.constraints import FixBondLength
from ase.optimize import BFGS
from mace.calculators import mace_polar

HERE = os.path.dirname(os.path.abspath(__file__))
KCAL = 23.0605


def build():
    """粗略起始构型：Mg 同时连着 Cl、甲氧基 O、乙醛 O。"""
    P, S = [], []

    def add(s, x, y, z):
        S.append(s); P.append([x, y, z])

    # Mg 在原点；Cl 沿 -x；甲氧基 O 沿 +x；乙醛 O 沿 +z
    add("Mg", 0.00, 0.00, 0.00)
    add("Cl", -2.28, 0.00, 0.00)
    add("O", 1.98, 0.00, 0.00)              # 甲氧基 O
    # 甲氧基 C：从 O 指向 +x 外侧，稍微偏上，使它的一个 H 朝向羰基碳
    add("C", 3.28, 0.00, 0.62)
    # 甲氧基的三个 H：一个朝羰基碳（迁移 H），两个朝外
    add("H", 3.86, 0.62, -0.35)             # 迁移 H（先放远）
    add("H", 3.62, -0.92, 0.62)
    add("H", 3.24, 0.28, 1.72)
    # 乙醛：O 配位到 Mg（+z），C 在 O 上方
    add("O", 0.00, 0.00, 2.05)
    add("C", 0.00, 0.00, 3.27)              # 羰基 C
    add("H", 0.10, -0.98, 3.63)             # 醛基 H
    add("C", 0.62, 0.92, 3.88)              # 甲基 C
    for d in ((1.30, 1.62, 4.30), (0.00, 1.55, 4.55), (1.28, 1.27, 3.12)):
        add("H", *d)
    return Atoms(symbols=S, positions=np.array(P))


def walk(work, i, j, target, step=0.15, fmax=0.30, maxsteps=60):
    from ase.constraints import FixBondLength as FBL
    for _ in range(120):
        cur = float(np.linalg.norm(work.positions[j] - work.positions[i]))
        gap = target - cur
        if abs(gap) < 3e-3:
            break
        d = float(np.clip(gap, -step, step))
        u = (work.positions[j] - work.positions[i]) / cur
        work.positions[i] -= u * d / 2
        work.positions[j] += u * d / 2
        work.set_constraint(FBL(i, j))
        opt = BFGS(work, logfile=None, maxstep=0.05)
        try:
            opt.run(fmax=fmax, steps=maxsteps)
        except Exception:
            pass
    return work


def main():
    calc = mace_polar(model="polar-1-m", device="cpu", default_dtype="float64",
                      pbc_handling="realspace")
    at = build()
    at.info["charge"] = 0
    at.info["spin"] = 1
    at.calc = calc
    write(os.path.join(HERE, "mpv_start.xyz"), at)

    # 先整体弛豫（不加约束），作为参考态
    opt = BFGS(at, logfile=None, maxstep=0.05)
    try:
        opt.run(fmax=0.10, steps=400)
    except Exception:
        pass
    e0 = float(at.get_potential_energy())
    write(os.path.join(HERE, "mpv_reactant.xyz"), at)
    D = at.get_all_distances(); np.fill_diagonal(D, 9)
    print("reactant relaxed: E=%.6f  min_d=%.2f  steps=%d" % (e0, D.min(), opt.get_number_of_steps()))

    H = 4            # 迁移 H 的索引
    Cmet = 3         # 甲氧基碳
    Ccar = 8         # 羰基碳
    print("d(Cmet-H)=%.2f  d(H-Ccar)=%.2f" % (
        at.get_distance(Cmet, H), at.get_distance(H, Ccar)))

    work = at.copy()
    work.calc = calc
    prof = []
    for target in np.arange(3.00, 1.04, -0.10):
        # 直接把 H 挪到目标距离，再用约束优化（每点只做一次 BFGS）
        cur = float(np.linalg.norm(work.positions[Ccar] - work.positions[H]))
        if cur < 0.6:
            break
        u = (work.positions[Ccar] - work.positions[H]) / cur
        d = float(target) - cur
        work.positions[H] -= u * d   # u 指向羰基碳，d>0 时需远离
        work.set_constraint(FixBondLength(H, Ccar))
        opt = BFGS(work, logfile=None, maxstep=0.05)
        try:
            opt.run(fmax=0.15, steps=150)
        except Exception:
            pass
        e = float(work.get_potential_energy())
        d1 = float(work.get_distance(Cmet, H))
        d2 = float(work.get_distance(H, Ccar))
        prof.append({"d_HC": round(float(target), 2), "d1": round(d1, 2), "d2": round(d2, 2),
                     "dE_kcal": float((e - e0) * KCAL)})
        print("  d(H-Ccar)=%.2f  d(Cmet-H)=%.2f  dE=%7.2f  (%d steps)" % (
            d2, d1, prof[-1]["dE_kcal"], opt.get_number_of_steps()), flush=True)

    best = max(prof, key=lambda p: p["dE_kcal"])
    json.dump({"e0": e0, "profile": prof, "barrier_kcal": best["dE_kcal"],
               "at": best}, open(os.path.join(HERE, "mpv_barrier.json"), "w"), indent=1)
    print("\nMPV barrier = %.1f kcal/mol at d(H-Ccarbonyl)=%.2f" % (
        best["dE_kcal"], best["d2"]))
    kB, h, R, T = 1.380649e-23, 6.62607015e-34, 8.314462618e-3, 298.15
    k = (kB * T / h) * np.exp(-best["dE_kcal"] / (R * T))
    print("Eyring k (per second, unimolecular reference) = %.2e 1/s" % k)


if __name__ == "__main__":
    main()
