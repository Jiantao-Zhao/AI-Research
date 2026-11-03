#!/usr/bin/env python3
"""2D relaxed scan over (C...C forming) x (Mg-C breaking) for the Agem-like adduct.

Each grid point is reached by *walking* in <=0.18 A increments, minimising after
every move, so substituents are never left behind (the failure mode of a hard
one-shot projection). Results are flushed after every point so a viewer can
follow the surface live.

Grid:  cc(C1-C20) 3.0 -> 2.0 (0.2)  x  mgc(Mg0-C1) 2.1 -> 3.3 (0.2)
Order: serpentine over rows, so consecutive points are always 0.2 A apart.
"""
import json, os, time
import numpy as np
import torch
torch.serialization.add_safe_globals([slice])
from ase.io import read, write
from ase.constraints import FixBondLength
from ase.optimize import BFGS
from mace.calculators import mace_polar

HERE = os.path.dirname(os.path.abspath(__file__))
KCAL = 23.0605
TAG = "Agem_2d"
CC = [3.0, 2.7, 2.4, 2.1]
MGC = [2.1, 2.4, 2.7, 3.0, 3.3]
I, J, MG = 1, 20, 0
PAIRS = [(I, J), (MG, I)]


def shift_pair(atoms, a, b, delta):
    v = atoms.positions[b] - atoms.positions[a]
    n = np.linalg.norm(v)
    if n < 1e-6:
        return
    u = v / n
    atoms.positions[a] -= u * delta / 2
    atoms.positions[b] += u * delta / 2


def walk(work, targets, step=0.25, fmax=0.60, maxsteps=20):
    """Approach both target distances in small increments, minimising in between."""
    for _ in range(80):
        cur = [float(np.linalg.norm(work.positions[b] - work.positions[a])) for a, b in PAIRS]
        gaps = [t - c for c, t in zip(cur, targets)]
        if max(abs(g) for g in gaps) < 5e-3:
            break
        for (a, b), g in zip(PAIRS, gaps):
            shift_pair(work, a, b, float(np.clip(g, -step, step)))
        work.set_constraint([FixBondLength(a, b) for a, b in PAIRS])
        opt = BFGS(work, logfile=None, maxstep=0.05)
        try:
            opt.run(fmax=fmax, steps=maxsteps)
        except Exception:
            pass
    return work


def dump(state):
    tmp = os.path.join(HERE, "scan2d_partial.json.tmp")
    with open(tmp, "w") as f:
        json.dump(state, f)
    os.replace(tmp, os.path.join(HERE, "scan2d_partial.json"))


def main():
    calc = mace_polar(model="polar-1-m", device="cpu", default_dtype="float64",
                      pbc_handling="realspace")
    seed = read(os.path.join(HERE, "Agem_like_min.xyz"))
    seed.info["charge"] = 0
    seed.info["spin"] = 1
    seed.calc = calc
    e0 = float(seed.get_potential_energy())

    geo = os.path.join(HERE, "geoms_2d")
    os.makedirs(geo, exist_ok=True)
    state = {"tag": TAG, "e0": e0, "cc": CC, "mgc": MGC,
             "total": len(CC) * len(MGC), "done": 0,
             "started": time.strftime("%Y-%m-%d %H:%M:%S"),
             "points": [], "running": None, "finished": False}
    dump(state)

    work = seed.copy()
    work.calc = calc
    t_start = time.time()

    for r, mgc in enumerate(MGC):
        row = CC if r % 2 == 0 else CC[::-1]
        for cc in row:
            state["running"] = {"cc": cc, "mgc": mgc}
            dump(state)
            t0 = time.time()
            walk(work, [cc, mgc])
            work.set_constraint([FixBondLength(a, b) for a, b in PAIRS])   # 保持约束！
            opt = BFGS(work, logfile=None, maxstep=0.08)
            ok = True
            try:
                opt.run(fmax=0.10, steps=150)
            except Exception:
                ok = False
            e = float(work.get_potential_energy())
            D = work.get_all_distances(); np.fill_diagonal(D, 9.0)
            mind = float(D.min())
            cc_act = float(work.get_distance(I, J))
            mgc_act = float(work.get_distance(MG, I))
            if mind < 0.80 or abs(cc_act - cc) > 0.06 or abs(mgc_act - mgc) > 0.06:
                ok = False
            fn = "%s_c%.1f_m%.1f.xyz" % (TAG, cc, mgc)
            write(os.path.join(geo, fn), work)
            state["points"].append({
                "cc": cc, "mgc": mgc, "E_eV": e,
                "dE_kcal": float((e - e0) * KCAL),
                "steps": int(opt.get_number_of_steps()), "min_d": round(mind, 2),
                "cc_act": round(cc_act, 2), "mgc_act": round(mgc_act, 2),
                "ok": ok, "xyz": fn, "wall": round(time.time() - t0, 1)})
            state["done"] += 1
            dump(state)
            print("[%2d/%d] cc=%.1f(%.2f) mgc=%.1f(%.2f) dE=%8.2f min_d=%.2f ok=%s (%ss)" % (
                state["done"], state["total"], cc, cc_act, mgc, mgc_act,
                state["points"][-1]["dE_kcal"], mind, ok, state["points"][-1]["wall"]), flush=True)

    state["finished"] = True
    state["running"] = None
    state["wall_total"] = round(time.time() - t_start, 1)
    dump(state)
    print("DONE in %s s" % state["wall_total"], flush=True)


if __name__ == "__main__":
    main()
