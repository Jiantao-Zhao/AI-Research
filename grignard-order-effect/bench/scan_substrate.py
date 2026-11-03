#!/usr/bin/env python3
"""对任意醛/酮重复 Agem-like 的加成势垒扫描。

用法: scan_substrate.py --smiles "O=Cc1ccccc1" --name BzCHO --tag Agem_Bz
"""
import argparse, json, os, sys, time
import numpy as np
import torch
torch.serialization.add_safe_globals([slice])
from ase.io import read, write
from ase.constraints import FixBondLength
from ase.optimize import BFGS
from mace.calculators import mace_polar

HERE = os.path.dirname(os.path.abspath(__file__))
XTB = os.path.join(HERE, "..", "model", "xtb")
KCAL = 23.0605


def mol_from_smiles(smiles, seed=0xF00D):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m)
    c = m.GetConformer()
    pos = np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    sym = [a.GetSymbol() for a in m.GetAtoms()]
    oi = [i for i, s in enumerate(sym) if s == "O"][0]
    ci = [a.GetIdx() for a in m.GetAtomWithIdx(oi).GetNeighbors()][0]
    order = [oi, ci] + [i for i in range(len(sym)) if i not in (oi, ci)]
    return list(np.array(sym)[order]), pos[order]


def thf_groups(atoms, mg=0, cut=2.45):
    from ase.data import covalent_radii, atomic_numbers
    D = atoms.get_all_distances(); z = [atomic_numbers[a.symbol] for a in atoms]
    n = len(atoms); adj = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] < 1.25 * (covalent_radii[z[i]] + covalent_radii[z[j]]):
                adj[i].append(j); adj[j].append(i)
    outs = []
    for o in [i for i in range(n) if atoms[i].symbol == "O" and D[mg, i] < cut]:
        seen, st = {o}, [o]
        while st:
            a = st.pop()
            for b in adj[a]:
                if atoms[b].symbol == "Mg": continue
                if b not in seen: seen.add(b); st.append(b)
        outs.append(sorted(seen))
    return outs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--model", default="polar-1-m")
    a = ap.parse_args()
    tag = a.tag or a.name

    calc = mace_polar(model=a.model, device="cpu", default_dtype="float64",
                      pbc_handling="realspace")
    src = read(os.path.join(XTB, "ch3mgcl_2thf_opt.xyz"))
    groups = thf_groups(src)
    drop = set(groups[0])
    o_pos = np.array([src[i].position for i in groups[0] if src[i].symbol == "O"][0])
    mg = src[0].position
    keep = [i for i in range(len(src)) if i not in drop]
    sym = [src[i].symbol for i in keep]
    pos = np.array([src[i].position for i in keep])

    asym, apos = mol_from_smiles(a.smiles)
    o0, c0 = apos[0], apos[1]
    v = c0 - o0; v /= np.linalg.norm(v)
    d = (o_pos - mg); d /= np.linalg.norm(d)
    aa = np.cross(v, d); s = np.linalg.norm(aa); c = float(np.dot(v, d))
    if s < 1e-8:
        R = np.eye(3) if c > 0 else -np.eye(3)
    else:
        aa /= s
        K = np.array([[0, -aa[2], aa[1]], [aa[2], 0, -aa[0]], [-aa[1], aa[0], 0]])
        R = np.eye(3) + s * K + (1 - c) * K @ K
    apos = (apos - o0) @ R.T + o_pos
    atoms = __import__("ase").Atoms(symbols=sym + asym, positions=np.vstack([pos, apos]))
    atoms.info["charge"] = 0; atoms.info["spin"] = 1; atoms.calc = calc

    i = 1; j = len(keep) + 1
    opt = BFGS(atoms, logfile=None, maxstep=0.10)
    try: opt.run(fmax=0.10, steps=400)
    except Exception: pass
    e0 = float(atoms.get_potential_energy())
    write(os.path.join(HERE, f"{tag}_min.xyz"), atoms)
    print("built %d atoms; relaxing done; C%d-C%d" % (len(atoms), i, j), flush=True)

    work = atoms.copy(); work.calc = calc
    prof = []
    for r in np.arange(4.40, 1.79, -0.20):
        cur = float(work.get_distance(i, j))
        u = (work.positions[j] - work.positions[i]) / cur
        work.positions[i] -= u * (float(r) - cur) / 2
        work.positions[j] += u * (float(r) - cur) / 2
        work.set_constraint(FixBondLength(i, j))
        o = BFGS(work, logfile=None, maxstep=0.08)
        try: o.run(fmax=0.12, steps=150)
        except Exception: pass
        e = float(work.get_potential_energy())
        work.set_constraint()
        prof.append({"r": round(float(r), 2), "dE_kcal": float((e - e0) * KCAL),
                     "r_MgC": round(float(work.get_distance(0, i)), 2)})
        print("  r=%.2f dE=%7.2f" % (r, prof[-1]["dE_kcal"]), flush=True)

    best = max(prof, key=lambda p: p["dE_kcal"])
    json.dump({"tag": tag, "name": a.name, "smiles": a.smiles, "e0": e0,
               "profile": prof, "barrier_kcal": best["dE_kcal"], "r": best["r"]},
              open(os.path.join(HERE, f"{tag}_barrier.json"), "w"), indent=1)
    print("\n%s: ΔE‡ = %.2f kcal/mol at r=%.2f" % (a.name, best["dE_kcal"], best["r"]))


if __name__ == "__main__":
    main()
