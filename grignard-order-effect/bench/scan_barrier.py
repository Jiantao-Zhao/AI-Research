#!/usr/bin/env python3
"""Zero-shot MACE-POLAR barrier test for Grignard addition to acetaldehyde.

Reconstruction: the JACS 2020 SI contains no coordinates, so we rebuild the two
mononuclear "geminal" complexes using the protocol stated in that paper
("one THF coordinated to Mg is replaced by a substrate molecule").

  Agem-like : CH3MgCl(THF)(CH3CHO)   [literature analogue: 8.7 +/- 0.8 kcal/mol]
  Bgem-like : Mg(CH3)2(THF)(CH3CHO)  [literature analogue: 6.5 +/- 0.7 kcal/mol]

Barrier = relaxed scan along the forming C(methyl)-C(carbonyl) bond, electronic
energy only (no ZPE/entropy), so numbers are NOT directly comparable to the
reported Helmholtz DeltaA(dagger); they are an upper-bound-ish proxy.
"""
import json, os, sys, time
import warnings
warnings.filterwarnings("ignore")
import numpy as np

import torch
torch.serialization.add_safe_globals([slice])
from ase import Atoms
from ase.io import read, write
from ase.constraints import FixBondLength
from ase.optimize import BFGS
from mace.calculators import mace_polar

HERE = os.path.dirname(os.path.abspath(__file__))
XTB = os.path.join(HERE, "..", "model", "xtb")
KCAL = 23.0605


def acetaldehyde(seed=0xF00D):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles("CC=O"))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m)
    c = m.GetConformer()
    pos = np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    sym = [a.GetSymbol() for a in m.GetAtoms()]
    # order: put O first, then the carbonyl C (the C bonded to O)
    oi = [i for i, s in enumerate(sym) if s == "O"][0]
    ci = [a.GetIdx() for a in m.GetAtomWithIdx(oi).GetNeighbors()][0]
    order = [oi, ci] + [i for i in range(len(sym)) if i not in (oi, ci)]
    return np.array(sym)[order].tolist(), pos[order]


def orient_aldehyde(pos, sym, o_pos, out_dir, mg_o=2.12, co=1.22):
    """Place the aldehyde so O sits at o_pos and the C=O axis points along out_dir."""
    p = pos.copy()
    o_local = p[0]
    c_local = p[1]
    v = c_local - o_local
    v /= np.linalg.norm(v)
    d = np.asarray(out_dir, float)
    d /= np.linalg.norm(d)
    # rotation taking v -> d
    a = np.cross(v, d)
    s = np.linalg.norm(a)
    cth = float(np.dot(v, d))
    if s < 1e-8:
        R = np.eye(3) if cth > 0 else -np.eye(3)
    else:
        a /= s
        K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
        R = np.eye(3) + s * K + (1 - cth) * K @ K
    p = (p - o_local) @ R.T
    return p + np.asarray(o_pos, float)


def thf_indices(atoms, mg_idx=0, cut=2.45):
    """Whole THF molecules coordinated to Mg, found by graph traversal."""
    from ase.data import covalent_radii, atomic_numbers
    D = atoms.get_all_distances()
    z = [atomic_numbers[a.symbol] for a in atoms]
    n = len(atoms)
    adj = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] < 1.25 * (covalent_radii[z[i]] + covalent_radii[z[j]]):
                adj[i].append(j)
                adj[j].append(i)
    o_idx = [i for i in range(n) if atoms[i].symbol == "O" and D[mg_idx, i] < cut]
    out = []
    for o in o_idx:
        seen, stack = {o}, [o]
        while stack:
            a = stack.pop()
            for b in adj[a]:
                if atoms[b].symbol == "Mg":   # do not walk through Mg into the rest
                    continue
                if b not in seen:
                    seen.add(b)
                    stack.append(b)
        out.append(sorted(seen))
    return o_idx, out


def build_complex(src, which_thf=0):
    at = read(src)
    o_idx, groups = thf_indices(at)
    if len(groups) < 2:
        raise RuntimeError(f"{src}: found {len(groups)} THF, expected 2")
    drop = set(groups[which_thf])
    o_pos = at[o_idx[which_thf]].position
    mg = at[0].position
    out_dir = o_pos - mg
    keep = [i for i in range(len(at)) if i not in drop]
    sym = [at[i].symbol for i in keep]
    pos = np.array([at[i].position for i in keep])
    asym, apos = acetaldehyde()
    apos = orient_aldehyde(apos, asym, np.array(o_pos), out_dir)
    sym += asym
    pos = np.vstack([pos, apos])
    new = Atoms(symbols=sym, positions=pos)
    # index bookkeeping: methyl C = the C bonded to Mg (first C in the original)
    mg_i = 0
    D = new.get_all_distances()
    c_methyl = [i for i in range(len(new)) if new[i].symbol == "C" and D[0, i] < 2.35][0]
    c_carbonyl = len(keep) + 1  # O is first of the appended block
    return new, c_methyl, c_carbonyl


def set_distance(atoms, i, j, r):
    v = atoms.positions[j] - atoms.positions[i]
    n = np.linalg.norm(v)
    if n < 1e-9:
        raise RuntimeError("zero bond vector")
    shift = (r - n) / 2.0
    u = v / n
    atoms.positions[i] -= u * shift
    atoms.positions[j] += u * shift


def scan(calc, atoms, i, j, distances, fmax=0.08, maxsteps=250, tag="",
         geo_dir=None, resume=False):
    """Relaxed scan along i-j. Every point is written to geo_dir as extxyz."""
    atoms.info["charge"] = 0
    atoms.info["spin"] = 1
    atoms.calc = calc
    e0 = atoms.get_potential_energy()
    prof = []
    work = atoms.copy()
    work.calc = calc
    if geo_dir:
        os.makedirs(geo_dir, exist_ok=True)
    for r in distances:
        geo = os.path.join(geo_dir, f"{tag}_r{r:.2f}.xyz") if geo_dir else None
        if resume and geo and os.path.exists(geo):
            work = read(geo)
            work.info["charge"] = 0
            work.info["spin"] = 1
            work.calc = calc
            e = float(work.get_potential_energy())
            r_mgc = float(work.get_distance(0, i))
            prof.append({"r": round(float(r), 3), "E_eV": float(e),
                         "dE_kcal": float((e - e0) * KCAL), "r_MgC": round(r_mgc, 3),
                         "steps": 0, "xyz": os.path.basename(geo), "cached": True})
            print(f"    [{tag}] r={r:.2f} A  dE={prof[-1]['dE_kcal']:7.2f} kcal/mol "
                  f"(cached)", flush=True)
            continue
        set_distance(work, i, j, r)
        work.set_constraint(FixBondLength(i, j))
        opt = BFGS(work, logfile=None, maxstep=0.12)
        try:
            opt.run(fmax=fmax, steps=maxsteps)
        except Exception as exc:
            print(f"    [{tag}] r={r}: opt failed ({exc})", flush=True)
        e = work.get_potential_energy()
        r_mgc = float(work.get_distance(0, i))
        work.set_constraint()
        if geo:
            write(geo, work)
        prof.append({"r": round(float(r), 3), "E_eV": float(e),
                     "dE_kcal": float((e - e0) * KCAL),
                     "r_MgC": round(r_mgc, 3),
                     "steps": int(opt.get_number_of_steps()),
                     "xyz": os.path.basename(geo) if geo else None,
                     "cached": False})
        print(f"    [{tag}] r={r:.2f} A  dE={prof[-1]['dE_kcal']:7.2f} kcal/mol "
              f"r(Mg-C)={r_mgc:.2f} A ({opt.get_number_of_steps()} steps)", flush=True)
    work.set_constraint()
    return prof


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=HERE, help="where scan_results.json goes")
    ap.add_argument("--geoms", default=None, help="dir for per-point extxyz (default OUTDIR/geoms)")
    ap.add_argument("--max-points", type=int, default=0, help="limit points (debug)")
    ap.add_argument("--resume", action="store_true", help="reuse existing point files")
    ap.add_argument("--only", default=None, help="comma-separated species names")
    args = ap.parse_args()

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    geo_dir = os.path.abspath(args.geoms) if args.geoms else os.path.join(outdir, "geoms")

    calc = mace_polar(model="polar-1-m", device="cpu", default_dtype="float64",
                      pbc_handling="realspace")
    results = {}
    cases = [
        ("Agem_like", os.path.join(XTB, "ch3mgcl_2thf_opt.xyz"), 8.7),
        ("Bgem_like", os.path.join(XTB, "mgme2_2thf_opt.xyz"), 6.5),
    ]
    only = set(args.only.split(",")) if args.only else None
    for name, src, lit in cases:
        if only and name not in only:
            continue
        t0 = time.time()
        print(f"\n=== {name}  (literature analogue {lit} kcal/mol) ===", flush=True)
        atoms, ci, cj = build_complex(src)
        print(f"  built: {len(atoms)} atoms; forming bond C{ci}-C{cj}", flush=True)
        atoms.info["charge"] = 0
        atoms.info["spin"] = 1
        atoms.calc = calc
        opt = BFGS(atoms, logfile=None, maxstep=0.10)
        opt.run(fmax=0.08, steps=400)
        print(f"  precursor relaxed in {opt.get_number_of_steps()} steps", flush=True)
        write(os.path.join(outdir, f"{name}_min.xyz"), atoms)
        dists = list(np.round(np.arange(4.40, 1.79, -0.20), 2))
        if args.max_points:
            dists = dists[:args.max_points]
        prof = scan(calc, atoms, ci, cj, dists, tag=name,
                    geo_dir=geo_dir, resume=args.resume)
        emax = max(prof, key=lambda p: p["dE_kcal"])
        results[name] = {
            "literature_kcal": lit,
            "max_dE_kcal": emax["dE_kcal"],
            "r_at_max": emax["r"],
            "profile": prof,
            "wall_s": round(time.time() - t0, 1),
            "n_atoms": len(atoms),
            "n_points": len(prof),
            "geoms_dir": geo_dir,
        }
        print(f"  -> barrier (DeltaE, relaxed scan) = {emax['dE_kcal']:.2f} kcal/mol "
              f"at r = {emax['r']} A   [{results[name]['wall_s']} s]", flush=True)
        json.dump(results, open(os.path.join(outdir, "scan_results.json"), "w"), indent=1)
    print("\nDONE")
