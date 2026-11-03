#!/usr/bin/env python3
"""Build explicitly THF-solvated clusters using RDKit for the THF geometry."""
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

OUT = os.path.dirname(os.path.abspath(__file__))


def thf_template(seed=0xF00D):
    m = Chem.AddHs(Chem.MolFromSmiles("C1CCOC1"))
    AllChem.EmbedMolecule(m, randomSeed=seed, useRandomCoords=False)
    AllChem.MMFFOptimizeMolecule(m)
    conf = m.GetConformer()
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    sym = [a.GetSymbol() for a in m.GetAtoms()]
    o_idx = [i for i, s in enumerate(sym) if s == "O"][0]
    o = pos[o_idx]
    ring = [i for i, s in enumerate(sym) if s in ("C", "O")]
    centroid = pos[ring].mean(axis=0)
    v = centroid - o
    v /= np.linalg.norm(v)
    return pos - o, sym, v  # coords relative to O, plus local "into ring" direction


def rotation_from_to(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1 / (1 + c))


_T = thf_template()


def thf(origin, direction, mg_o=2.15):
    """THF whose O sits at origin+mg_o*direction, ring pointing along direction."""
    d = np.asarray(direction, float)
    d /= np.linalg.norm(d)
    coords, sym, v = _T
    R = rotation_from_to(v, d)
    pts = coords @ R.T + np.asarray(origin, float) + d * mg_o
    return [(s, p) for s, p in zip(sym, pts)]


def methyl(mg, direction, r_mg_c=2.10, r_ch=1.09, tilt=70.0):
    d = np.asarray(direction, float)
    d /= np.linalg.norm(d)
    c = np.asarray(mg, float) + d * r_mg_c
    tmp = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(d, tmp); e1 /= np.linalg.norm(e1)
    e2 = np.cross(d, e1)
    th = np.radians(tilt)
    hs = []
    for k in range(3):
        phi = 2 * np.pi * k / 3
        hd = np.sin(th) * (np.cos(phi) * e1 + np.sin(phi) * e2) + np.cos(th) * d
        hs.append(("H", c + hd * r_ch))
    return [("C", c)] + hs


def write(name, atoms):
    path = os.path.join(OUT, name + ".xyz")
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\n{name}\n")
        for s, p in atoms:
            f.write(f"{s:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")
    return path


def build(only=None):
    X, Y, Z = np.array([1., 0, 0]), np.array([0, 1., 0]), np.array([0, 0, 1.])
    o = np.zeros(3)
    out = {}

    a = [("Mg", o)] + methyl(o, X) + [("Cl", np.array([-2.25, 0, 0.]))]
    a += thf(o, Y) + thf(o, -Y)
    out["ch3mgcl_2thf"] = a

    a = [("Mg", o)] + methyl(o, X) + methyl(o, -X) + thf(o, Y) + thf(o, -Y)
    out["mgme2_2thf"] = a

    a = [("Mg", o), ("Cl", np.array([2.20, 0, 0.])), ("Cl", np.array([-2.20, 0, 0.]))]
    a += thf(o, Y) + thf(o, -Y)
    out["mgcl2_2thf"] = a

    # tetrahedral MgCl2(THF)4: two Cl along +-x, four THF in yz plane tetrahedral
    a = [("Mg", o), ("Cl", np.array([2.35, 0, 0.])), ("Cl", np.array([-2.35, 0, 0.]))]
    for ang in (45, 135, 225, 315):
        t = np.radians(ang)
        a += thf(o, np.array([0, np.cos(t), np.sin(t)]))
    out["mgcl2_4thf"] = a

    # chloro-bridged dimer: Mg2Cl2 rhombus in xy; CH3 terminal along +-x, THF along +-z
    mg1, mg2 = np.array([1.70, 0, 0.]), np.array([-1.70, 0, 0.])
    a = [("Mg", mg1), ("Mg", mg2),
         ("Cl", np.array([0, 1.25, 0.])), ("Cl", np.array([0, -1.25, 0.]))]
    a += methyl(mg1, X) + methyl(mg2, -X) + thf(mg1, Z) + thf(mg2, Z)
    out["dimer_2thf"] = a

    out["thf_free"] = thf(np.array([0, 0, -2.15]), Z)

    for k, v in out.items():
        if only and k not in only:
            continue
        write(k, v)
    print("wrote:", ", ".join(k for k in out if not only or k in only))


if __name__ == "__main__":
    import sys
    build(only=set(sys.argv[1:]) or None)
