#!/usr/bin/env python3
"""Build starting geometries for the CH3MgCl Schlenk system (no external deps beyond numpy)."""
import os
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)


def methyl_on(anchor, axis, r_mg_c=2.10, r_ch=1.09, tilt_deg=70.0):
    """Place CH3 carbon at anchor + axis*r_mg_c, three H's tilted off the axis."""
    axis = np.asarray(axis, float)
    axis /= np.linalg.norm(axis)
    c = np.asarray(anchor, float) + axis * r_mg_c
    # orthonormal frame
    tmp = np.array([0.0, 0.0, 1.0]) if abs(axis[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(axis, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    th = np.radians(tilt_deg)
    hs = []
    for k in range(3):
        phi = 2 * np.pi * k / 3
        d = np.sin(th) * (np.cos(phi) * e1 + np.sin(phi) * e2) + np.cos(th) * axis
        hs.append(c + d * r_ch)
    return c, hs


def write_xyz(name, atoms, comment):
    path = os.path.join(OUT, name + ".xyz")
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for sym, xyz in atoms:
            f.write(f"{sym:2s} {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f}\n")
    return path


def monomer_ch3mgcl():
    atoms = [("Mg", np.zeros(3))]
    c, hs = methyl_on(np.zeros(3), [1, 0, 0])
    atoms.append(("C", c))
    atoms += [("H", h) for h in hs]
    atoms.append(("Cl", np.array([-2.25, 0.0, 0.0])))
    return atoms


def dimethyl_mg():
    atoms = [("Mg", np.zeros(3))]
    for ax in ([1, 0, 0], [-1, 0, 0]):
        c, hs = methyl_on(np.zeros(3), ax)
        atoms.append(("C", c))
        atoms += [("H", h) for h in hs]
    return atoms


def mgcl2():
    return [("Mg", np.zeros(3)), ("Cl", np.array([2.20, 0, 0])),
            ("Cl", np.array([-2.20, 0, 0]))]


def dimer_bridged():
    """(CH3MgCl)2 with two bridging Cl: planar Mg2Cl2 rhombus, terminal CH3."""
    atoms = []
    mg1 = np.array([1.70, 0.0, 0.0])
    mg2 = np.array([-1.70, 0.0, 0.0])
    cl1 = np.array([0.0, 1.25, 0.0])
    cl2 = np.array([0.0, -1.25, 0.0])
    atoms += [("Mg", mg1), ("Mg", mg2), ("Cl", cl1), ("Cl", cl2)]
    for mg, sgn in ((mg1, 1), (mg2, -1)):
        c, hs = methyl_on(mg, [sgn, 0, 0])
        atoms.append(("C", c))
        atoms += [("H", h) for h in hs]
    return atoms


if __name__ == "__main__":
    specs = {
        "ch3mgcl_monomer": monomer_ch3mgcl(),
        "mgch32_monomer": dimethyl_mg(),
        "mgcl2_monomer": mgcl2(),
        "ch3mgcl_dimer": dimer_bridged(),
    }
    for name, atoms in specs.items():
        p = write_xyz(name, atoms, name)
        print("wrote", p, len(atoms), "atoms")
