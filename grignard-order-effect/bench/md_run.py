"""显式 THF 溶剂中的加成自由能面（MACE-POLAR + 伞形采样 + WHAM）。

阶段：
  compress  低密度建盒 → 逐步压缩 → 最小化
  equil     NVT 300 K 平衡
  umbrella  沿 C1–C20 的伞形采样
  wham      一维 WHAM → PMF
"""
import argparse, os, sys, time, json
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.optimize import BFGS
from ase import units
from mace.calculators import mace_polar
import md_build_box as mb

C1, C20 = 1, 20
BOX0, NTHF = 18.0, 22


class Umbrella:
    """谐波约束 d(C1,C20) 的包装器。"""

    def __init__(self, calc, i, j, r0, k, box):
        self.calc, self.i, self.j, self.r0, self.k, self.box = calc, i, j, r0, k, box

    def _vec(self, atoms):
        d = atoms.positions[self.j] - atoms.positions[self.i]
        return d - self.box * np.round(d / self.box)

    def get_potential_energy(self, atoms=None, **kw):
        e = self.calc.get_potential_energy(atoms)
        d = np.linalg.norm(self._vec(atoms))
        return e + 0.5 * self.k * (d - self.r0) ** 2

    def get_forces(self, atoms=None, **kw):
        f = np.array(self.calc.get_forces(atoms))
        v = self._vec(atoms); d = np.linalg.norm(v); u = v / d
        g = self.k * (d - self.r0) * u
        f[self.i] += g
        f[self.j] -= g
        return f

    def get_potential_energies(self, atoms=None, **kw):
        return self.calc.get_potential_energies(atoms)

    def reset(self):
        if hasattr(self.calc, "reset"):
            self.calc.reset()


def get_calc(device="cuda"):
    return mace_polar(model="polar-1-m", device=device, default_dtype="float32",
                      pbc_handling="pbc")


def stage_compress(workdir, device, target=15.2, nsteps=22):
    at = mb.build(n_thf=NTHF, box=BOX0, out=os.path.join(workdir, "box_start.xyz"))
    calc = get_calc(device)
    at.calc = calc
    s = (target / BOX0) ** (1.0 / nsteps)
    for k in range(nsteps):
        at.set_cell(at.cell * s, scale_atoms=True)
        o = BFGS(at, logfile=None, maxstep=0.08)
        o.run(fmax=0.35, steps=60)
        print(f"  [compress {k+1}/{nsteps}] L={at.cell[0,0]:.2f} "
              f"E={at.get_potential_energy():.3f} eV", flush=True)
    at.set_constraint(None)
    o = BFGS(at, logfile=None, maxstep=0.05)
    o.run(fmax=0.12, steps=300)
    write(os.path.join(workdir, "box_min.xyz"), at)
    rho = NTHF * 72 / 6.022e23 / (at.cell[0, 0] ** 3 * 1e-24)
    print(f"最终盒长 {at.cell[0,0]:.2f} Å，密度 {rho:.2f} g/cm3", flush=True)
    return at


def stage_equil(workdir, device, ps=10.0, T=300.0):
    at = read(os.path.join(workdir, "box_min.xyz"))
    at.info["charge"] = 0; at.info["spin"] = 1
    at.calc = get_calc(device)
    MaxwellBoltzmannDistribution(at, temperature_K=T)
    Stationary(at)
    dt = 0.5 * units.fs
    dyn = Langevin(at, dt, temperature_K=T, friction=0.01)
    n = int(ps * 1000 / 0.5)
    t0 = time.time()
    for i in range(n):
        dyn.run(1)
        if (i + 1) % 4000 == 0:
            T_inst = at.get_kinetic_energy() / (1.5 * units.kB * len(at))
            print(f"  [equil {(i+1)*0.5/1000:.0f} ps] T={T_inst:.0f} K "
                  f"E={at.get_potential_energy():.3f} eV", flush=True)
    write(os.path.join(workdir, "box_equil.xyz"), at)
    print(f"平衡完成，耗时 {time.time()-t0:.0f} s", flush=True)
    return at


def stage_umbrella(workdir, device, rmin=1.85, rmax=3.95, dr=0.19, ps=3.0,
                   k=2.0, T=300.0, stride=4):
    at0 = read(os.path.join(workdir, "box_equil.xyz"))
    at0.info["charge"] = 0; at0.info["spin"] = 1
    box = float(at0.cell[0, 0])
    base = get_calc(device)
    windows = np.arange(rmin, rmax + 1e-9, dr)
    out = {"windows": [], "k": k, "T": T, "stride": stride}
    os.makedirs(os.path.join(workdir, "umb"), exist_ok=True)
    t0 = time.time()
    for w, r0 in enumerate(windows):
        at = at0.copy()
        at.info["charge"] = 0; at.info["spin"] = 1
        # 把 C1–C20 拉到窗口中心（小幅）
        d = at.get_distance(C1, C20)
        u = (at.positions[C20] - at.positions[C1]) / d
        at.positions[C1] -= u * (r0 - d) / 2
        at.positions[C20] += u * (r0 - d) / 2
        at.calc = Umbrella(base, C1, C20, r0, k, box)
        at.set_constraint(None)
        o = BFGS(at, logfile=None, maxstep=0.05)
        try:
            o.run(fmax=0.4, steps=80)
        except Exception:
            pass
        MaxwellBoltzmannDistribution(at, temperature_K=T)
        Stationary(at)
        dyn = Langevin(at, 0.5 * units.fs, temperature_K=T, friction=0.01)
        n = int(ps * 1000 / 0.5)
        ds = []
        for i in range(n):
            dyn.run(1)
            if i % stride == 0 and i > n // 4:
                ds.append(float(at.get_distance(C1, C20)))
        np.save(os.path.join(workdir, "umb", f"d_{w:02d}.npy"), np.array(ds))
        out["windows"].append({"w": w, "r0": float(r0), "n": len(ds),
                               "mean": float(np.mean(ds)), "std": float(np.std(ds))})
        print(f"  [umb {w+1}/{len(windows)}] r0={r0:.2f} <d>={np.mean(ds):.3f} "
              f"sd={np.std(ds):.3f}", flush=True)
    out["wall_s"] = time.time() - t0
    json.dump(out, open(os.path.join(workdir, "umb", "meta.json"), "w"), indent=1)
    print(f"伞形采样完成，耗时 {out['wall_s']:.0f} s", flush=True)


def wham(workdir, k, T, nbins=120, tol=1e-8, maxiter=20000):
    meta = json.load(open(os.path.join(workdir, "umb", "meta.json")))
    r0s = np.array([w["r0"] for w in meta["windows"]])
    data = [np.load(os.path.join(workdir, "umb", f"d_{w['w']:02d}.npy")) for w in meta["windows"]]
    allx = np.concatenate(data)
    lo, hi = allx.min(), allx.max()
    edges = np.linspace(lo, hi, nbins + 1)
    hist = np.array([np.histogram(x, bins=edges)[0] for x in data], float)
    N = hist.sum(axis=1)
    centers = 0.5 * (edges[1:] + edges[:-1])
    kT = T * units.kB  # eV
    bias = 0.5 * k * (centers[None, :] - r0s[:, None]) ** 2
    f = np.zeros(len(r0s))
    nb = hist.sum(axis=0)
    for it in range(maxiter):
        denom = (N[:, None] * np.exp(f[:, None] - bias) / kT).sum(axis=0)
        denom = np.where(nb > 0, denom, 1.0)
        newf = -kT * np.log(np.where(nb > 0, nb / denom, 1.0))
        newf -= newf[0]
        if np.max(np.abs(newf - f)) < tol:
            f = newf; break
        f = newf
    pmf = -kT * np.log(np.where(nb > 0, nb / np.maximum(
        (N[:, None] * np.exp(f[:, None] - bias) / kT).sum(axis=0), 1e-300), 1.0))
    pmf -= pmf[nb > 0].min()
    pmf_kcal = pmf * 23.0605
    return centers, pmf_kcal, nb


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["compress", "equil", "umbrella", "wham", "all"])
    ap.add_argument("--dir", default="run")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--ps", type=float, default=8.0)
    a = ap.parse_args()
    os.makedirs(a.dir, exist_ok=True)
    os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..") if False else None
    if a.stage in ("compress", "all"):
        stage_compress(a.dir, a.device)
    if a.stage in ("equil", "all"):
        stage_equil(a.dir, a.device)
    if a.stage in ("umbrella", "all"):
        stage_umbrella(a.dir, a.device, ps=a.ps)
    if a.stage in ("wham", "all"):
        c, p, nb = wham(a.dir, 2.0, 300.0)
        json.dump({"r": c.tolist(), "pmf": p.tolist(), "counts": nb.tolist()},
                  open(os.path.join(a.dir, "pmf.json"), "w"))
        print("PMF 已写出")
