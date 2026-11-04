"""多卡并行伞形采样：一个 rank 分到一组窗口。

用法：
  srun --gpus-per-task=1 -n 4 python umbrella_rank.py --rank $SLURM_PROCID --nranks 4
"""
import argparse, json, os, time
import numpy as np
from ase.io import read
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.optimize import BFGS
from ase import units
from mace.calculators import mace_polar

C1, C20 = 1, 20


def spring(r0):
    """内层用硬弹簧，外层用软弹簧。"""
    if r0 <= 2.60:
        return 20.0
    if r0 <= 3.20:
        return 5.0
    return 2.0


class Umbrella:
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank", type=int, default=0)
    ap.add_argument("--nranks", type=int, default=1)
    ap.add_argument("--dir", default="run")
    ap.add_argument("--rmin", type=float, default=1.90)
    ap.add_argument("--rmax", type=float, default=3.90)
    ap.add_argument("--dr", type=float, default=0.167)
    ap.add_argument("--ps", type=float, default=3.0)
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--stride", type=int, default=4)
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()

    umb = os.path.join(a.dir, "umb")
    os.makedirs(umb, exist_ok=True)
    at0 = read(os.path.join(a.dir, "box_equil.xyz"))
    at0.info["charge"] = 0; at0.info["spin"] = 1
    box = float(at0.cell[0, 0])
    calc = mace_polar(model="polar-1-m", device=a.device, default_dtype="float32",
                      pbc_handling="pbc")
    windows = np.arange(a.rmin, a.rmax + 1e-9, a.dr)

    for w, r0 in enumerate(windows):
        if w % a.nranks != a.rank:
            continue
        dfin = os.path.join(umb, f"d_{w:02d}.npy")
        if os.path.exists(dfin):
            print(f"[r{a.rank}] win {w} 已存在，跳过", flush=True)
            continue
        k = spring(r0)
        at = at0.copy()
        at.info["charge"] = 0; at.info["spin"] = 1
        d = at.get_distance(C1, C20)
        u = (at.positions[C20] - at.positions[C1]) / d
        shift = np.clip(r0 - d, -1.2, 1.2)
        at.positions[C1] -= u * shift / 2
        at.positions[C20] += u * shift / 2
        at.calc = Umbrella(calc, C1, C20, r0, k, box)

        o = BFGS(at, logfile=None, maxstep=0.04)
        try:
            o.run(fmax=0.35, steps=100)
        except Exception:
            pass
        MaxwellBoltzmannDistribution(at, temperature_K=a.T)
        Stationary(at)
        dyn = Langevin(at, 0.5 * units.fs, temperature_K=a.T, friction=0.01)
        n = int(a.ps * 1000 / 0.5)
        ds = []
        for i in range(n):
            dyn.run(1)
            if i % a.stride == 0 and i > n // 5:
                ds.append(float(at.get_distance(C1, C20)))
        np.save(dfin, np.array(ds))
        json.dump({"w": w, "r0": float(r0), "k": k, "n": len(ds),
                   "mean": float(np.mean(ds)), "std": float(np.std(ds))},
                  open(os.path.join(umb, f"w_{w:02d}.json"), "w"), indent=1)
        print(f"[r{a.rank}] win {w} r0={r0:.2f} k={k:.0f} <d>={np.mean(ds):.3f} "
              f"sd={np.std(ds):.3f}", flush=True)


if __name__ == "__main__":
    main()
