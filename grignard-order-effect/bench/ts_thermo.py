#!/usr/bin/env python3
"""Step 1+2: refine the barrier, verify it is a true saddle point, add thermochemistry.

1. fine relaxed scan (0.05 A) around the coarse maximum
2. Dimer saddle-point optimisation from the scan maximum -> true TS
3. numerical Hessian: confirm exactly one imaginary frequency
4. RRHO and quasi-RRHO free-energy corrections at 298.15 K

Everything at the MACE-POLAR-1-M level, gas phase, single conformation.
"""
import json, os, sys, time
import warnings
warnings.filterwarnings("ignore")
import numpy as np

import torch
torch.serialization.add_safe_globals([slice])
from ase.io import read, write
from ase.constraints import FixBondLength
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from mace.calculators import mace_polar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scan_barrier as sb

HERE = os.path.dirname(os.path.abspath(__file__))
KCAL = 23.0605
T = 298.15
KB_EV = 8.617333262e-5
H_EV_S = 4.135667696e-15
C_CM = 2.99792458e10


def fine_scan(calc, atoms, i, j, dists, fmax=0.05, maxsteps=300, tag=""):
    work = atoms.copy()
    work.calc = calc
    work.info["charge"], work.info["spin"] = 0, 1
    e0 = work.get_potential_energy()
    prof, geoms = [], []
    for r in dists:
        sb.set_distance(work, i, j, r)
        work.set_constraint(FixBondLength(i, j))
        opt = BFGS(work, logfile=None, maxstep=0.10)
        try:
            opt.run(fmax=fmax, steps=maxsteps)
        except Exception as exc:
            print(f"    [{tag}] r={r}: {exc}", flush=True)
        e = work.get_potential_energy()
        prof.append({"r": round(float(r), 3), "dE_kcal": float((e - e0) * KCAL),
                     "r_MgC": round(float(work.get_distance(0, i)), 3)})
        g = work.copy()
        g.set_constraint()
        geoms.append(g)
        print(f"    [{tag}] r={r:.2f}  dE={prof[-1]['dE_kcal']:6.2f}  "
              f"r(Mg-C)={prof[-1]['r_MgC']:.2f}", flush=True)
    work.set_constraint()
    return prof, geoms, e0


def dimer_ts(atoms, vec, fmax=0.02, steps=400, dimer_log=None):
    from ase.mep.dimer import DimerControl, MinModeAtoms, MinModeTranslate
    dc = DimerControl(initial_eigenmode_method="displacement",
                      displacement_method="vector",
                      maximum_translation=0.05, logfile=None)
    d_at = MinModeAtoms(atoms, dc)
    d_at.displace(displacement_vector=vec)
    opt = MinModeTranslate(d_at, logfile=f"{dimer_log}")
    opt.run(fmax=fmax, steps=steps)
    return d_at.atoms.copy(), int(opt.get_number_of_steps())


def freqs(atoms, name):
    """返回复数频率（cm^-1）。虚频是复数的虚部，不能取 .real。

    ASE 的 get_frequencies() 对负本征值返回纯虚数；`.real` 会把虚频变成 0。
    """
    at = atoms.copy()
    at.calc = atoms.calc          # Atoms.copy() drops the calculator
    if at.calc is None:
        raise RuntimeError(f"{name}: no calculator attached")
    at.pbc = False
    at.set_cell([0, 0, 0])
    vib = Vibrations(at, name=name, delta=0.01)
    vib.run()
    f = np.asarray(vib.get_frequencies(), dtype=complex)
    vib.clean()
    return f[np.argsort(f.real)]


def split_modes(f, thresh=1.0):
    """把复数频率拆成 (真实振动模, 虚频幅值)。"""
    imag_mask = np.abs(f.imag) > thresh
    real_modes = f.real[~imag_mask]
    imag_mag = np.abs(f.imag[imag_mask])
    return np.sort(real_modes), np.sort(imag_mag)[::-1]


def rrho_thermo(freq_cm, scale_low=None):
    """Vibrational free-energy correction (eV) at 1 atm, ideal gas, nonlinear.

    freq_cm : array, wavenumbers (imaginary modes must be excluded/zeroed)
    scale_low : if not None, all frequencies below this value are replaced by it
                (simple quasi-harmonic treatment)
    """
    from ase import Atoms
    f = np.array([x for x in freq_cm if x > 0], dtype=float)
    if scale_low is not None:
        f = np.where(f < scale_low, scale_low, f)
    kT = KB_EV * T
    # ZPE
    zpe = 0.5 * (H_EV_S * C_CM * f * 1e0).sum() * 1e-3 * 0  # placeholder
    nu = f * C_CM                      # Hz
    e_vib = (0.5 * H_EV_S * nu + H_EV_S * nu / (np.exp(H_EV_S * nu / kT) - 1)).sum()
    # entropy
    x = H_EV_S * nu / kT
    S_vib = 8.617333262e-5 * (x / (np.exp(x) - 1) - np.log(1 - np.exp(-x))).sum()
    return float(e_vib - T * S_vib), float(e_vib), float(S_vib)


if __name__ == "__main__":
    calc = mace_polar(model="polar-1-m", device="cpu", default_dtype="float64",
                      pbc_handling="realspace")
    out = {}
    cases = [("Agem_like", 8.7, 0.8), ("Bgem_like", 6.5, 0.7)]
    for name, lit, err in cases:
        t0 = time.time()
        print(f"\n=== {name} (literature {lit} ± {err} kcal/mol) ===", flush=True)
        prec = read(os.path.join(HERE, f"{name}_min.xyz"))
        prec.calc = calc
        prec.info["charge"], prec.info["spin"] = 0, 1
        syms = prec.get_chemical_symbols()
        # rebuild index bookkeeping identical to scan_barrier.build_complex
        _at, i, j = sb.build_complex(
            os.path.join(HERE, "..", "model", "xtb",
                         "ch3mgcl_2thf_opt.xyz" if name.startswith("A") else "mgme2_2thf_opt.xyz"))
        print(f"  forming bond C{i}-C{j}", flush=True)

        # ---- 1. fine scan (cached, so restarts are cheap)
        cache = os.path.join(HERE, f"{name}_fine.json")
        if os.path.exists(cache):
            d = json.load(open(cache))
            prof, e0 = d["prof"], d["e0"]
            geoms = [read(os.path.join(HERE, f"{name}_fine_{n:02d}.xyz"))
                     for n in range(len(prof))]
            print(f"  loaded cached fine scan ({len(prof)} points)", flush=True)
        else:
            ramp = list(np.round(np.arange(4.40, 2.59, -0.20), 2))
            fine = list(np.round(np.arange(2.60, 2.19, -0.05), 2))
            prof, geoms, e0 = fine_scan(calc, prec, i, j, ramp + fine, tag=name)
            if max(p_["dE_kcal"] for p_ in prof) > 200:
                raise RuntimeError("energy blew up during ramp")
            for n_, g_ in enumerate(geoms):
                write(os.path.join(HERE, f"{name}_fine_{n_:02d}.xyz"), g_)
            json.dump({"prof": prof, "e0": e0}, open(cache, "w"), indent=1)
        near = [n for n, p_ in enumerate(prof) if p_["r"] <= 2.601]
        k = max(near, key=lambda n: prof[n]["dE_kcal"])
        e_scan = prof[k]["dE_kcal"]
        print(f"  fine-scan max: {e_scan:.2f} kcal/mol at r={prof[k]['r']}", flush=True)
        # ---- 2. dimer, seeded by the scan eigenvector
        g0, g1 = geoms[max(k - 1, 0)], geoms[min(k + 1, len(geoms) - 1)]
        dv = (g1.get_positions() - g0.get_positions()).ravel()
        dv /= np.linalg.norm(dv)
        seed = geoms[k].copy()
        seed.calc = calc
        seed.info["charge"], seed.info["spin"] = 0, 1
        e_seed = seed.get_potential_energy()
        try:
            ts, nsteps = dimer_ts(seed, dv.reshape(-1, 3),
                                  dimer_log=os.path.join(HERE, f"{name}_dimer.log"))
            ts.calc = calc
            e_ts = ts.get_potential_energy()
            print(f"  dimer converged in {nsteps} steps, ΔE = {(e_ts - e0) * KCAL:.2f} kcal/mol",
                  flush=True)
            write(os.path.join(HERE, f"{name}_TS.xyz"), ts)
        except Exception as exc:
            print(f"  dimer FAILED ({exc}); falling back to scan maximum", flush=True)
            ts, e_ts = geoms[k].copy(), geoms[k].get_potential_energy()

        # ---- 3. frequencies
        f_ts = freqs(ts, f"{name}_ts")
        f_min = freqs(prec, f"{name}_min")
        real_ts, imag_ts = split_modes(f_ts)
        real_min, imag_min = split_modes(f_min)
        n_imag_ts = int(len(imag_ts))
        print(f"  TS imaginary modes: {n_imag_ts} -> {np.round(imag_ts, 1)}", flush=True)
        print(f"  lowest TS real modes: {np.round(real_ts[:3], 1)}", flush=True)
        print(f"  reactant imaginary modes: {len(imag_min)} -> {np.round(imag_min, 1)}", flush=True)
        if n_imag_ts != 1:
            print(f"  WARNING: 一阶鞍点应恰有 1 个虚频，实测 {n_imag_ts} 个", flush=True)

        # ---- 4. thermochemistry (exclude imaginary modes)
        dG_rrho, Ev, S = rrho_thermo(real_ts)
        dG_rrho_min, _, _ = rrho_thermo(real_min)
        dG_q, _, _ = rrho_thermo(real_ts, scale_low=100.0)
        dG_q_min, _, _ = rrho_thermo(real_min, scale_low=100.0)
        dE = (e_ts - e0) * KCAL
        out[name] = {
            "literature_kcal": lit, "lit_err": err,
            "dE_scan_kcal": round(float(e_scan), 2),
            "dE_ts_kcal": round(float(dE), 2),
            "n_imag_TS": n_imag_ts,
            "imag_cm": [round(float(x), 1) for x in imag_ts],
            "n_imag_reactant": int(len(imag_min)),
            "dG_rrho_kcal": round(float(dE + (dG_rrho - dG_rrho_min) * KCAL), 2),
            "dG_qrrho_kcal": round(float(dE + (dG_q - dG_q_min) * KCAL), 2),
            "lowest_modes_TS": [round(float(x), 1) for x in real_ts[:5]],
            "wall_s": round(time.time() - t0, 1),
        }
        print(f"  ΔE‡={dE:.2f}  ΔG‡(RRHO)={out[name]['dG_rrho_kcal']:.2f}  "
              f"ΔG‡(qRRHO)={out[name]['dG_qrrho_kcal']:.2f} kcal/mol", flush=True)
        json.dump(out, open(os.path.join(HERE, "ts_thermo.json"), "w"), indent=1)
    print("\nDONE")
