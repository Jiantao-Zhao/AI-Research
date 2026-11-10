"""MPV 的活化熵：反应物 vs 过渡态的振动频率。"""
import os, json
import numpy as np
import torch
torch.serialization.add_safe_globals([slice])
from ase.io import read
from ase.constraints import FixBondLength
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from mace.calculators import mace_polar
from ts_thermo import rrho_thermo, split_modes

HERE = os.path.dirname(os.path.abspath(__file__))
KCAL = 23.0605
H, CCAR = 4, 8


def freqs(at, tag, calc):
    w = at.copy(); w.calc = calc
    w.pbc = False; w.set_cell([0, 0, 0])
    v = Vibrations(w, name=tag, delta=0.01)
    v.run()
    f = np.asarray(v.get_frequencies(), dtype=complex)
    v.clean()
    return f


def main():
    calc = mace_polar(model="polar-1-m", device="cpu", default_dtype="float64",
                      pbc_handling="realspace")
    reac = read(os.path.join(HERE, "mpv_reactant.xyz"))
    reac.info["charge"] = 0; reac.info["spin"] = 1
    reac.calc = calc
    e0 = float(reac.get_potential_energy())

    # 扫 H···C 找极大
    best = None
    for r in (1.20, 1.25, 1.30, 1.35, 1.45, 1.60):
        w = reac.copy(); w.calc = calc
        d = w.get_distance(H, CCAR)
        u = (w.positions[CCAR] - w.positions[H]) / d
        w.positions[H] += u * (r - d)
        w.set_constraint(FixBondLength(H, CCAR))
        BFGS(w, logfile=None, maxstep=0.05).run(fmax=0.15, steps=200)
        de = (float(w.get_potential_energy()) - e0) * KCAL
        print(f"  d(H..C)={r:.2f}  ΔE={de:5.2f}", flush=True)
        if best is None or de > best[0]:
            best = (de, r, w)
    dE, rmax, ts = best
    print(f"\n势垒 ΔE‡ = {dE:.2f} kcal/mol at d={rmax:.2f}", flush=True)

    f_r = freqs(reac, "mpv_r", calc)
    f_t = freqs(ts, "mpv_t", calc)
    real_r, imag_r = split_modes(f_r)
    real_t, imag_t = split_modes(f_t)
    print(f"反应物虚频 {len(imag_r)}，过渡态虚频 {len(imag_t)} -> {np.round(imag_t,1)}", flush=True)

    dG_r, Ev_r, S_r = rrho_thermo(real_r)
    dG_t, Ev_t, S_t = rrho_thermo(real_t)
    dS = (S_t - S_r) * 1000 / 1000  # eV/K
    dS_cal = (S_t - S_r) * 23060.5   # eV/K -> cal/mol/K
    out = {"barrier_kcal": float(dE), "r": float(rmax),
           "n_imag_reactant": int(len(imag_r)), "n_imag_ts": int(len(imag_t)),
           "imag_ts_cm": [round(float(x), 1) for x in imag_t],
           "S_reactant_cal": float(S_r * 23060.5), "S_ts_cal": float(S_t * 23060.5),
           "dS_cal": float(dS_cal),
           "dG_rrho_kcal": float(dE + (dG_t - dG_r) * KCAL)}
    print(f"S(反应物) = {S_r*23060.5:.1f} cal/mol/K")
    print(f"S(过渡态) = {S_t*23060.5:.1f} cal/mol/K")
    print(f"ΔS‡ = {dS_cal:.1f} cal/mol/K   （模型原来的假设是 -25）")
    print(f"ΔG‡ = {out['dG_rrho_kcal']:.1f} kcal/mol")
    json.dump(out, open(os.path.join(HERE, "mpv_entropy.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
