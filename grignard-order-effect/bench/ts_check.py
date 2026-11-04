"""用修正后的频率分析复核已保存的"TS"结构：真鞍点还是极小值？"""
import os, json
import numpy as np
import torch
torch.serialization.add_safe_globals([slice])
from ase.io import read
from ase.vibrations import Vibrations
from mace.calculators import mace_polar
from ts_thermo import split_modes

HERE = os.path.dirname(os.path.abspath(__file__))
calc = mace_polar(model="polar-1-m", device="cpu", default_dtype="float64",
                  pbc_handling="realspace")

out = {}
for name in ("Agem_like_TS", "Bgem_like_TS", "Agem_Bz_min"):
    p = os.path.join(HERE, f"{name}.xyz")
    if not os.path.exists(p):
        continue
    at = read(p)
    at.info["charge"] = 0; at.info["spin"] = 1
    at.calc = calc
    at.pbc = False; at.set_cell([0, 0, 0])
    vib = Vibrations(at, name=f"chk_{name}", delta=0.01)
    vib.run()
    f = np.asarray(vib.get_frequencies(), dtype=complex)
    vib.clean()
    real, imag = split_modes(f)
    out[name] = {"n_imag": int(len(imag)),
                 "imag_cm": [round(float(x), 1) for x in imag],
                 "lowest_real": [round(float(x), 1) for x in real[:5]]}
    print(f"{name}: 虚频 {len(imag)} 个 -> {np.round(imag,1)}")
    print(f"   最低实模 {np.round(real[:5],1)}", flush=True)

json.dump(out, open(os.path.join(HERE, "ts_check.json"), "w"), indent=1)
