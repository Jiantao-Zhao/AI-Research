import time, numpy as np
from ase.io import read
from mace.calculators import mace_polar

at = read("run/box_min.xyz")
at.info["charge"]=0; at.info["spin"]=1
print("atoms", len(at), "cell", round(at.cell[0,0],2), flush=True)
rng=np.random.default_rng(0)

for pbc in ("realspace","pbc"):
    calc = mace_polar(model="polar-1-m", device="cuda", default_dtype="float32", pbc_handling=pbc)
    at.calc = calc
    at.get_forces()
    t0=time.time(); n=10
    for _ in range(n):
        at.positions += rng.normal(0,0.002,(len(at),3))
        at.get_forces()
    print(f"{pbc:10s} {(time.time()-t0)/n*1000:8.1f} ms/step", flush=True)
