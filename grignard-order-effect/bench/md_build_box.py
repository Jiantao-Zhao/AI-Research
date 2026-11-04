"""建显式 THF 溶剂盒子：Monte Carlo 插入，18 Å 立方。"""
import numpy as np
from ase import Atoms
from ase.io import read, write


def random_rotation(rng):
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-z*w),   2*(x*z+y*w)],
        [2*(x*y+z*w),   1-2*(x*x+z*z), 2*(y*z-x*w)],
        [2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x*x+y*y)]])


def build(n_thf=26, box=18.0, seed=7, out="box.xyz"):
    rng = np.random.default_rng(seed)
    mol = read("agm.xyz")
    mol.info["charge"] = 0; mol.info["spin"] = 1
    idx = list(range(6, 19))
    thf = mol[idx]
    tp0 = thf.get_positions() - thf.get_positions().mean(axis=0)

    molp = mol.get_positions() - mol.get_positions().mean(axis=0) + box/2
    pos = list(molp); sym = list(mol.get_chemical_symbols())

    def mindist(cand):
        P = np.array(pos)
        best = 9.0
        for p_ in cand:
            d = P - p_
            d -= box*np.round(d/box)
            best = min(best, np.sqrt((d*d).sum(1)).min())
        return best

    placed = 0
    for _ in range(6000):
        if placed >= n_thf:
            break
        c = rng.uniform(2.0, box-2.0, 3)
        cand = tp0 @ random_rotation(rng).T + c
        if mindist(cand) > 2.05:
            pos += list(cand); sym += thf.get_chemical_symbols(); placed += 1

    at = Atoms(symbols=sym, positions=pos)
    at.set_cell([box, box, box]); at.set_pbc(True)
    at.info["charge"] = 0; at.info["spin"] = 1
    write(out, at)
    print("放入 THF:", placed, " 总原子:", len(at))
    d = at.get_all_distances(); np.fill_diagonal(d, 9)
    print("最近距离 %.2f Å" % d.min())
    return at


if __name__ == "__main__":
    build()
