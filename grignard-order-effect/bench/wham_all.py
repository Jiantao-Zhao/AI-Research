"""收集所有窗口的直方图，做一维 WHAM。"""
import glob, json, os, sys
import numpy as np
from ase import units
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["Hiragino Sans GB", "PingFang HK", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False


def run(workdir="run", T=300.0, nbins=120, tol=1e-9, maxiter=20000):
    """标准一维 WHAM。

    p(b) = sum_i n_ib / sum_i N_i exp(f_i - V_i(b)) / kT
    f_i  = -kT ln sum_b p(b) exp(-V_i(b)/kT)
    两者交替迭代。
    """
    umb = os.path.join(workdir, "umb")
    metas = [json.load(open(m)) for m in glob.glob(os.path.join(umb, "w_*.json"))]
    if not metas:                       # 兼容旧格式：单个 meta.json
        meta = json.load(open(os.path.join(umb, "meta.json")))
        kk = meta.get("k", 2.0)
        metas = [{"w": w["w"], "r0": w["r0"], "k": kk} for w in meta["windows"]]
    metas.sort(key=lambda m: m["r0"])
    r0s = np.array([m["r0"] for m in metas])
    ks = np.array([m["k"] for m in metas])
    data = [np.load(os.path.join(umb, f"d_{m['w']:02d}.npy")) for m in metas]
    allx = np.concatenate(data)
    edges = np.linspace(allx.min(), allx.max(), nbins + 1)
    centers = 0.5 * (edges[1:] + edges[:-1])
    hist = np.array([np.histogram(x, bins=edges)[0] for x in data], float)
    N = hist.sum(axis=1)
    nb = hist.sum(axis=0)
    kT = T * units.kB
    bias = 0.5 * ks[:, None] * (centers[None, :] - r0s[:, None]) ** 2   # (nwin, nbins)
    expb = np.exp(-bias / kT)
    f = np.zeros(len(r0s))
    p = np.zeros_like(centers)
    for it in range(maxiter):
        denom = (N[:, None] * np.exp(f[:, None] / kT) * expb).sum(axis=0)
        p = np.where(denom > 0, nb / np.maximum(denom, 1e-300), 0.0)
        newf = np.array([-kT * np.log(max((p * expb[i]).sum(), 1e-300))
                         for i in range(len(r0s))])
        newf -= newf[0]
        if np.max(np.abs(newf - f)) < tol:
            f = newf
            break
        f = newf
    denom = (N[:, None] * np.exp(f[:, None] / kT) * expb).sum(axis=0)
    p = np.where(denom > 0, nb / np.maximum(denom, 1e-300), 0.0)
    pmf = -kT * np.log(np.maximum(p, 1e-300)) * 23.0605
    m = nb > 0
    return centers[m], pmf[m], nb, len(r0s)


if __name__ == "__main__":
    wd = sys.argv[1] if len(sys.argv) > 1 else "run"
    r, p, nb, nw = run(wd)
    far = r >= max(r.min(), 3.3)
    ref = p[far].min() if far.sum() >= 2 else p.min()
    p = p - ref
    i = int(np.argmax(np.where(r < 3.1, p, -1e9)))
    print(f"窗口数 {nw}，有效格点 {len(r)}")
    print(f"势垒 {p[i]:.1f} kcal/mol at r={r[i]:.2f} Å")
    print("文献 8.7 ± 0.8")
    json.dump({"r": r.tolist(), "pmf": p.tolist(),
               "barrier": float(p[i]), "r_at": float(r[i]), "nwin": nw},
              open(os.path.join(wd, "pmf.json"), "w"), indent=1)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(r, p, "o-", color="#c0392b", ms=3.5, lw=1.4)
    ax.axhline(8.7, color="#2c6fb5", ls="--", lw=1.2)
    ax.set_xlabel("C1–C20 距离 / Å"); ax.set_ylabel("PMF / kcal mol⁻¹")
    ax.set_title("显式 THF 中的加成自由能面", fontsize=10.5)
    ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig("pmf_solvent.png", dpi=170, facecolor="white")
    print("图已保存 pmf_solvent.png")
