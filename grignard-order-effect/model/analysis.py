#!/usr/bin/env python3
"""Sweeps and figures for the protocol-resolved microkinetic model."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import microkinetic as mk

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

BASE = mk.Params(t_add_S=3600.0, k_pin=2.0e-5, k_red=3.0e-4)
PROTOS = ("preformed", "barbier", "inverse")


def run_all(p):
    out = {}
    for pr in PROTOS:
        r = mk.simulate(p, pr)
        out[pr] = {"P": r["yields"]["P"], "Red": r["yields"]["Red"],
                   "Pin": r["yields"]["Pin"], "W": r["yields"]["W"],
                   "conv": r["S_conv"], "t_end": float(r["t"][-1]),
                   "exposure_integral": r["exposure_integral"]}
    out["exposure_s"] = {pr: mk.feed_end(p, pr) for pr in PROTOS}
    out["dP_bar_minus_pre"] = out["barbier"]["P"] - out["preformed"]["P"]
    out["dP_inv_minus_pre"] = out["inverse"]["P"] - out["preformed"]["P"]
    return out


def fig1_timecourses(p=BASE):
    runs = {pr: mk.simulate(p, pr) for pr in PROTOS}
    tmax = max(r["t"][-1] for r in runs.values()) / 60.0
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.7), sharey=True, sharex=True)
    colors = {"RX": "#c0392b", "S": "#2980b9", "G": "#27ae60", "D": "#8e44ad", "P": "#111111"}
    for ax, proto in zip(axes, PROTOS):
        r = runs[proto]
        t = r["t"] / 60.0
        for key, c in colors.items():
            ax.plot(t, r["series"][key], lw=1.9, color=c,
                    label={"G": "RMgX", "D": "R$_2$Mg", "S": "carbonyl"}.get(key, key))
        ax.axvline(mk.feed_end(p, proto) / 60.0, color="grey", ls=":", lw=1)
        ax.set_title(f"{proto}", fontsize=11)
        ax.set_xlabel("time / min")
        ax.grid(alpha=.3)
        ax.text(0.03, 0.06,
                f"reduction {100*r['yields']['Red']:.1f}%   product {100*r['yields']['P']:.1f}%",
                transform=ax.transAxes, fontsize=8.5,
                bbox=dict(fc="white", ec="0.7", lw=.6))
    axes[0].set_ylabel("concentration / M")
    axes[0].set_xlim(0, tmax)
    axes[0].legend(fontsize=8, frameon=False, loc="upper left")
    fig.suptitle("Same charge, same dosing rate; the difference is the free-substrate / free-site overlap",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_timecourses.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig2_map():
    kreds = np.logspace(-4.5, -2.5, 17)
    ks = np.logspace(-0.5, 2, 17)
    Z = np.full((len(ks), len(kreds)), np.nan)
    rows = []
    for i, K in enumerate(ks):
        for j, kr in enumerate(kreds):
            p = BASE.copy(K_S=float(K), k_red=float(kr))
            d = run_all(p)
            Z[i, j] = d["dP_bar_minus_pre"]
            rows.append({"K_S": float(K), "k_red": float(kr), "dP": d["dP_bar_minus_pre"],
                         "P_pre": d["preformed"]["P"], "P_bar": d["barbier"]["P"],
                         "P_inv": d["inverse"]["P"]})
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    lim = 0.25
    im = ax.imshow(Z, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=0.0,
                   extent=[np.log10(kreds[0]), np.log10(kreds[-1]),
                           np.log10(ks[0]), np.log10(ks[-1])])
    ax.set_xlabel(r"$\log_{10}(k_{\rm red})$   [substrate reduction on Mg]")
    ax.set_ylabel(r"$\log_{10}(K_S)$   [substrate surface affinity, L/mol]")
    cb = fig.colorbar(im, ax=ax)
    cs = ax.contour(np.log10(kreds), np.log10(ks), Z, levels=[-0.20, -0.10, -0.05, -0.02],
                    colors="k", linewidths=.7)
    ax.clabel(cs, fmt="%.2f", fontsize=7)
    cb.set_label(r"$\Delta$ product yield (Barbier $-$ preformed)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_phase_map.png"), dpi=200)
    plt.close(fig)
    return rows, Z


def fig3_curves():
    kreds = np.logspace(-4.5, -2.5, 21)
    curves = {pr: [] for pr in PROTOS}
    rows = []
    for kr in kreds:
        p = BASE.copy(k_red=float(kr))
        d = run_all(p)
        for pr in PROTOS:
            curves[pr].append(d[pr]["P"])
        rows.append({"k_red": float(kr), **{f"P_{pr}": d[pr]["P"] for pr in PROTOS},
                     "Red_preformed": d["preformed"]["Red"],
                     "Red_barbier": d["barbier"]["Red"]})
    fig, ax = plt.subplots(figsize=(5.8, 3.9))
    for pr, c in zip(PROTOS, ("#111111", "#c0392b", "#2980b9")):
        ax.plot(kreds, curves[pr], "o-", ms=3, lw=1.6, color=c, label=pr)
    ax.set_xscale("log")
    ax.set_xlabel(r"$k_{\rm red}$ / M$^{-1}$s$^{-1}$   (how easily the metal reduces the substrate)")
    ax.set_ylabel("product yield")
    ax.grid(alpha=.3)
    ax.legend(frameon=False, fontsize=8, title="addition order")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig3_protocol_curves.png"), dpi=200)
    plt.close(fig)
    return rows


def movchan_check():
    """Normalise the effect by the difference in metal-exposure time."""
    p = BASE
    d = run_all(p)
    exp = d["exposure_s"]
    out = {}
    for pr in PROTOS:
        out[pr] = {"P": d[pr]["P"], "exposure_min": exp[pr] / 60.0,
                   "P_per_exposure": d[pr]["P"] / (exp[pr] / 60.0)}
    return out


if __name__ == "__main__":
    base = run_all(BASE)
    print(json.dumps(base, indent=2))
    fig1_timecourses()
    rows, Z = fig2_map()
    r3 = fig3_curves()
    mc = movchan_check()
    json.dump({"baseline": base, "map": rows, "curves": r3, "exposure_check": mc},
              open(os.path.join(HERE, "sweep_results.json"), "w"), indent=1)
    print("map dP range: %.4f .. %.4f" % (np.nanmin(Z), np.nanmax(Z)))
    print("curves dP at kred=1e-3: %.4f" % (r3[-1]["P_barbier"] - r3[-1]["P_preformed"]))
