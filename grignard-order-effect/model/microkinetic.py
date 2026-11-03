#!/usr/bin/env python3
"""
Protocol-resolved microkinetic model for the Barbier / Grignard addition.

Question the model answers: at fixed overall stoichiometry, under what conditions does
the *feed trajectory* (order of addition) change the product distribution?

Mechanism (elementary mass action; homogeneous steps in M, s units):

  F   RX + Mg(site)     -> RMgX                       k_form   reagent synthesis
  R   S  + Mg(site)     -> Red                        k_red    substrate reduction
  W   RMgX + RX         -> W                          k_wurtz  Wurtz coupling
  A1  RMgX + S          -> P                          k_G      addition via monomer
  A2  R2Mg + 2 S        -> 2 P                        k_D      addition via Schlenk product
  S1  2 RMgX <-> R2Mg + MgX2                          k_sf, k_sr
  Pn  2 S + Mg(site)    -> Pin                        k_pin    pinacol-type coupling
  MPV P  + S            -> Red                        k_mpv    Meerwein-Ponndorf-Verley:
                                                                 alkoxide transfers H to a
                                                                 still-unreacted aldehyde;
                                                                 the product is oxidised back
                                                                 to S, the aldehyde is reduced

Surface: competitive Langmuir adsorption of RX and S on Mg. The free-site fraction is
theta0 = 1 / (1 + K_RX*[RX] + K_S*[S]). The substrate therefore poisons the surface that
makes the reagent -- this is the only link from feed order to chemistry, and it is enough
to generate order dependence.

The model does not claim literature-accurate rate constants. All constants are stated
assumptions whose *ratios* are the physically meaningful inputs; the results are reported
as functions of those ratios.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict

import numpy as np
from scipy.integrate import solve_ivp

KB = 1.380649e-23
HPLANK = 6.62607015e-34
R_GAS = 8.314462618e-3
T_REF = 298.15


def eyring(dg_kj: float, T: float = T_REF) -> float:
    """Eyring rate constant (s^-1) from a barrier in kJ/mol."""
    return (KB * T / HPLANK) * np.exp(-dg_kj / (R_GAS * T))


@dataclass
class Params:
    # --- charge (mol/L) ---
    RX_total: float = 1.00
    S_total: float = 1.00
    sites: float = 1.00           # accessible Mg surface (mol/L), held constant (Mg in excess)
    # --- schedule (s) ---
    t_feed: float = 3600.0        # duration of RX addition
    t_age: float = 900.0          # preformed: hold after RX addition, before S
    t_add_S: float = 120.0        # preformed: duration of S addition
    t_hold: float = 1800.0        # fixed post-addition stir before workup
    # --- kinetics ---
    k_form: float = 8.0e-3        # L mol^-1 s^-1 (per mol/L sites)
    k_red: float = 1.0e-4         # L mol^-1 s^-1
    k_wurtz: float = 4.0e-5       # L mol^-1 s^-1
    k_G: float = 3.0e-3           # L mol^-1 s^-1
    k_D: float = 3.0e-2           # L mol^-1 s^-1  (Schlenk product reacts faster)
    k_sf: float = 2.0e0           # L mol^-1 s^-1  2 RMgX -> R2Mg + MgX2
    k_sr: float = 2.0e0           # L mol^-1 s^-1  R2Mg + MgX2 -> 2 RMgX
    k_pin: float = 2.0e-4         # L^2 mol^-2 s^-1
    k_mpv: float = 0.0            # L mol^-1 s^-1  (0 = channel off)
    # 若给出 k_app（s^-1），生成速率用它替代 k_form*sites。
    # 物理含义：k_app = k_s * (A/V)，k_s 为文献绝对速率常数（~1e-4 cm/s），
    # A/V 为每单位体积的镁表面积（cm^2/cm^3 = cm^-1）。
    k_app: float = None
    # --- adsorption ---
    K_RX: float = 0.5             # L/mol
    K_S: float = 12.0             # L/mol   (carbonyl binds Mg far more strongly)

    def copy(self, **kw):
        d = asdict(self)
        d.update(kw)
        return Params(**d)


def theta_free(rx: float, s: float, p: Params) -> float:
    return 1.0 / (1.0 + p.K_RX * max(rx, 0.0) + p.K_S * max(s, 0.0))


def feed(t: float, p: Params, protocol: str):
    if protocol == "preformed":
        r_rx = p.RX_total / p.t_feed if t < p.t_feed else 0.0
        t0 = p.t_feed + p.t_age
        r_s = p.S_total / p.t_add_S if t0 <= t < t0 + p.t_add_S else 0.0
        return r_rx, r_s
    if protocol == "barbier":
        on = t < p.t_feed
        return (p.RX_total / p.t_feed if on else 0.0,
                p.S_total / p.t_feed if on else 0.0)
    if protocol == "inverse":
        r_rx = p.RX_total / p.t_feed if t < p.t_feed else 0.0
        return r_rx, 0.0
    if protocol == "barbier_batch":
        return p.t_feed, 0.0
    if protocol == "preformed_batch":
        t0 = p.t_feed + p.t_age
        r_s = p.S_total / p.t_add_S if t0 <= t < t0 + p.t_add_S else 0.0
        return 0.0, r_s
    if protocol == "preformed_bolus":
        t0 = p.t_feed + p.t_age
        r_s = p.S_total / p.t_add_S if t0 <= t < t0 + p.t_add_S else 0.0
        return 0.0, r_s
    raise ValueError(protocol)


def rhs(t, y, p: Params, protocol: str):
    RX, S, G, D, M, P, W, Red, Pin, _I, _J = y
    RX, S = max(RX, 0.0), max(S, 0.0)
    th0 = theta_free(RX, S, p)

    kf = p.k_app if p.k_app is not None else p.k_form * p.sites
    r_form = kf * RX * th0
    r_red = p.k_red * S * th0 * p.sites
    r_wurtz = p.k_wurtz * RX * G
    r_G = p.k_G * G * S
    r_D = p.k_D * D * S
    r_sf = p.k_sf * G * G
    r_sr = p.k_sr * D * M
    r_pin = p.k_pin * S * S * th0 * p.sites
    r_mpv = p.k_mpv * P * S
    f_rx, f_s = feed(t, p, protocol)

    return [
        f_rx - r_form - r_wurtz,
        f_s - r_red - r_G - 2.0 * r_D - 2.0 * r_pin,
        r_form - 2.0 * r_sf + 2.0 * r_sr - r_G - r_wurtz,
        r_sf - r_sr - r_D,
        r_sf - r_sr,
        r_G + 2.0 * r_D - r_mpv,
        r_wurtz,
        r_red + r_mpv,
        r_pin,
        S * th0 * p.sites,          # substrate-on-metal exposure integral
        P * S,                      # alkoxide x free aldehyde integral (MPV driving)
    ]


NAMES = ["RX", "S", "G", "D", "M", "P", "W", "Red", "Pin", "I_exp", "I_mpv"]


def feed_end(p: Params, protocol: str) -> float:
    """Time at which the last addition ends; workup follows t_hold later."""
    if protocol == "preformed":
        return p.t_feed + p.t_age + p.t_add_S
    if protocol == "barbier_batch":
        return p.t_feed
    if protocol in ("preformed_batch", "preformed_bolus"):
        return p.t_feed + p.t_age + p.t_add_S
    if protocol in ("barbier", "inverse"):
        return p.t_feed
    raise ValueError(protocol)


def simulate(p: Params, protocol: str, n_out: int = 600):
    rx0 = p.RX_total if protocol in ("barbier_batch", "preformed_batch", "preformed_bolus") else 0.0
    s0 = p.S_total if protocol in ("inverse", "barbier_batch") else 0.0
    y0 = [rx0, s0] + [0.0] * 9
    t_end = feed_end(p, protocol) + p.t_hold
    t_eval = np.linspace(0.0, t_end, n_out)
    sol = solve_ivp(rhs, (0.0, t_end), y0, args=(p, protocol), method="LSODA",
                    t_eval=t_eval, rtol=1e-9, atol=1e-12, max_step=10.0)
    series = {n: sol.y[i] for i, n in enumerate(NAMES)}
    fin = {n: float(v[-1]) for n, v in series.items()}
    carbon_out = 2 * fin["W"] + fin["G"] + 2 * fin["D"] + fin["P"]
    sub_out = fin["P"] + fin["Red"] + 2 * fin["Pin"]
    return {
        "protocol": protocol,
        "ok": bool(sol.success),
        "t": sol.t,
        "series": series,
        "final": fin,
        "closure": {"carbon_in": p.RX_total, "carbon_found": carbon_out,
                    "S_in": p.S_total, "S_found": sub_out},
        "yields": {"P": fin["P"] / p.S_total, "Red": fin["Red"] / p.S_total,
                   "Pin": 2 * fin["Pin"] / p.S_total, "W": fin["W"] / p.RX_total},
        "S_conv": 1.0 - fin["S"] / p.S_total,
        "exposure_integral": fin["I_exp"],
        "mpv_integral": fin["I_mpv"],
    }


def schlenk_tau(p: Params) -> float:
    """Relaxation time of 2G <-> D + M at the equilibrium composition of the reagent alone."""
    # solve the isolated Schlenk equilibrium for G0 = 0.5 M total Mg
    from scipy.optimize import brentq
    Gtot = 0.5

    def resid(g):
        # 2g <-> d + m with d = m = (Gtot-g)/2 , K = ksf/ksr
        d = (Gtot - g) / 2.0
        return p.k_sf * g * g - p.k_sr * d * d

    g = brentq(resid, 1e-8, Gtot - 1e-8)
    d = (Gtot - g) / 2.0
    J = np.array([[-4 * p.k_sf * g, 2 * p.k_sr * d],
                  [2 * p.k_sf * g, -p.k_sr * d]])
    lam = max(abs(np.linalg.eigvals(J)))
    return 1.0 / lam if lam > 0 else np.inf


def protocol_contrast(p: Params) -> dict:
    pre = simulate(p, "preformed")
    bar = simulate(p, "barbier")
    return {
        "sel_P_preformed": pre["yields"]["P"],
        "sel_P_barbier": bar["yields"]["P"],
        "delta_P": bar["yields"]["P"] - pre["yields"]["P"],
        "delta_W": bar["yields"]["W"] - pre["yields"]["W"],
        "delta_Red": bar["yields"]["Red"] - pre["yields"]["Red"],
        "tau_schlenk": schlenk_tau(p),
        "t_feed_over_tau": p.t_feed / schlenk_tau(p),
    }


if __name__ == "__main__":
    p = Params()
    out = {pr: simulate(p, pr)["yields"] for pr in ("preformed", "barbier", "inverse")}
    print(json.dumps(out, indent=2))
    print(json.dumps(protocol_contrast(p), indent=2))
