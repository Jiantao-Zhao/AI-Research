#!/usr/bin/env python3
"""Level-0 speciation engine test: GFN2-xTB opt + RRHO free energies in implicit THF.

This is a feasibility probe, not production thermochemistry: implicit solvation only,
no explicit THF, no configurational averaging. Its job is to (i) prove the pipeline and
(ii) quantify how sensitive the Schlenk free energy is to the solvation treatment.
"""
import json, os, re, subprocess, sys

XTB = "/opt/anaconda3/envs/qc/bin/xtb"
HERE = os.path.dirname(os.path.abspath(__file__))
SPECIES = ["ch3mgcl_monomer", "mgch32_monomer", "mgcl2_monomer", "ch3mgcl_dimer"]
SOLV = sys.argv[1] if len(sys.argv) > 1 else "thf"


def run(args, tag):
    p = subprocess.run(args, cwd=HERE, capture_output=True, text=True)
    log = os.path.join(HERE, f"{tag}.log")
    with open(log, "w") as f:
        f.write(p.stdout + "\n" + p.stderr)
    return p.stdout + p.stderr


def parse_g(text):
    m = re.search(r"TOTAL FREE ENERGY\s+(-?\d+\.\d+)", text)
    e = re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)", text)
    imag = re.search(r"# imaginary freq\.\s+(\d+)", text)
    return (
        float(m.group(1)) if m else None,
        float(e.group(1)) if e else None,
        int(imag.group(1)) if imag else 0,
    )


def optimize(src, tag):
    out = run(
        [XTB, src, "--opt", "tight", "--alpb", SOLV, "--chrg", "0", "--uhf", "0"],
        tag + "_opt",
    )
    return out


def freq(xyz, tag):
    out = run([XTB, xyz, "--ohess", "--alpb", SOLV, "--chrg", "0", "--uhf", "0"], tag + "_freq")
    return parse_g(out), out


results = {}
for sp in SPECIES:
    src = sp + ".xyz"
    optimize(src, sp)
    opt_xyz = os.path.join(HERE, "xtbopt.xyz")
    if not os.path.exists(opt_xyz):
        results[sp] = {"error": "no xtbopt"}
        continue
    target = os.path.join(HERE, sp + "_opt.xyz")
    os.replace(opt_xyz, target)
    (g, e, imag), out = freq(target, sp)
    results[sp] = {"G": g, "E": e, "n_imag": imag}

with open(os.path.join(HERE, f"schlenk_{SOLV}.json"), "w") as f:
    json.dump(results, f, indent=2)
print(json.dumps(results, indent=2))
