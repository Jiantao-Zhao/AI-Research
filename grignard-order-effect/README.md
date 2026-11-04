# Addition-Order Effects in Barbier/Grignard Reactions

A protocol-resolved kinetic framework that treats "order of addition" as the ratio of
reagent-supply rate to substrate-consumption rate, tested against published yield data.

**Short paper:** [`paper/文章.md`](paper/文章.md) · **Supplementary index:** [`paper/补充材料.md`](paper/补充材料.md)

## What this is

Same charge, same stoichiometry, different order of addition — and the yield differs by
25–39 percentage points. This repository contains a kinetic model that reproduces that gap,
identifies it as predominantly supply-limited incomplete conversion rather than side
reaction, and delimits where the framework fails.

## Headline numbers

| Quantity | Value |
|---|---|
| Experimental gap (lit.) | 25–39 percentage points |
| Model gap at A/V = 20 cm⁻¹ | −38.5 points |
| Model gap at A/V = 18.2 cm⁻¹ (fitted) | −43.7 points |
| Share from incomplete conversion | ~84% (model-inferred) |
| Share from MPV | ~16% (model-inferred) |
| One A/V per group | RMS 4.1 (mechanichal, 3 pts), 12.5 (solution, 4 pts) |
| Domain of validity | aldehydes + primary bromides/chlorides |

## Limitations (read before citing)

- Barriers come from constrained scans. The dimer/frequency check initially used `.real` on
  complex frequencies, so imaginary modes registered as zero. **Re-run after the fix:** the
  addition-path stationary point has one significant imaginary mode (356.9 cm⁻¹, the C–C
  forming stretch), i.e. it *is* a first-order saddle along that coordinate. However the
  low-frequency part of the spectrum is unreliable — the reference minimum itself shows three
  imaginary modes (up to 76.9 cm⁻¹) — so ZPE and entropy corrections from these Hessians are
  not usable.
- All barriers are **single conformer, single path**. Two conformers of the same system at
  r = 4.40 Å differ by 3.84 kcal/mol, larger than the literature uncertainty (±0.8).
- The MPV entropy correction is an assumed bimolecular value; the predicted gap is sensitive
  to it (−0.73 at ΔS‡ = −20, −0.385 at −25, −0.33 at −30).
- A/V is partly **inferred from yields** in the condition-dependence sweep, so those
  residuals are not an independent validation. The independent checks are the geometric
  estimate (400 µm → 17, 200 µm → 35 cm⁻¹) and the monotonic ordering of A/V across conditions.
- The 84/16 split is **model output, not measured**. The source paper reports neither
  unreacted-aldehyde recovery nor quantitative yields of the MPV by-products.
- The model's fit does not test the quantum-chemical rate constants: k_G sits in an
  insensitive regime, and k_MPV contributes only ~16% of the gap.

## Layout

```
paper/             short article and supplementary index
model/             reaction network, parameter sweeps, HTML reports
bench/             barrier scans, 2D PES, MPV barrier, interactive dashboards
protocol-lit/      literature dataset extracted from the source paper
lit/               reference PDFs (not redistributed)
notes/             early-stage planning documents
```

Entry points:

- `model/microkinetic.py` — reaction network
- `model/av_model.py` — surface-area sweep and yield inversion
- `bench/scan_barrier.py` — addition-barrier relaxed scans
- `bench/mpv_barrier.py` — MPV hydride-transfer scan
- `bench/scan2d.py` — 2D potential-energy surface

## Reproducing

Requires `mace-torch`, `ase`, `numpy`, `scipy`, `matplotlib`. Calculations used
MACE-POLAR-1-M on CPU (`default_dtype="float64"`, `pbc_handling="realspace"`).

```bash
python model/microkinetic.py          # baseline protocol comparison
python model/av_model.py              # surface-area sweep + report
python bench/scan_barrier.py --outdir out/
```

## Data and licences

Experimental yields from Deitmann et al., *Angew. Chem. Int. Ed.* 2023, 62, e202305775.
Grignard-formation rate constant from *J. Phys. Chem. A* 2001 (DOI 10.1021/jp013357d).
Publisher PDFs are not redistributed.

Code: MIT. Text and figures: CC BY 4.0.
