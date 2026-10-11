# Supplementary Material — Dengue SEIR/SEI Computational Framework

This folder is the complete reproducibility package for the manuscript.
It contains the Python source code, the single-command pipeline that
regenerates all four figures, and the generated figure files. No external
data files are required.

## Contents

| File | Description |
|------|-------------|
| `model.py` | Coupled human-vector SEIR/SEI ODE system, temperature forcing (biting rate + EIP), seasonal temperature driver, RK45 solver |
| `r0.py` | Next-generation matrix construction and spectral-radius R0 computation |
| `sensitivity.py` | Latin hypercube sampling (n = 2,000) with PRCC analysis against R0 |
| `scenarios.py` | Five scenarios: Baseline, Vector control, Warming +1°C, Warming +2°C, Import 4x |
| `run_all.py` | Single-command pipeline that regenerates all four figures (300 dpi PNG previews) |
| `make_figures.py` | Publication-grade figure builder — same pipeline, journal styling (300 dpi TIFF + vector PDF, Okabe–Ito palette, no in-figure titles) |
| `Figure1_r0_temperature.png` | R0 vs temperature curve |
| `Figure2_prcc_sensitivity.png` | PRCC global-sensitivity bar chart |
| `Figure3_scenario_peaks.png` | Scenario peak prevalence comparison |
| `Figure4_incidence_timeseries.png` | Incidence time series |
| `figures_publication/` | The four journal figures as 300 dpi TIFF and vector PDF (as submitted) |

## Requirements

- Python 3.8 or newer
- NumPy
- SciPy
- Matplotlib (Pillow is pulled in automatically as a Matplotlib dependency)

## Reproducibility

All random sampling uses a fixed seed (`20260816`). Results are exactly
reproducible across machines and Python versions.

## Usage

From this `Supplementary/` directory:

```bash
python run_all.py
```

This regenerates the four figures as 300 dpi PNG files in the current
directory.

For the journal-grade figures (300 dpi TIFF + vector PDF, matching the
figures in the manuscript):

```bash
python make_figures.py
```

Output lands in `../02_Figures/` when run from the submission layout, or
alongside the sources with a adjusted `OUT` path; the copies actually used
in the manuscript are in `figures_publication/`.

## Verified numerical outputs

The pipeline asserts the following key outputs (within tolerance):

- R0 at 28°C = 2.93
- Baseline peak infectious = 1,650
- Vector control peak = 120 (-93%)
- Warming +1°C peak = 1,796
- Warming +2°C peak = 1,904
- Import 4x peak = 1,721

## Key parameters

| Parameter | Symbol | Baseline value |
|-----------|--------|----------------|
| Human population | Nh | 10,000 |
| Vector-to-human ratio | Nv/Nh | 10 |
| Human natural mortality | mu_h | 1/(70x365) /day |
| Human incubation rate | sigma_h | 1/5.5 /day |
| Human recovery rate | gamma_h | 1/5 /day |
| Vector-to-human transmission | beta_vh | 0.70 |
| Human-to-vector transmission | beta_hv | 0.60 |
| Adult mosquito mortality | mu_v | 1/12 /day |
| Daily biting rate (28°C) | b | 0.21 /day |
| Extrinsic incubation period (28°C) | 1/sigma_v | 4.25 days |

## Temperature forcing

- Biting rate: b(T) = 0.0043*T + 0.0943, clamped to [0.18, 0.23]
- EIP rate: thermodynamic model from Li et al. (2016) [18], calibrated for
  Aedes albopictus in Guangdong
- Seasonal temperature: T(t) = 22 + 7.8*cos(2*pi*(t-200)/365)

## Data sources

No primary data were collected. All parameters were taken from published
literature (cited in the main manuscript: references [14], [15], [16],
[17], [18]).

## Repository

The complete code is also available at:
https://github.com/GustaveKateb0709/aedes-dengue-seir.git
