# Supplementary Code — Temperature-Driven Dengue SEIR/SEI Framework

Reproducibility package for the manuscript on temperature-forced dengue
transmission (coupled human–vector SEIR/SEI model, next-generation-matrix R0,
Latin-hypercube PRCC global sensitivity analysis, and five intervention /
climate scenarios).

The package is self-contained: no data files, no network access, and no
external style assets are required. Everything below can be reproduced from a
clean copy of this archive.

---

## 1. Contents

| Path | Description |
|------|-------------|
| `src/model.py` | Coupled SEIR/SEI ODE system, temperature-dependent biting rate and EIP, seasonal temperature forcing, imported-case pulse, RK45 integration |
| `src/r0.py` | Next-generation matrix and spectral-radius R0 |
| `src/sensitivity.py` | Latin hypercube sampling (n = 2,000) + PRCC against R0 |
| `src/scenarios.py` | Five scenarios: Baseline, Vector control, Warming +1 °C, Warming +2 °C, Import 4x |
| `src/config.py` | Central configuration (random seed, environment override) |
| `src/style_presets.py` | Self-contained Matplotlib publication preset (Okabe–Ito palette) |
| `src/figure_export.py` | Self-contained 300 dpi TIF/PDF exporter with JSON provenance manifests |
| `src/make_figures.py` | Builds Figures 1–5 (300 dpi LZW-RGB TIFF + vector PDF) |
| `src/run_all.py` | End-to-end pipeline: all analyses + five 300 dpi PNG figures |
| `tests/test_verification.py` | 21 numerical checks pinning the values reported in the manuscript |
| `requirements.txt` | Exact package versions used for the published results |

---

## 2. Environment

* Python 3.13.12 (CPython) — the version used for the manuscript run.
  Python ≥ 3.9 is expected to work; the test suite will confirm.
* macOS arm64 was used; the code is platform independent.

Pinned dependencies (`requirements.txt`):

```
numpy==2.4.6
scipy==1.18.0
matplotlib==3.11.1
pillow==12.3.0
```

---

## 3. Installation

```bash
# optional but recommended
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

python -m pip install -r requirements.txt
```

No further installation step is needed: the modules in `src/` are imported
directly (no `pip install .`, no `PYTHONPATH` setup — each entry-point script
adds its own directory to `sys.path`).

---

## 4. One-command run

```bash
python src/run_all.py
```

This runs the complete pipeline and prints, in order:

1. R0 at 20/25/28/30/32 °C,
2. the LHS-PRCC sensitivity analysis (n = 2,000),
3. the five scenario simulations,
4. five 300 dpi PNG figures written to the **current working directory**:
   `Figure1_r0_temperature.png`, `Figure2_prcc_sensitivity.png`,
   `Figure3_scenario_peaks.png`, `Figure4_incidence_timeseries.png`,
   `Figure5_prcc_designs.png`.

Runtime is a few seconds on a laptop.

---

## 5. Verification

```bash
python tests/test_verification.py
```

The suite runs without pytest (it never imports it) and prints one line per
check plus a summary. It exits with status 0 only if all 21 checks pass.
`pytest tests/test_verification.py` works too, if pytest is installed.

Expected tail of the output:

```
21 passed, 0 failed of 21 tests
RESULT: ALL TESTS PASSED
```

### Verified values

| Quantity | Expected | Tolerance |
|----------|----------|-----------|
| R0 at 28 °C | 2.928 | ±0.01 |
| R0 at 25 °C | 2.270 | ±0.01 |
| R0 at 30 °C | 3.262 | ±0.01 |
| R0 at 32 °C | 3.492 | ±0.01 |
| R0 at 19 °C | 0.940 (< 1) | ±0.01 |
| EIP at 28 °C, 1/σ_v | 4.251 days | ±1e-3 |
| σ_v at 28 °C | 0.2352 day⁻¹ | ±1e-3 |
| Baseline peak infectious | 1,650 | ±15 |
| Vector control peak | 120 | ±10 |
| Warming +1 °C peak | 1,796 | ±15 |
| Warming +2 °C peak | 1,904 | ±15 |
| Import 4x peak | 1,721 | ±15 |
| Baseline peak day | 233.5 | ±2 days |
| b(15 °C) | 0.180 (clamped) | exact |

The individual analysis scripts also self-check when run directly:

```bash
cd src
python r0.py          # prints R0(15..35 °C) and asserts the key values
python scenarios.py   # prints the five-scenario table and asserts the peaks
python sensitivity.py # prints the PRCC table
```

---

## 6. Reproducing the publication figures

The manuscript figures are the 300 dpi LZW-RGB TIFF / PDF pair produced by:

```bash
python src/make_figures.py
```

Output directory resolution (highest priority first):

1. `--outdir <dir>` on the command line,
2. the `DENGUE_FIG_DIR` environment variable,
3. automatic detection: inside the full manuscript tree
   (`.../04_Code_and_Data/src`) the figures go to `.../02_Figures`; in an
   extracted archive they go to `<archive root>/figures`.

Each figure is written as `<name>.tiff`, `<name>.pdf` and a
`<name>.export.json` manifest recording settings, library versions and the
data provenance (the numeric values plotted).

Reproduction status verified for this release (Matplotlib 3.11.1, DejaVu/Arial
font stack as installed):

* **TIFF files are byte-identical** (same MD5) to the files supplied with the
  manuscript.
* **PDF files differ only in the embedded `/CreationDate` timestamp** (6 bytes);
  all other bytes, including the Type-42 TrueType font subsets, are identical.
  Re-running the script therefore never changes figure content, only the
  creation time.
* If your system lacks Arial, Matplotlib falls back to DejaVu Sans (declared in
  the preset), which can shift glyph shapes slightly; pixel dimensions, page
  size and 300 dpi resolution are unaffected.

`make_figures.py` overwrites existing outputs, so it can be re-run in place.

---

## 7. Changing the random seed

The default seed `20260816` reproduces every number in the manuscript. It can
be overridden without editing any file:

```bash
python src/run_all.py --seed 42          # command-line flag
DENGUE_SEED=42 python src/run_all.py     # environment variable
DENGUE_SEED=42 python -c "import sys; sys.path.insert(0,'src'); import model; print(model.SEED)"
```

Resolution order is `--seed` → `DENGUE_SEED` → default. The seed affects only
the Latin-hypercube draws used by the PRCC analysis; R0, the ODE solutions and
the scenario peaks are deterministic and unchanged (the test suite enforces
this).

---

## 8. Model summary (for orientation)

State vector `[S_h, E_h, I_h, R_h, S_v, E_v, I_v]`, solved with
`scipy.integrate.solve_ivp` (RK45, `max_step = 1 day`, `rtol = 1e-8`,
`atol = 1e-10`) over 1,095 days.

* Biting rate: `b(T) = 0.0043 T + 0.0943`, clamped to `[0.18, 0.23]`
* EIP rate: thermodynamic model of Li et al. (2016) for *Aedes albopictus*
* Seasonal temperature: `T(t) = 22 + 7.8 cos(2π (t − 200)/365)`
* Imported cases: Gaussian pulse peaking at day 210 of each year
* Baseline parameters: `N_h = 10,000`, `N_v/N_h = 10`, `μ_h = 1/(70·365)`,
  `σ_h = 1/5.5`, `γ_h = 1/5`, `β_vh = 0.70`, `β_hv = 0.60`, `μ_v = 1/12`
  (all rates in day⁻¹)

---

## 9. Licence and data availability

No primary data were collected; all parameters come from the literature cited
in the manuscript. Code is released for review and reproduction of the
published results.
