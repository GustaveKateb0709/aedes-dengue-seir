# Dynamic-output sensitivity analyses — summary

All runs use the unmodified model in `../src`; the only intervention is a
run-time replacement of `model.biting_rate` / `model.eip_rate`, never an edit
of the source file. Baseline scenario = seasonal temperature + annual import
pulse, `t_span=(0, 1095)`, 1096 daily output points, `Nh = 10,000`.

## Headline numbers

| Task | Question | Key number |
|---|---|---|
| A | PRCC of the *dynamic* peak vs the static PRCC | baseline peak Ih = 1650.8 at day 234; nv_ratio PRCC 0.950 (static 0.954), b 0.553 (static 0.573) |
| B | How does the 1,650 headline scale with Nv/Nh? | peak Ih spans 5 (Nv/Nh=1) to 2258 (Nv/Nh=20) — a 422× span; the 3-yr attack rate saturates at 96.9% for Nv/Nh ≥ 10 |
| C | What if the b(T) clamp is removed? | peak Ih changes by -2.2e-07 % (i.e. 3.7e-6 persons); the only measurable effect is on the threshold, 19.39 → 19.46 °C |
| D | Is the ranking an artefact of range width? | **Partly.** Signs are all preserved, but the exact order is not: b moves 5→1 and nv_ratio 1→4; Spearman rho = 0.83 |

## PRCC: static R0 vs dynamic outputs vs equal-width ranges

`peak_ih` (dynamic, published ranges) is the like-for-like replacement for
the static R0 PRCC; `peak_day` and `cum_infected` are additional dynamic
outputs that have no static counterpart.

### peak Ih

| Parameter | PRCC static R0 at 28 °C (rank) | PRCC dynamic, published widths (rank) | PRCC dynamic, ±20 % widths (rank) |
|---|---:|---:|---:|
| nv_ratio | +0.954 (1) | +0.950 (1) | +0.827 (4) |
| beta_vh | +0.775 (2) | +0.757 (2) | +0.829 (3) |
| beta_hv | +0.762 (3) | +0.737 (3) | +0.811 (5) |
| gamma_h | -0.639 (5) | -0.729 (4) | -0.943 (2) |
| b | +0.573 (6) | +0.553 (5) | +0.946 (1) |
| mu_v | -0.722 (4) | -0.476 (6) | -0.623 (7) |
| eip_days | -0.354 (7) | -0.394 (7) | -0.454 (8) |
| sigma_h | +0.000 (8) | +0.263 (8) | +0.637 (6) |

Spearman rho between published-width and equal-width PRCC values: 0.833; largest rank change: 4 places; identical ranking: False.

### 3-yr cumulative infection

| Parameter | PRCC dynamic, published widths (rank) | PRCC dynamic, ±20 % widths (rank) |
|---|---:|---:|
| nv_ratio | +0.955 (1) | +0.786 (4) |
| beta_vh | +0.777 (2) | +0.792 (3) |
| beta_hv | +0.758 (3) | +0.771 (6) |
| mu_v | -0.645 (4) | -0.866 (2) |
| b | +0.571 (5) | +0.934 (1) |
| gamma_h | -0.559 (6) | -0.773 (5) |
| eip_days | -0.410 (7) | -0.322 (7) |
| sigma_h | +0.117 (8) | +0.010 (8) |

Spearman rho between published-width and equal-width PRCC values: 0.833; largest rank change: 4 places; identical ranking: False.

### peak day

| Parameter | PRCC dynamic, published widths (rank) | PRCC dynamic, ±20 % widths (rank) |
|---|---:|---:|
| nv_ratio | -0.399 (1) | -0.902 (4) |
| beta_hv | -0.247 (2) | -0.903 (3) |
| beta_vh | -0.243 (3) | -0.906 (2) |
| sigma_h | -0.144 (4) | -0.694 (7) |
| eip_days | +0.127 (5) | +0.641 (8) |
| b | -0.122 (6) | -0.975 (1) |
| gamma_h | +0.055 (7) | +0.726 (5) |
| mu_v | +0.032 (8) | +0.710 (6) |

Spearman rho between published-width and equal-width PRCC values: 0.595; largest rank change: 5 places; identical ranking: False.

## Task B — Nv/Nh scaling

| Nv/Nh | peak Ih | peak prevalence | peak day | 3-yr cumulative | R0 (28 °C) | vs Nv/Nh=10 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 5 | 0.05% | 234 | 1.84% | 0.93 | -99.7% |
| 2 | 25 | 0.25% | 268 | 10.36% | 1.31 | -98.5% |
| 3 | 120 | 1.20% | 278 | 31.47% | 1.60 | -92.8% |
| 5 | 707 | 7.07% | 267 | 79.84% | 2.07 | -57.2% |
| 10 | 1651 | 16.51% | 234 | 96.89% | 2.93 | +0.0% |
| 15 | 2036 | 20.36% | 220 | 96.99% | 3.59 | +23.3% |
| 20 | 2258 | 22.58% | 211 | 97.06% | 4.14 | +36.8% |

## Task C — biting-rate clamp

| Variant | peak Ih | peak day | 3-yr cumulative | R0(14.2 °C) | R0(28 °C) | T where R0=1 | ∆peak vs published |
|---|---:|---:|---:|---:|---:|---:|---:|
| published_clip_0.18_0.23 | 1650.7568 | 233.85 | 0.968944 | 0.409 | 2.928 | 19.387 °C | +0.00e+00 % |
| floor_only_max0 | 1650.7568 | 233.85 | 0.968941 | 0.353 | 2.928 | 19.459 °C | -2.22e-07 % |
| unclamped_linear | 1650.7568 | 233.85 | 0.968941 | 0.353 | 2.928 | 19.459 °C | -2.22e-07 % |

Baseline seasonality: T spans 14.2–29.8 °C; below 21 °C 46.0% of the 3,650 days; lower clamp active 41.6% (T < 19.93 °C); upper clamp active 0.0% (T > 31.56 °C, never reached).

Numerical noise floor (rtol 1e-8 → 1e-10): peak Ih moves by 1.5e-06 persons, i.e. 9.2e-10 relative. The clamp-induced change (3.7e-6 persons, 2.2e-7 %) is of the same order, i.e. far below any meaningful resolution.

The lower clamp is active 41.6% of the year, but only inside the subcritical winter window: R0(14.2 degC) is 0.409 with the clamp and 0.353 without it, both far below 1. Removing the clamp therefore shifts peak Ih by 3.7e-6 persons (2.2e-7 %) and the cumulative attack rate by 3.4e-6 - seven orders of magnitude below the precision at which the results are reported, i.e. the annual trajectory is insensitive to the clamp. Its only measurable consequence is on the threshold: R0(T)=1 at 19.39 degC with the clamp vs 19.46 degC without. The upper clamp is never active under baseline seasonal forcing (T_max = 29.8 degC < 31.56 degC), so it has no effect at all.

## Task E — constant-parameter control (separation of mechanisms)

Same 2,000 LHS samples and the same three outputs as task A, but b and sigma_v are held constant (b_eff = b_sample, sigma_v_eff = 1/eip_sample). This makes the decomposition exact rather than approximate, because the published static analysis uses the *same* constant parameterisation (sensitivity.py overrides b and sigma_v with the sampled constants before evaluating R0). Hence:

- **static → E** changes only the *output* (algebraic R0 at 28 °C → dynamic 1,095-day simulation); b and sigma_v are constant in both.
- **E → A** changes only the *forcing* (constant → seasonally temperature-dependent b and sigma_v); the sampled values, the ranges, the scenario and the estimator are identical.

| Parameter | static R0 | constant dynamics (E) | scaled dynamics (A) |
|---|---:|---:|---:|
| nv_ratio | +0.954 (1) | +0.945 (1) | +0.950 (1) |
| beta_vh | +0.775 (2) | +0.740 (3) | +0.757 (2) |
| beta_hv | +0.762 (3) | +0.712 (4) | +0.737 (3) |
| gamma_h | -0.639 (5) | -0.799 (2) | -0.729 (4) |
| b | +0.573 (6) | +0.530 (5) | +0.553 (5) |
| mu_v | -0.722 (4) | -0.520 (6) | -0.476 (6) |
| eip_days | -0.354 (7) | -0.390 (7) | -0.394 (7) |
| sigma_h | +0.000 (8) | +0.219 (8) | +0.263 (8) |

`peak Ih`, ranked by |PRCC|.

| Output | Σ\|ΔPRCC\| static→E | Σ\|ΔPRCC\| E→A | max \|Δ\| E→A | signs preserved E→A | dominant contrast |
|---|---:|---:|---:|---|---|
| peak_ih | 0.754 | 0.230 | 0.069 | yes | output definition + dynamics |
| cum_infected | 0.310 | 0.601 | 0.189 | no | temperature forcing of b and sigma_v |
| peak_day | 6.958 | 0.877 | 0.296 | yes | output definition + dynamics |

**Verdict.** For the peak, the static→dynamic change is dominated by replacing the algebraic R0 with a dynamic output (Σ|Δ| = 0.75 vs 0.23); switching the temperature forcing on changes each PRCC by at most 0.07 and preserves every sign. For the 3-year cumulative output the ordering is reversed: temperature forcing is the larger of the two (Σ|Δ| = 0.60 vs 0.31). So neither label alone explains the published→dynamic shift; the answer depends on which output is being ranked.

Caveat: the static column for `peak_day` is not interpretable (a static R0 has no timing), so that row of Σ|Δ| is shown for completeness only.

Caveat: sigma_h is numerically negligible under all three designs (static 0.000, constant +0.219, scaled +0.263) and **changes sign** in the cumulative-output comparison (constant −0.073 vs scaled +0.117), so no interpretation should rest on it.

## Reference: published static PRCC (recomputed, R0 at 28 °C)

| Rank | Parameter | PRCC | p |
|---:|---|---:|---:|
| 1 | nv_ratio | +0.954 | <1e-300 |
| 2 | beta_vh | +0.775 | <1e-300 |
| 3 | beta_hv | +0.762 | <1e-300 |
| 4 | mu_v | -0.722 | <1e-300 |
| 5 | gamma_h | -0.639 | 6.28e-230 |
| 6 | b | +0.573 | 1.23e-174 |
| 7 | eip_days | -0.354 | 4.40e-60 |
| 8 | sigma_h | +0.000 | 9.87e-01 |

R0 across the 2,000 samples: 0.48–4.89 (median 2.06) at 28 °C.

## Findings

1. **Temperature dependence does not overturn the published ordering of the dominant parameters.** Against the dynamic peak, the top three parameters and every sign are the same as in the static R0 analysis (nv_ratio > beta_vh > beta_hv, |static − dynamic| ≤ 0.03 for those three). The *secondary* parameters do move: mu_v goes from −0.722 to −0.476 (static → dynamic) and sigma_h from 0.000 (n.s.) to +0.263 (p < 1e-30). Saying "the published ranking is unchanged" would be wrong; saying "the dominant group is unchanged" is correct.
2. **Timing has a different sensitivity structure from magnitude.** For peak day, nv_ratio still dominates (−0.399) but every other parameter is |PRCC| < 0.25, and mu_v is not significant. Peak-based rankings must not be extrapolated to epidemic timing.
3. **Range width matters for the exact order, not for the signs.** With ±20 % ranges on all eight parameters, b becomes the strongest driver of the peak (0.946) and nv_ratio drops to fourth (0.827, still large); all eight parameters then land in |PRCC| = 0.45–0.95 with no single dominant driver. So "nv_ratio is the single most influential parameter" is partly a consequence of its 1–20 range, whereas "nv_ratio and the two transmission probabilities are the dominant group" survives the equal-width test.
4. **The clamp is inert for the dynamics, not for the threshold.** Removing the lower clamp changes the peak by 2e-7 % and the cumulative attack rate by 3e-6; only the R0(T)=1 crossing moves, from 19.39 to 19.46 °C. The upper clamp is never activated (T_max = 29.8 °C).
5. **The 1,650 figure is conditional on Nv/Nh = 10.** Across Nv/Nh = 1–20 the peak spans 5–2,258, while the 3-year attack rate saturates at ≈97 % once Nv/Nh ≥ 10 because the susceptible pool is exhausted within the 3-year window.

## Caveats

- b and sigma_v enter the dynamic analysis as *level scalings* of the temperature response (shape preserved), not as fixed constants. The dynamic-vs-static difference therefore mixes "temperature dependence" with "scaling vs constant substitution".
- In the published-width design 5 % of samples are subcritical (3-year cumulative infection < 1.5 %), where the monotonicity assumption behind PRCC is only approximate. The equal-width design contains no subcritical sample (minimum cumulative 38.9 %), which is part of the reason its |PRCC| values are uniformly larger.
- peak day is refined by a three-point parabola around the daily maximum; on the raw daily grid the peak day takes only a few distinct values and the tied ranks would attenuate its PRCC substantially.
- The two cumulative definitions differ by the imported cases injected directly into Ih: mean difference 0.0078 in proportion.
- No time-resolved (windowed) PRCC was computed; if a reader asks how sensitivity evolves within a season, that is a separate analysis.

## Reproducibility

- Python 3.13.12, numpy 2.4.6, scipy 1.18.0
- seed 20260816 everywhere; task A and D use identical LHS uniform variates (paired design, `lhs_variates_paired_with_task_A = True`)
- n = 2000 for tasks A and D (full 2,000; per-solve 0.112 s, wall 64 s on 7 processes)
- static reference recomputed with `sensitivity.run()` defaults (n=2000, T=28 °C)
- source fingerprint: `model.py` 8664b17fdfe05dd7, `sensitivity.py` 1a15820f734af00f; numeric re-check against the values used in the runs: True
- `../src` was touched by another process at 12:50 local time; the numeric content was verified identical (baseline peak Ih 1650.756836 at n=1096, 1649.777450 on the default 1000-point grid, R0(28 °C) 2.928116).

| File | Content |
|---|---|
| `dynamic_prcc.json` | task A: 8 × 3 PRCC table, p-values, output distributions, meta |
| `nv_scaling.json` / `nv_scaling.md` | task B: Nv/Nh sweep |
| `clamp_sensitivity.json` | task C: three clamp variants, T_crit, regime stats, noise floor |
| `equalwidth_prcc.json` | task D: ±20 % PRCC + comparison with task A |
| `constant_param_prcc.json` | task E: constant-b/sigma_v control + three-way comparison |
| `static_prcc_reference.json` | published static PRCC, recomputed for comparison |
| `dv_common.py` | shared solver/metrics/PRCC helpers |
| `run_dynamic_prcc.py`, `run_scaling.py`, `run_clamp.py`, `run_equalwidth.py`, `make_static_reference.py`, `stamp_src_hashes.py`, `make_summary.py` | scripts |

