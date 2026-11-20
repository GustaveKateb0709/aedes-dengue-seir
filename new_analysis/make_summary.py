"""Assemble SUMMARY.md from the four result files (no new computation)."""

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with open(os.path.join(HERE, name), encoding='utf-8') as fh:
        return json.load(fh)


def C_fmt(p):
    return '<1e-300' if p < 1e-300 else f'{p:.2e}'


def main():
    dA = load('dynamic_prcc.json')
    dB = load('nv_scaling.json')
    dC = load('clamp_sensitivity.json')
    dD = load('equalwidth_prcc.json')
    dS = load('static_prcc_reference.json')

    L = []
    a = L.append

    a('# Dynamic-output sensitivity analyses — summary')
    a('')
    a('All runs use the unmodified model in `../src`; the only intervention is a')
    a('run-time replacement of `model.biting_rate` / `model.eip_rate`, never an edit')
    a('of the source file. Baseline scenario = seasonal temperature + annual import')
    a('pulse, `t_span=(0, 1095)`, 1096 daily output points, `Nh = 10,000`.')
    a('')

    # ── A/B/C/D headline table ─────────────────────────────────────────────
    a('## Headline numbers')
    a('')
    a('| Task | Question | Key number |')
    a('|---|---|---|')
    r0 = dA['baseline_reference_run']
    a(f"| A | PRCC of the *dynamic* peak vs the static PRCC | baseline peak Ih = "
      f"{r0['peak_ih']:.1f} at day {r0['peak_day']:.0f}; nv_ratio PRCC 0.950 "
      f"(static 0.954), b 0.553 (static 0.573) |")
    b10 = next(r for r in dB['rows'] if r['nv_ratio'] == 10)
    b1 = next(r for r in dB['rows'] if r['nv_ratio'] == 1)
    b20 = next(r for r in dB['rows'] if r['nv_ratio'] == 20)
    a(f"| B | How does the 1,650 headline scale with Nv/Nh? | peak Ih spans "
      f"{b1['peak_ih']:.0f} (Nv/Nh=1) to {b20['peak_ih']:.0f} (Nv/Nh=20) — a "
      f"{b20['peak_ih'] / max(b1['peak_ih'], 1e-9):.0f}× span; the 3-yr attack "
      f"rate saturates at {b10['cum_infected']:.1%} for Nv/Nh ≥ 10 |")
    pub = dC['variants'][0]
    unc = dC['variants'][2]
    a(f"| C | What if the b(T) clamp is removed? | peak Ih changes by "
      f"{unc['peak_pct_change_vs_published']:.1e} % (i.e. 3.7e-6 persons); the "
      f"only measurable effect is on the threshold, "
      f"{pub['critical_temperature_R0_eq_1']:.2f} → "
      f"{unc['critical_temperature_R0_eq_1']:.2f} °C |")
    cmp_peak = dD['comparison_vs_task_A']['peak_ih']
    def rk(out_key, name, key):
        for p in dD['comparison_vs_task_A'][out_key]['per_param']:
            if p['param'] == name:
                return p[key]
    a(f"| D | Is the ranking an artefact of range width? | **Partly.** Signs are "
      f"all preserved, but the exact order is not: b moves "
      f"{rk('peak_ih', 'b', 'rank_published_widths')}→"
      f"{rk('peak_ih', 'b', 'rank_equal_widths')} and nv_ratio "
      f"{rk('peak_ih', 'nv_ratio', 'rank_published_widths')}→"
      f"{rk('peak_ih', 'nv_ratio', 'rank_equal_widths')}; Spearman rho = "
      f"{cmp_peak['spearman_rank_correlation_of_prcc']:.2f} |")
    a('')

    # ── PRCC comparison ────────────────────────────────────────────────────
    a('## PRCC: static R0 vs dynamic outputs vs equal-width ranges')
    a('')
    a('`peak_ih` (dynamic, published ranges) is the like-for-like replacement for')
    a('the static R0 PRCC; `peak_day` and `cum_infected` are additional dynamic')
    a('outputs that have no static counterpart.')
    a('')
    static_map = {n: (p, i + 1) for i, (n, p) in enumerate(
        zip(dS['param_names'], dS['prcc']))}

    for out_key, title in (('peak_ih', 'peak Ih'),
                           ('cum_infected', '3-yr cumulative infection'),
                           ('peak_day', 'peak day')):
        a(f'### {title}')
        a('')
        dmap = {r['param']: r for r in dA['outputs'][out_key]['prcc_table']}
        emap = {r['param']: r for r in dD['outputs'][out_key]['prcc_table']}
        order = [r['param'] for r in dA['outputs'][out_key]['prcc_table']]
        if out_key == 'peak_ih':
            a('| Parameter | PRCC static R0 at 28 °C (rank) | PRCC dynamic, '
              'published widths (rank) | PRCC dynamic, ±20 % widths (rank) |')
            a('|---|---:|---:|---:|')
            for name in order:
                sp, sr = static_map[name]
                a(f"| {name} | {sp:+.3f} ({sr}) | {dmap[name]['prcc']:+.3f} "
                  f"({dmap[name]['rank']}) | {emap[name]['prcc']:+.3f} "
                  f"({emap[name]['rank']}) |")
        else:
            a('| Parameter | PRCC dynamic, published widths (rank) | '
              'PRCC dynamic, ±20 % widths (rank) |')
            a('|---|---:|---:|')
            for name in order:
                a(f"| {name} | {dmap[name]['prcc']:+.3f} ({dmap[name]['rank']}) | "
                  f"{emap[name]['prcc']:+.3f} ({emap[name]['rank']}) |")
        a('')
        a(f"Spearman rho between published-width and equal-width PRCC values: "
          f"{dD['comparison_vs_task_A'][out_key]['spearman_rank_correlation_of_prcc']:.3f}; "
          f"largest rank change: "
          f"{dD['comparison_vs_task_A'][out_key]['max_abs_rank_change']} places; "
          f"identical ranking: "
          f"{dD['comparison_vs_task_A'][out_key]['identical_ranking']}.")
        a('')

    # ── B table ────────────────────────────────────────────────────────────
    a('## Task B — Nv/Nh scaling')
    a('')
    a('| Nv/Nh | peak Ih | peak prevalence | peak day | 3-yr cumulative | R0 (28 °C) | vs Nv/Nh=10 |')
    a('|---:|---:|---:|---:|---:|---:|---:|')
    for r in dB['rows']:
        a(f"| {r['nv_ratio']:.0f} | {r['peak_ih']:.0f} | {r['peak_ih'] / 10000:.2%} | "
          f"{r['peak_day']:.0f} | {r['cum_infected']:.2%} | {r['R0_at_28C']:.2f} | "
          f"{r['peak_pct_change_vs_nv10']:+.1f}% |")
    a('')

    # ── C table ────────────────────────────────────────────────────────────
    a('## Task C — biting-rate clamp')
    a('')
    a('| Variant | peak Ih | peak day | 3-yr cumulative | R0(14.2 °C) | R0(28 °C) | T where R0=1 | ∆peak vs published |')
    a('|---|---:|---:|---:|---:|---:|---:|---:|')
    for r in dC['variants']:
        tc = r['critical_temperature_R0_eq_1']
        a(f"| {r['variant']} | {r['peak_ih']:.4f} | {r['peak_day']:.2f} | "
          f"{r['cum_infected']:.6f} | {r['R0_probe']['14.2']:.3f} | "
          f"{r['R0_probe']['28.0']:.3f} | "
          f"{'n/a' if tc is None else f'{tc:.3f} °C'} | "
          f"{r['peak_pct_change_vs_published']:+.2e} % |")
    a('')
    fs = dC['temperature_regime_stats']
    a(f"Baseline seasonality: T spans {fs['T_min']:.1f}–{fs['T_max']:.1f} °C; "
      f"below 21 °C {fs['frac_T_below_21']:.1%} of the 3,650 days; "
      f"lower clamp active {fs['frac_b_lower_clamp']:.1%} (T < "
      f"{dC['clamp_thresholds']['T_lower_degC']:.2f} °C); upper clamp active "
      f"{fs['frac_b_upper_clamp']:.1%} (T > "
      f"{dC['clamp_thresholds']['T_upper_degC']:.2f} °C, never reached).")
    a('')
    nf = dC['numerical_noise_floor']
    a(f"Numerical noise floor (rtol 1e-8 → 1e-10): peak Ih moves by "
      f"{abs(nf['peak_ih_default_minus_tight']):.1e} persons, i.e. "
      f"{nf['peak_ih_relative']:.1e} relative. The clamp-induced change "
      f"(3.7e-6 persons, 2.2e-7 %) is of the same order, i.e. far below any "
      f"meaningful resolution.")
    a('')
    a(dC['interpretation'])
    a('')

    # ── task E: constant-parameter control ─────────────────────────────────
    dE_path = os.path.join(HERE, 'constant_param_prcc.json')
    if os.path.exists(dE_path):
        with open(dE_path, encoding='utf-8') as fh:
            dE = json.load(fh)
        a('## Task E — constant-parameter control (separation of mechanisms)')
        a('')
        a('Same 2,000 LHS samples and the same three outputs as task A, but b '
          'and sigma_v are held constant (b_eff = b_sample, sigma_v_eff = '
          '1/eip_sample). This makes the decomposition exact rather than '
          'approximate, because the published static analysis uses the *same* '
          'constant parameterisation (sensitivity.py overrides b and sigma_v '
          'with the sampled constants before evaluating R0). Hence:')
        a('')
        a('- **static → E** changes only the *output* (algebraic R0 at 28 °C → '
          'dynamic 1,095-day simulation); b and sigma_v are constant in both.')
        a('- **E → A** changes only the *forcing* (constant → seasonally '
          'temperature-dependent b and sigma_v); the sampled values, the '
          'ranges, the scenario and the estimator are identical.')
        a('')
        a('| Parameter | static R0 | constant dynamics (E) | scaled dynamics (A) |')
        a('|---|---:|---:|---:|')
        cmap = {r['param']: r for r in dE['outputs']['peak_ih']['prcc_table']}
        amap = {r['param']: r for r in dA['outputs']['peak_ih']['prcc_table']}
        for name in [r['param'] for r in dA['outputs']['peak_ih']['prcc_table']]:
            sp, sr = static_map[name]
            a(f"| {name} | {sp:+.3f} ({sr}) | {cmap[name]['prcc']:+.3f} "
              f"({cmap[name]['rank']}) | {amap[name]['prcc']:+.3f} "
              f"({amap[name]['rank']}) |")
        a('')
        a('`peak Ih`, ranked by |PRCC|.')
        a('')
        a('| Output | Σ\\|ΔPRCC\\| static→E | Σ\\|ΔPRCC\\| E→A | max \\|Δ\\| E→A | signs preserved E→A | dominant contrast |')
        a('|---|---:|---:|---:|---|---|')
        for out in ('peak_ih', 'cum_infected', 'peak_day'):
            dec = dE['decomposition'][out]
            cmp_ea = dE['comparisons'][out]['vs_scaled_dynamic_task_A']
            a(f"| {out} | {dec['sum_abs_delta_static_to_constant']:.3f} | "
              f"{dec['sum_abs_delta_constant_to_scaled']:.3f} | "
              f"{cmp_ea['max_abs_delta_prcc']:.3f} | "
              f"{'yes' if cmp_ea['signs_identical'] else 'no'} | "
              f"{dec['larger_contribution']} |")
        a('')
        a('**Verdict.** For the peak, the static→dynamic change is dominated by '
          'replacing the algebraic R0 with a dynamic output '
          f"(Σ|Δ| = {dE['decomposition']['peak_ih']['sum_abs_delta_static_to_constant']:.2f} "
          'vs '
          f"{dE['decomposition']['peak_ih']['sum_abs_delta_constant_to_scaled']:.2f}); "
          'switching the temperature forcing on changes each PRCC by at most '
          f"{dE['comparisons']['peak_ih']['vs_scaled_dynamic_task_A']['max_abs_delta_prcc']:.2f} "
          'and preserves every sign. For the 3-year cumulative output the '
          'ordering is reversed: temperature forcing is the larger of the two '
          f"(Σ|Δ| = {dE['decomposition']['cum_infected']['sum_abs_delta_constant_to_scaled']:.2f} "
          'vs '
          f"{dE['decomposition']['cum_infected']['sum_abs_delta_static_to_constant']:.2f}). "
          'So neither label alone explains the published→dynamic shift; the '
          'answer depends on which output is being ranked.')
        a('')
        a('Caveat: the static column for `peak_day` is not interpretable (a '
          'static R0 has no timing), so that row of Σ|Δ| is shown for '
          'completeness only.')
        a('')
        a('Caveat: sigma_h is numerically negligible under all three designs '
          '(static 0.000, constant +0.219, scaled +0.263) and **changes sign** '
          'in the cumulative-output comparison (constant −0.073 vs scaled '
          '+0.117), so no interpretation should rest on it.')
        a('')

    # ── published static reference ─────────────────────────────────────────
    a('## Reference: published static PRCC (recomputed, R0 at 28 °C)')
    a('')
    a('| Rank | Parameter | PRCC | p |')
    a('|---:|---|---:|---:|')
    for i, name in enumerate(dS['param_names']):
        a(f"| {i + 1} | {name} | {dS['prcc'][i]:+.3f} | {C_fmt(dS['p_value'][i])} |")
    a('')
    a(f"R0 across the 2,000 samples: {dS['r0_min']:.2f}–{dS['r0_max']:.2f} "
      f"(median {dS['r0_median']:.2f}) at 28 °C.")
    a('')

    # ── findings ───────────────────────────────────────────────────────────
    a('## Findings')
    a('')
    a('1. **Temperature dependence does not overturn the published ordering of '
      'the dominant parameters.** Against the dynamic peak, the top three '
      'parameters and every sign are the same as in the static R0 analysis '
      '(nv_ratio > beta_vh > beta_hv, |static − dynamic| ≤ 0.03 for those '
      'three). The *secondary* parameters do move: mu_v goes from −0.722 to '
      f"−0.476 (static → dynamic) and sigma_h from 0.000 (n.s.) to +0.263 "
      '(p < 1e-30). Saying "the published ranking is unchanged" would be '
      'wrong; saying "the dominant group is unchanged" is correct.')
    a('2. **Timing has a different sensitivity structure from magnitude.** For '
      'peak day, nv_ratio still dominates (−0.399) but every other parameter '
      'is |PRCC| < 0.25, and mu_v is not significant. Peak-based rankings '
      'must not be extrapolated to epidemic timing.')
    a('3. **Range width matters for the exact order, not for the signs.** With '
      '±20 % ranges on all eight parameters, b becomes the strongest driver of '
      'the peak (0.946) and nv_ratio drops to fourth (0.827, still large); all '
      'eight parameters then land in |PRCC| = 0.45–0.95 with no single '
      'dominant driver. So "nv_ratio is the single most influential parameter" '
      'is partly a consequence of its 1–20 range, whereas "nv_ratio and the '
      'two transmission probabilities are the dominant group" survives the '
      'equal-width test.')
    a('4. **The clamp is inert for the dynamics, not for the threshold.** '
      'Removing the lower clamp changes the peak by 2e-7 % and the cumulative '
      'attack rate by 3e-6; only the R0(T)=1 crossing moves, from 19.39 to '
      '19.46 °C. The upper clamp is never activated (T_max = 29.8 °C).')
    a('5. **The 1,650 figure is conditional on Nv/Nh = 10.** Across Nv/Nh = '
      '1–20 the peak spans 5–2,258, while the 3-year attack rate saturates at '
      '≈97 % once Nv/Nh ≥ 10 because the susceptible pool is exhausted within '
      'the 3-year window.')
    a('')

    a('## Caveats')
    a('')
    a('- b and sigma_v enter the dynamic analysis as *level scalings* of the '
      'temperature response (shape preserved), not as fixed constants. The '
      'dynamic-vs-static difference therefore mixes "temperature dependence" '
      'with "scaling vs constant substitution".')
    a('- In the published-width design 5 % of samples are subcritical (3-year '
      'cumulative infection < 1.5 %), where the monotonicity assumption behind '
      'PRCC is only approximate. The equal-width design contains no '
      'subcritical sample (minimum cumulative 38.9 %), which is part of the '
      'reason its |PRCC| values are uniformly larger.')
    a('- peak day is refined by a three-point parabola around the daily '
      'maximum; on the raw daily grid the peak day takes only a few distinct '
      'values and the tied ranks would attenuate its PRCC substantially.')
    a('- The two cumulative definitions differ by the imported cases injected '
      'directly into Ih: mean difference 0.0078 in proportion.')
    a('- No time-resolved (windowed) PRCC was computed; if a reader asks '
      'how sensitivity evolves within a season, that is a separate analysis.')
    a('')

    # ── methods notes ──────────────────────────────────────────────────────
    a('## Reproducibility')
    a('')
    a(f"- Python {dA['meta']['python']}, numpy {dA['meta']['numpy']}, "
      f"scipy {dA['meta']['scipy']}")
    a(f"- seed {dA['meta']['random_seed']} everywhere; task A and D use identical "
      f"LHS uniform variates (paired design, `lhs_variates_paired_with_task_A = "
      f"{dD['lhs_variates_paired_with_task_A']}`)")
    a(f"- n = {dA['meta']['n_samples']} for tasks A and D "
      f"(full 2,000; per-solve {dA['meta']['seconds_per_ode_solve']:.3f} s, "
      f"wall {dA['meta']['wall_seconds_total']:.0f} s on "
      f"{dA['meta']['parallel_processes']} processes)")
    a(f"- static reference recomputed with `sensitivity.run()` defaults "
      f"(n={dS['n_samples']}, T={dS['T']:.0f} °C)")
    sha = dA['meta']['src_sha256']
    a(f"- source fingerprint: `model.py` {sha['model.py'][:16]}, "
      f"`sensitivity.py` {sha['sensitivity.py'][:16]}; numeric re-check against "
      f"the values used in the runs: "
      f"{dA['meta']['src_numeric_recheck_matches_run_values']}")
    a(f"- `../src` was touched by another process at 12:50 local time; the "
      f"numeric content was verified identical (baseline peak Ih "
      f"{dA['meta']['src_numeric_recheck']['peak_ih_n1096']:.6f} at n=1096, "
      f"{dA['meta']['src_numeric_recheck']['peak_ih_default_grid_n1000']:.6f} "
      f"on the default 1000-point grid, R0(28 °C) "
      f"{dA['meta']['src_numeric_recheck']['R0_at_28C']:.6f}).")
    a('')
    a('| File | Content |')
    a('|---|---|')
    a('| `dynamic_prcc.json` | task A: 8 × 3 PRCC table, p-values, output distributions, meta |')
    a('| `nv_scaling.json` / `nv_scaling.md` | task B: Nv/Nh sweep |')
    a('| `clamp_sensitivity.json` | task C: three clamp variants, T_crit, regime stats, noise floor |')
    a('| `equalwidth_prcc.json` | task D: ±20 % PRCC + comparison with task A |')
    a('| `constant_param_prcc.json` | task E: constant-b/sigma_v control + three-way comparison |')
    a('| `static_prcc_reference.json` | published static PRCC, recomputed for comparison |')
    a('| `dv_common.py` | shared solver/metrics/PRCC helpers |')
    a('| `run_dynamic_prcc.py`, `run_scaling.py`, `run_clamp.py`, `run_equalwidth.py`, '
      '`make_static_reference.py`, `stamp_src_hashes.py`, `make_summary.py` | scripts |')
    a('')

    out = os.path.join(HERE, 'SUMMARY.md')
    with open(out, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')
    print(f'[written] {out}')


if __name__ == '__main__':
    main()
