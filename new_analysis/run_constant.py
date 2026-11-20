"""Control run: PRCC with b and sigma_v held CONSTANT (no temperature forcing).

Purpose
-------
Task A replaced the temperature-dependent rates by level scalings, so its
difference from the published static PRCC mixes two effects.  This run removes
the ambiguity with a third design that differs from task A in *exactly one*
respect:

    task A  (scaled)    b_eff(T)      = b_model(T) * (b_sample / b(28 C))
                        sigma_v_eff(T)= sigma_v_model(T) * (sv_sample / sv(28 C))
    task E  (constant)  b_eff         = b_sample              (constant)
                        sigma_v_eff   = 1 / eip_days_sample   (constant)

Everything else is identical: the same 2,000 LHS samples (seed 20260816),
the same baseline scenario, the same 1,096-point daily output grid, the same
PRCC estimator, and the same param->model mapping for the remaining six
parameters.  All temperature dependence in the *other* parameters is absent by
construction (none of them is temperature-dependent in this model), so the
constant-vs-scaled contrast isolates "does b/sigma_v respond to the seasonal
temperature cycle?".

The comparison against the published static R0 PRCC therefore decomposes as

    static -> constant   : output definition (algebraic R0 vs dynamic peak)
                           + dynamical non-linearity
    constant -> scaled   : temperature forcing of b and sigma_v   (exact pair)
"""

import json
import os
import sys
import time
import multiprocessing as mp

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

import dv_common as C
import model as M
import sensitivity as S

SEED = 20260816
N_SAMPLES = 2000


def build_jobs(samples):
    jobs = []
    for row in samples:
        overrides = {
            'nv_ratio': float(row[0]),
            'beta_hv': float(row[1]),
            'beta_vh': float(row[2]),
            'mu_v': float(row[3]),
            'gamma_h': float(row[4]),
            'sigma_h': float(row[7]),
        }
        jobs.append((float(row[5]), 1.0 / float(row[6]), overrides))
    return jobs


def five(x):
    x = x[np.isfinite(x)]
    return {'min': float(np.min(x)), 'p05': float(np.percentile(x, 5)),
            'median': float(np.median(x)), 'p95': float(np.percentile(x, 95)),
            'max': float(np.max(x)), 'mean': float(np.mean(x))}


def load_json(name):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return None
    with open(path, encoding='utf-8') as fh:
        return json.load(fh)


def compare(this_rows, other_map, this_label, other_label):
    """Per-parameter comparison of two PRCC tables (same parameter set)."""
    per = []
    for r in this_rows:
        o = other_map[r['param']]
        per.append({
            'param': r['param'],
            'prcc_%s' % this_label: r['prcc'],
            'prcc_%s' % other_label: o['prcc'],
            'delta': r['prcc'] - o['prcc'],
            'rank_%s' % this_label: r['rank'],
            'rank_%s' % other_label: o['rank'],
            'rank_change': r['rank'] - o['rank'],
        })
    a = [p['prcc_%s' % this_label] for p in per]
    b = [p['prcc_%s' % other_label] for p in per]
    signs_same = all(np.sign(x) == np.sign(y) if (x != 0 and y != 0) else True
                     for x, y in zip(a, b))
    return {
        'per_param': per,
        'max_abs_delta_prcc': float(max(abs(p['delta']) for p in per)),
        'mean_abs_delta_prcc': float(np.mean([abs(p['delta']) for p in per])),
        'sum_abs_delta_prcc': float(np.sum([abs(p['delta']) for p in per])),
        'signs_identical': bool(signs_same),
        'max_abs_rank_change': int(max(abs(p['rank_change']) for p in per)),
        'identical_ranking': bool(all(p['rank_change'] == 0 for p in per)),
        'spearman_of_prcc': float(stats.spearmanr(a, b).statistic),
    }


def main():
    t_start = time.perf_counter()
    C.reset_rates()

    t0 = time.perf_counter()
    C.run_case(M.baseline_params())
    per_run = time.perf_counter() - t0

    samples = S.sample_lhs(n_samples=N_SAMPLES, seed=SEED)
    jobs = build_jobs(samples)

    n_proc = max(1, min(7, mp.cpu_count() - 1))
    print(f"[run] {len(jobs)} constant-parameter solves on {n_proc} processes ...")
    if n_proc > 1:
        with mp.Pool(n_proc) as pool:
            results = pool.map(C.worker_constant, jobs, chunksize=16)
    else:
        results = [C.worker_constant(j) for j in jobs]
    wall = time.perf_counter() - t_start
    n_fail = sum(0 if r['ok'] else 1 for r in results)
    print(f"[run] done in {wall:.1f} s ({n_fail} failed)")

    peak_ih = np.array([r['peak_ih'] for r in results])
    peak_day = np.array([r['peak_day'] for r in results])
    cum = np.array([r['cum_infected'] for r in results])

    blocks = {}
    for out_name, arr in (('peak_ih', peak_ih), ('peak_day', peak_day),
                          ('cum_infected', cum)):
        rows, dropped = C.prcc_table(samples, arr, S.PARAM_NAMES)
        blocks[out_name] = {'prcc_table': rows, 'n_dropped_nonfinite': dropped}
        print(f"\n--- constant-parameter PRCC vs {out_name} (n={rows[0]['n_eff']}) ---")
        for r in rows:
            sig = '***' if r['p_value'] < 0.001 else (
                '**' if r['p_value'] < 0.01 else ('*' if r['p_value'] < 0.05 else 'n.s.'))
            print(f"{r['rank']:<3}{r['param']:<12}{r['prcc']:>9.3f}"
                  f"{C.fmt_p(r['p_value']):>12}  {sig}")

    # ── comparisons ────────────────────────────────────────────────────────
    dA = load_json('dynamic_prcc.json')
    dS = load_json('static_prcc_reference.json')

    comparisons = {}
    for out_name in blocks:
        this_rows = blocks[out_name]['prcc_table']
        entry = {}
        if dA is not None:
            a_map = {r['param']: r for r in dA['outputs'][out_name]['prcc_table']}
            entry['vs_scaled_dynamic_task_A'] = compare(
                this_rows, a_map, 'constant', 'scaled')
        if dS is not None:
            s_rows = [{'param': n, 'prcc': p,
                       'rank': i + 1} for i, (n, p) in enumerate(
                zip(dS['param_names'], dS['prcc']))]
            s_map = {r['param']: r for r in s_rows}
            entry['vs_static_R0'] = compare(this_rows, s_map, 'constant', 'static')
        comparisons[out_name] = entry

    # ── which contrast dominates? ──────────────────────────────────────────
    verdict = {}
    if dA is not None and dS is not None:
        for out_name in blocks:
            c = comparisons[out_name]
            d_temp = c['vs_scaled_dynamic_task_A']['sum_abs_delta_prcc']
            d_dyn = c['vs_static_R0']['sum_abs_delta_prcc']
            verdict[out_name] = {
                'sum_abs_delta_static_to_constant': d_dyn,
                'sum_abs_delta_constant_to_scaled': d_temp,
                'larger_contribution': ('temperature forcing of b and sigma_v'
                                        if d_temp > d_dyn
                                        else 'output definition + dynamics'),
                'ratio_temp_over_dynamics': float(d_temp / d_dyn) if d_dyn else None,
            }

    payload = {
        'task': 'E - control run, PRCC with b and sigma_v held constant',
        'method': (
            'Identical to task A (LHS n=2000, seed 20260816, baseline '
            'scenario, t_span=(0,1095), daily output, Marino PRCC) except '
            'that b_eff = b_sample and sigma_v_eff = 1/eip_days_sample are '
            'constants with no temperature dependence.'
        ),
        'design_contrast': (
            'Task E differs from task A in exactly one respect: whether b and '
            'sigma_v follow the seasonal temperature cycle. Their sampled '
            'values (and hence their nominal ranges) are identical, and every '
            'other parameter enters through the same mapping, so the '
            'constant-vs-scaled contrast isolates the effect of temperature '
            'forcing on the sensitivity ranking. The static-vs-constant '
            'contrast instead isolates the effect of replacing the algebraic '
            'R0 by a dynamic output (plus the model non-linearity).'
        ),
        'reported_outputs': ['peak_ih', 'cum_infected', 'peak_day'],
        'outputs': blocks,
        'output_distributions': {'peak_ih': five(peak_ih),
                                 'peak_day': five(peak_day),
                                 'cum_infected': five(cum)},
        'comparisons': comparisons,
        'decomposition': verdict,
        'n_failed_solves': n_fail,
        'meta': C.make_meta(
            script='run_constant.py', seed=SEED, n_samples=N_SAMPLES,
            wall_seconds=wall, per_run_seconds=per_run,
            extra={'parallel_processes': n_proc,
                   'b_ref_28C': C.B_REF, 'sigma_v_ref_28C': C.SV_REF}
        ),
    }
    C.write_json(os.path.join(HERE, 'constant_param_prcc.json'), payload)


if __name__ == '__main__':
    main()
