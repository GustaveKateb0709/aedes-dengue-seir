"""Task D: PRCC with equal-width (+/-20 %) parameter ranges.

Motivation: the published ranges differ in width by an order of magnitude (Nv/Nh 1-20 vs b 0.18-0.23), so the PRCC ranking may be driven by
range width rather than by model structure.

Here every parameter is given the same relative half-width (+/-20 % around its
baseline value) and the same dynamic outputs are recomputed with the same LHS
construction, the same sample size (n=2000) and the same seed.  Because the
LHS draws depend only on the seed and the parameter order, the uniform
variates u are *identical* to those of task A; the two designs are paired
sample-for-sample and differ only in the ranges.

Range construction
------------------
    lo = 0.8 * theta_base,  hi = 1.2 * theta_base

with theta_base the model baseline (b and eip_days are centred on their 28 degC
reference values, i.e. the same reference point used by task A).
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
HALF_WIDTH = 0.20

BASE_VALUES = {
    'nv_ratio': float(M.NV_RATIO),
    'beta_hv': float(M.BETA_HV),
    'beta_vh': float(M.BETA_VH),
    'mu_v': float(M.MU_V),
    'gamma_h': float(M.GAMMA_H),
    'b': C.B_REF,
    'eip_days': C.EIP_REF_DAYS,
    'sigma_h': float(M.SIGMA_H),
}


def equal_width_ranges():
    return {k: (v * (1.0 - HALF_WIDTH), v * (1.0 + HALF_WIDTH))
            for k, v in BASE_VALUES.items()}


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
        b_scale = float(row[5]) / C.B_REF
        sv_scale = (1.0 / float(row[6])) / C.SV_REF
        jobs.append((b_scale, sv_scale, overrides))
    return jobs


def five(x):
    x = x[np.isfinite(x)]
    return {'min': float(np.min(x)), 'p05': float(np.percentile(x, 5)),
            'median': float(np.median(x)), 'p95': float(np.percentile(x, 95)),
            'max': float(np.max(x)), 'mean': float(np.mean(x))}


def main():
    t_start = time.perf_counter()
    C.reset_rates()

    t0 = time.perf_counter()
    C.run_case(M.baseline_params())
    per_run = time.perf_counter() - t0

    ranges = equal_width_ranges()
    samples = C.sample_lhs_ranges(S.PARAM_NAMES, ranges,
                                  n_samples=N_SAMPLES, seed=SEED)

    # paired-design check: identical uniform variates as task A
    a_samples = S.sample_lhs(n_samples=N_SAMPLES, seed=SEED)
    u_a = (a_samples - np.array([S.PARAM_RANGES[n][0] for n in S.PARAM_NAMES])) / \
          np.array([S.PARAM_RANGES[n][1] - S.PARAM_RANGES[n][0] for n in S.PARAM_NAMES])
    u_d = (samples - np.array([ranges[n][0] for n in S.PARAM_NAMES])) / \
          np.array([ranges[n][1] - ranges[n][0] for n in S.PARAM_NAMES])
    paired = bool(np.allclose(u_a, u_d, atol=1e-12))
    print(f"[check] LHS variates identical to task A: {paired}")

    jobs = build_jobs(samples)
    n_proc = max(1, min(7, mp.cpu_count() - 1))
    print(f"[run] {len(jobs)} solves on {n_proc} processes ...")
    if n_proc > 1:
        with mp.Pool(n_proc) as pool:
            results = pool.map(C.worker, jobs, chunksize=16)
    else:
        results = [C.worker(j) for j in jobs]
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
        print(f"\n--- equal-width PRCC vs {out_name} (n={rows[0]['n_eff']}) ---")
        for r in rows:
            sig = '***' if r['p_value'] < 0.001 else (
                '**' if r['p_value'] < 0.01 else ('*' if r['p_value'] < 0.05 else 'n.s.'))
            print(f"{r['rank']:<3}{r['param']:<12}{r['prcc']:>9.3f}"
                  f"{C.fmt_p(r['p_value']):>12}  {sig}")

    # ── comparison with task A (published widths) ──────────────────────────
    comparison = {}
    a_path = os.path.join(HERE, 'dynamic_prcc.json')
    if os.path.exists(a_path):
        with open(a_path, encoding='utf-8') as fh:
            a_json = json.load(fh)
        for out_name in blocks:
            a_tab = a_json['outputs'][out_name]['prcc_table']
            a_map = {r['param']: r for r in a_tab}
            d_map = {r['param']: r for r in blocks[out_name]['prcc_table']}
            names = [r['param'] for r in blocks[out_name]['prcc_table']]
            a_rank = [a_map[n]['rank'] for n in names]
            d_rank = [d_map[n]['rank'] for n in names]
            a_val = [a_map[n]['prcc'] for n in names]
            d_val = [d_map[n]['prcc'] for n in names]
            comparison[out_name] = {
                'per_param': [
                    {'param': n,
                     'prcc_published_widths': a_map[n]['prcc'],
                     'prcc_equal_widths': d_map[n]['prcc'],
                     'rank_published_widths': a_map[n]['rank'],
                     'rank_equal_widths': d_map[n]['rank'],
                     'delta_prcc': d_map[n]['prcc'] - a_map[n]['prcc']}
                    for n in names
                ],
                'spearman_rank_correlation_of_prcc': float(
                    stats.spearmanr(a_val, d_val).statistic),
                'top3_published_widths': [n for n in
                                          sorted(names, key=lambda k: a_map[k]['rank'])][:3],
                'top3_equal_widths': [n for n in
                                      sorted(names, key=lambda k: d_map[k]['rank'])][:3],
                'identical_ranking': bool(a_rank == d_rank),
                'max_abs_rank_change': int(max(abs(x - y) for x, y in zip(a_rank, d_rank))),
            }
    else:
        comparison = {'note': 'dynamic_prcc.json not found; run task A first'}

    payload = {
        'task': 'D - PRCC with equal-width (+/-20 %) parameter ranges',
        'method': (
            'Every parameter sampled on [0.8*theta_base, 1.2*theta_base] with '
            'the same LHS construction, n=2000 and seed as task A; identical '
            'dynamic outputs (baseline scenario, 1095 days) and identical '
            'PRCC estimator.'
        ),
        'ranges': {k: {'lo': v[0], 'hi': v[1], 'baseline': BASE_VALUES[k],
                       'relative_half_width': HALF_WIDTH}
                   for k, v in ranges.items()},
        'published_ranges_for_reference': {k: {'lo': v[0], 'hi': v[1],
                                               'relative_width': v[1] - v[0]}
                                           for k, v in S.PARAM_RANGES.items()},
        'lhs_variates_paired_with_task_A': paired,
        'outputs': blocks,
        'output_distributions': {'peak_ih': five(peak_ih),
                                 'peak_day': five(peak_day),
                                 'cum_infected': five(cum)},
        'comparison_vs_task_A': comparison,
        'n_failed_solves': n_fail,
        'meta': C.make_meta(
            script='run_equalwidth.py', seed=SEED, n_samples=N_SAMPLES,
            wall_seconds=wall, per_run_seconds=per_run,
            extra={'parallel_processes': n_proc,
                   'half_width': HALF_WIDTH}
        ),
    }
    C.write_json(os.path.join(HERE, 'equalwidth_prcc.json'), payload)


if __name__ == '__main__':
    main()
