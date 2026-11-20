"""Task A: PRCC against *dynamic* outputs (baseline scenario, 3 years).

Same LHS design (n=2000, seed=20260816) and same 8 parameters as the published
static analysis, but every sample is now pushed through the ODE model instead
of the algebraic R0 at a fixed 28 degC.

Mapping of the two temperature-dependent parameters
---------------------------------------------------
The published static PRCC *replaces* b and sigma_v by constants, which is the
main limitation of that design.  Here the sampled values are applied as level
scalings of the temperature response, so the seasonal forcing stays in play:

    b_eff(T)      = b_model(T)      * (b_sample / b_model(28 degC))
    sigma_v_eff(T)= sigma_v_model(T)* ((1/eip_days_sample) / sigma_v_model(28 degC))

with b_model(28) = 0.2147 /day and sigma_v_model(28) = 0.23522 /day
(= 1 / 4.2515 days).  The sampled ranges therefore keep their published width
and interpretation ("value of b at the 28 degC reference point"), while the
dynamics remain temperature-forced.

Outputs
-------
peak_ih        maximum daily infectious humans
peak_day       day of the peak (parabolic sub-daily refinement)
cum_infected   1 - Sh(1095)/Nh
"""

import os
import sys
import time
import multiprocessing as mp

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'src')
sys.path.insert(0, SRC)
sys.path.insert(0, HERE)

import dv_common as C
import model as M
import sensitivity as S

SEED = 20260816
N_SAMPLES = 2000


def build_jobs(n_samples):
    samples = S.sample_lhs(n_samples=n_samples, seed=SEED)
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
    return samples, jobs


def main():
    t_start = time.perf_counter()

    # ── validation: scale = 1 must reproduce the published baseline ─────────
    C.set_rate_scales(1.0, 1.0)
    base = C.run_case(M.baseline_params())
    C.reset_rates()
    print(f"[check] baseline peak_ih={base['peak_ih']:.2f} "
          f"day={base['peak_day']:.2f} cum={base['cum_infected']:.4f}")

    t0 = time.perf_counter()
    C.run_case(M.baseline_params())
    per_run = time.perf_counter() - t0
    print(f"[timing] one ODE solve = {per_run:.3f} s; "
          f"{N_SAMPLES} solves ~ {per_run * N_SAMPLES / 60:.1f} min serial")

    samples, jobs = build_jobs(N_SAMPLES)

    n_proc = max(1, min(7, mp.cpu_count() - 1))
    print(f"[run] {len(jobs)} solves on {n_proc} processes ...")
    if n_proc > 1:
        with mp.Pool(n_proc) as pool:
            results = pool.map(C.worker, jobs, chunksize=16)
    else:
        results = [C.worker(j) for j in jobs]

    wall = time.perf_counter() - t_start
    n_fail = sum(0 if r['ok'] else 1 for r in results)
    print(f"[run] done in {wall:.1f} s ({n_fail} failed solves)")

    peak_ih = np.array([r['peak_ih'] for r in results])
    peak_day = np.array([r['peak_day'] for r in results])
    cum = np.array([r['cum_infected'] for r in results])
    cum_chk = np.array([r['cum_check'] for r in results])

    outputs = {
        'peak_ih': peak_ih,
        'peak_day': peak_day,
        'cum_infected': cum,
    }

    blocks = {}
    for out_name, arr in outputs.items():
        rows, dropped = C.prcc_table(samples, arr, S.PARAM_NAMES)
        blocks[out_name] = {'prcc_table': rows, 'n_dropped_nonfinite': dropped}
        print(f"\n--- PRCC vs {out_name} (n={rows[0]['n_eff']}, df={rows[0]['df']}) ---")
        print(f"{'rank':<5}{'parameter':<12}{'PRCC':>9}{'p':>12}  sig")
        for r in rows:
            sig = '***' if r['p_value'] < 0.001 else (
                '**' if r['p_value'] < 0.01 else ('*' if r['p_value'] < 0.05 else 'n.s.'))
            print(f"{r['rank']:<5}{r['param']:<12}{r['prcc']:>9.3f}"
                  f"{C.fmt_p(r['p_value']):>12}  {sig}")

    # ── descriptive statistics of the outputs ──────────────────────────────
    def five(x):
        x = x[np.isfinite(x)]
        return {
            'min': float(np.min(x)), 'p05': float(np.percentile(x, 5)),
            'median': float(np.median(x)), 'p95': float(np.percentile(x, 95)),
            'max': float(np.max(x)), 'mean': float(np.mean(x)),
        }

    payload = {
        'task': 'A - PRCC of dynamic model outputs (baseline scenario)',
        'method': (
            'LHS (sensitivity.sample_lhs, n=2000, seed=20260816) over the 8 '
            'published parameter ranges; each sample propagated through the '
            'full ODE model (model.solve, t_span=(0,1095), daily output). '
            'PRCC = Marino et al. (2008): rank-transform inputs and output, '
            'regress each on the remaining 7 ranked inputs, Pearson '
            'correlation of residuals. p-values from a two-sided t-test with '
            'df = n - 2 - (k-1) = n - 9.'
        ),
        'temperature_handling': (
            'b and sigma_v keep their temperature dependence; the sampled '
            'values act as level scalings relative to the 28 degC reference '
            'value (b_ref=%.6f /day, sigma_v_ref=%.6f /day = 1/%.4f d). '
            'This avoids the constant-parameter substitution used in the '
            'static analysis.' % (C.B_REF, C.SV_REF, C.EIP_REF_DAYS)
        ),
        'cumulative_metric': (
            "cum_infected = 1 - Sh(1095)/Nh (proportion ever infected, "
            "attack-rate analogue). Cross-check (Eh+Ih+Rh)(1095)/Nh differs "
            "only by the imported cases injected directly into Ih; mean "
            "absolute difference = %.2e." % float(np.mean(np.abs(cum - cum_chk)))
        ),
        'peak_definition': (
            'peak_ih / peak_day from a 1096-point daily grid, then refined by '
            'a 3-point parabola around the argmax (needed to break the 1-day '
            'quantisation of peak_day in the rank transform). Peak value '
            'bias of the raw grid relative to a 1000-point grid is < 0.5%.'
        ),
        'baseline_reference_run': base,
        'outputs': blocks,
        'output_distributions': {
            'peak_ih': five(peak_ih),
            'peak_day': five(peak_day),
            'cum_infected': five(cum),
            'cum_check_minus_cum_infected': {
                'mean': float(np.mean(cum_chk - cum)),
                'max_abs': float(np.max(np.abs(cum_chk - cum))),
            },
        },
        'n_failed_solves': n_fail,
        'meta': C.make_meta(
            script='run_dynamic_prcc.py', seed=SEED, n_samples=N_SAMPLES,
            wall_seconds=wall, per_run_seconds=per_run,
            extra={'parallel_processes': n_proc,
                   'b_ref_28C': C.B_REF, 'sigma_v_ref_28C': C.SV_REF,
                   'eip_ref_days': C.EIP_REF_DAYS}
        ),
    }
    C.write_json(os.path.join(HERE, 'dynamic_prcc.json'), payload)


if __name__ == '__main__':
    main()
