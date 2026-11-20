"""Shared helpers for dynamic-output sensitivity analyses.

All scripts in ``new_analysis/`` import this module. Nothing in ``../src`` is
modified: the temperature-dependent rates are replaced at run time by
multiplicative scalings of the *original* ``model.biting_rate`` /
``model.eip_rate``, so the shape of the temperature response is preserved and
only its level is resampled.  This is the key difference from the static PRCC
in ``sensitivity.py``, which replaces b and sigma_v by temperature-independent
constants and therefore bypasses the temperature forcing entirely.
"""

import json
import platform
import sys
import time
from datetime import datetime, timezone, timedelta

import numpy as np
from scipy import stats

import model as M

# ── Reference values (28 °C, the reference temperature of the static GSA) ────
B_REF = float(M.biting_rate(28.0))          # 0.2147 /day
SV_REF = float(M.eip_rate(28.0))            # 0.23522 /day
EIP_REF_DAYS = 1.0 / SV_REF                 # 4.2515 days
NH = M.NH

# Output grid: daily resolution over the 3-year baseline horizon.
N_POINTS = 1096

_ORIG_B = M.biting_rate
_ORIG_SV = M.eip_rate


def set_rate_scales(b_scale=1.0, sv_scale=1.0):
    """Apply multiplicative scalings to the temperature-dependent rates.

    b_eff(T)  = b_model(T)  * b_scale      (b_scale = b_sample / b(28 °C))
    sv_eff(T) = sv_model(T) * sv_scale     (sv_scale = (1/eip_sample) / sv(28 °C))

    The temperature *shape* is untouched; only the level is resampled.
    """
    if b_scale == 1.0 and sv_scale == 1.0:
        M.biting_rate = _ORIG_B
        M.eip_rate = _ORIG_SV
    else:
        M.biting_rate = (lambda T, s=b_scale: _ORIG_B(T) * s)
        M.eip_rate = (lambda T, s=sv_scale: _ORIG_SV(T) * s)


def reset_rates():
    M.biting_rate = _ORIG_B
    M.eip_rate = _ORIG_SV


# ── Simulation + metrics ────────────────────────────────────────────────────

def simulate(params, t_span=(0, 1095), n_points=N_POINTS):
    """Solve the ODE system and return t, y."""
    return M.solve(params, t_span=t_span, n_points=n_points)


def metrics_from_solution(t, y):
    """Extract the three dynamic outputs from a solved trajectory.

    peak_ih      : maximum daily infectious human count (parabolically
                   refined peak value, see below)
    peak_day     : day of that maximum (sub-daily, parabolic refinement of the
                   three grid points around the argmax; needed because the
                   daily grid alone quantises the peak day to 1 day and
                   produces massive ties in the rank transform)
    cum_infected : cumulative infection proportion = 1 - Sh(end)/Nh
    cum_check    : (Eh + Ih + Rh)(end)/Nh, cross-check that differs only by the
                   imported cases injected directly into Ih
    """
    ih = y[2]
    i = int(np.argmax(ih))
    peak_ih = float(ih[i])
    peak_day = float(t[i])

    if 0 < i < len(ih) - 1:
        y0, y1, y2 = float(ih[i - 1]), float(ih[i]), float(ih[i + 1])
        a = 0.5 * (y0 + y2) - y1
        bcoef = 0.5 * (y2 - y0)
        if a != 0.0:
            d = float(np.clip(-bcoef / (2.0 * a), -1.0, 1.0))
            dt = float(t[1] - t[0])
            peak_day = float(t[i] + d * dt)
            peak_ih = float(y1 + (y2 - y0) * d / 4.0)

    sh_end = float(y[0][-1])
    eh_end, ih_end, rh_end = float(y[1][-1]), float(y[2][-1]), float(y[3][-1])
    return {
        'peak_ih': peak_ih,
        'peak_day': peak_day,
        'cum_infected': 1.0 - sh_end / NH,
        'cum_check': (eh_end + ih_end + rh_end) / NH,
    }


def run_case(params, t_span=(0, 1095), n_points=N_POINTS):
    t, y = simulate(params, t_span=t_span, n_points=n_points)
    return metrics_from_solution(t, y)


def worker_constant(job):
    """Worker for the constant-parameter control: b and sigma_v do not depend
    on temperature at all (b_eff = b_sample, sigma_v_eff = 1/eip_sample)."""
    b_const, sv_const, overrides = job
    M.biting_rate = (lambda T, v=b_const: v)
    M.eip_rate = (lambda T, v=sv_const: v)
    p = M.baseline_params(**overrides)
    try:
        out = run_case(p)
        out['ok'] = True
    except Exception as exc:                              # pragma: no cover
        out = {'peak_ih': np.nan, 'peak_day': np.nan,
               'cum_infected': np.nan, 'cum_check': np.nan,
               'ok': False, 'error': repr(exc)}
    reset_rates()
    return out


def worker(job):
    """Multiprocessing worker: one (b_scale, sv_scale, overrides) case."""
    b_scale, sv_scale, overrides = job
    set_rate_scales(b_scale, sv_scale)
    p = M.baseline_params(**overrides)
    try:
        out = run_case(p)
        out['ok'] = True
    except Exception as exc:                              # pragma: no cover
        out = {'peak_ih': np.nan, 'peak_day': np.nan,
               'cum_infected': np.nan, 'cum_check': np.nan,
               'ok': False, 'error': repr(exc)}
    reset_rates()
    return out


# ── R0 with a custom b(T) ───────────────────────────────────────────────────

def r0_at(params, T, b_fn=None, sv_fn=None):
    """R0 = spectral radius of the NGM, written as sqrt(K02 * K20).

    Verified equal to ``r0.compute_r0`` to 1e-6 for T in 15..35 °C.
    """
    b = (M.biting_rate if b_fn is None else b_fn)(T)
    sv = (M.eip_rate if sv_fn is None else sv_fn)(T)
    m = params['nv_ratio']
    mu_v = params['mu_v']
    sigma_h = params['sigma_h']
    gamma_h = params['gamma_h']
    mu_h = params['mu_h']
    K02 = params['beta_vh'] * b * sv / ((sv + mu_v) * mu_v)
    K20 = params['beta_hv'] * b * m * sigma_h / ((sigma_h + mu_h) * (gamma_h + mu_h))
    return float(np.sqrt(max(K02 * K20, 0.0)))


def critical_temperature(params, b_fn=None, lo=10.0, hi=35.0, tol=1e-10):
    """Bisect R0(T) = 1 on [lo, hi]; return None if no sign change."""
    f_lo = r0_at(params, lo, b_fn) - 1.0
    f_hi = r0_at(params, hi, b_fn) - 1.0
    if f_lo * f_hi > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f_mid = r0_at(params, mid, b_fn) - 1.0
        if f_lo * f_mid <= 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
        if hi - lo < tol:
            break
    return float(0.5 * (lo + hi))


# ── LHS (mirrors sensitivity.sample_lhs) ────────────────────────────────────

def sample_lhs_ranges(names, ranges, n_samples=2000, seed=20260816):
    """Same LHS construction as ``sensitivity.sample_lhs``, arbitrary ranges.

    With the same seed the uniform draws u are identical, so two runs that
    differ only in the ranges are paired sample-for-sample.
    """
    rng = np.random.default_rng(seed)
    samples = np.zeros((n_samples, len(names)))
    for j, name in enumerate(names):
        lo, hi = ranges[name]
        perm = rng.permutation(n_samples)
        u = (perm + rng.uniform(0, 1, n_samples)) / n_samples
        samples[:, j] = lo + u * (hi - lo)
    return samples


# ── PRCC (Marino et al., 2008) ──────────────────────────────────────────────

def prcc_table(samples, output, param_names):
    """Rank-transform -> residualise on the other 7 inputs -> Pearson r.

    Delegates to ``sensitivity.prcc`` so the estimator is byte-identical to
    the one used for the published static analysis.
    """
    from sensitivity import prcc

    samples = np.asarray(samples, dtype=float)
    output = np.asarray(output, dtype=float)

    keep = np.isfinite(output) & np.all(np.isfinite(samples), axis=1)
    n_dropped = int((~keep).sum())
    r_vals, p_vals = prcc(samples[keep], output[keep])

    n_eff = int(keep.sum())
    # PRCC test: t = r * sqrt((n - 2 - (k-1)) / (1 - r^2)), df = n - 2 - (k-1)
    # (Marino et al., 2008).  The p-value returned by scipy's pearsonr uses
    # df = n - 2 and is anti-conservative here; both are reported.
    df = n_eff - 2 - (samples.shape[1] - 1)
    rows = []
    for j, name in enumerate(param_names):
        r = float(r_vals[j])
        t_stat = r * np.sqrt(df / max(1.0 - r * r, 1e-300))
        p_adj = float(2.0 * stats.t.sf(abs(t_stat), df))
        rows.append({
            'param': name,
            'prcc': r,
            'p_value': p_adj,
            'p_value_unadjusted': float(p_vals[j]),
            'abs_prcc': abs(r),
            't_stat': float(t_stat),
            'df': int(df),
            'n_eff': n_eff,
        })
    rows.sort(key=lambda d: -d['abs_prcc'])
    for rank, row in enumerate(rows, start=1):
        row['rank'] = rank
    return rows, n_dropped


# ── metadata ────────────────────────────────────────────────────────────────

def src_fingerprint():
    """SHA-256 of the unconverted source files, for provenance.

    Another teammate may touch ../src while these analyses run, so the hash
    is recorded next to every result together with a numeric check that the
    baseline is unchanged.
    """
    import hashlib
    import os as _os
    src = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'src')
    out = {}
    for name in ('model.py', 'r0.py', 'sensitivity.py', 'scenarios.py'):
        path = _os.path.join(src, name)
        if _os.path.exists(path):
            with open(path, 'rb') as fh:
                out[name] = hashlib.sha256(fh.read()).hexdigest()
    return out


def make_meta(script, seed, n_samples, wall_seconds, per_run_seconds,
              n_points=N_POINTS, extra=None):
    import numpy
    import scipy
    meta = {
        'script': script,
        'random_seed': int(seed),
        'n_samples': int(n_samples),
        'run_date_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'run_date_local': datetime.now().astimezone().isoformat(timespec='seconds'),
        'wall_seconds_total': round(float(wall_seconds), 3),
        'seconds_per_ode_solve': round(float(per_run_seconds), 5),
        'n_ode_solves': int(n_samples),
        't_span_days': [0, 1095],
        'n_output_points': int(n_points),
        'output_resolution_days': round(1095.0 / (n_points - 1), 4),
        'python': sys.version.split()[0],
        'numpy': numpy.__version__,
        'scipy': scipy.__version__,
        'platform': platform.platform(),
        'host_cpus': __import__('multiprocessing').cpu_count(),
        'model_source': '../src/model.py (unmodified)',
        'src_sha256': src_fingerprint(),
    }
    if extra:
        meta.update(extra)
    return meta


def write_json(path, payload):
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    print(f"[written] {path}")


def fmt_p(p):
    if p < 1e-300:
        return '<1e-300'
    return f"{p:.3e}"
