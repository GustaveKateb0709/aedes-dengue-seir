"""Global sensitivity analysis using LHS + PRCC.

Latin hypercube sampling (n=2000) across 8 parameters, with partial rank
correlation coefficient (PRCC) computed against R0 (Marino et al., 2008).
"""

import numpy as np
from scipy import stats

from model import baseline_params, SEED
from r0 import compute_r0


# ── Parameter definitions for sensitivity analysis ──────────────────────────

PARAM_NAMES = [
    'nv_ratio',    # Vector-to-human ratio
    'beta_hv',     # Human-to-vector transmission probability
    'beta_vh',     # Vector-to-human transmission probability
    'mu_v',        # Adult mosquito mortality (1/day)
    'gamma_h',     # Human recovery rate (1/day)
    'b',           # Daily biting rate (at reference T)
    'eip_days',    # Extrinsic incubation period (days)
    'sigma_h',     # Human incubation rate (1/day)
]

PARAM_LABELS = {
    'nv_ratio': 'Nv/Nh',
    'beta_hv': 'beta_hv',
    'beta_vh': 'beta_vh',
    'mu_v': 'mu_v',
    'gamma_h': 'gamma_h',
    'b': 'b',
    'eip_days': '1/sigma_v',
    'sigma_h': 'sigma_h',
}

PARAM_RANGES = {
    'nv_ratio': (1.0, 20.0),
    'beta_hv': (0.30, 0.70),
    'beta_vh': (0.40, 0.95),
    'mu_v': (1.0 / 14, 1.0 / 8),
    'gamma_h': (1.0 / 7, 1.0 / 4),
    'b': (0.18, 0.23),
    'eip_days': (4.0, 8.0),
    'sigma_h': (1.0 / 7, 1.0 / 4),
}


def sample_lhs(n_samples=2000, seed=SEED):
    """Latin hypercube sampling across the 8 parameters.

    Returns
    -------
    samples : ndarray, shape (n_samples, 7)
        Parameter values in the same order as PARAM_NAMES.
    """
    rng = np.random.default_rng(seed)
    n_params = len(PARAM_NAMES)
    samples = np.zeros((n_samples, n_params))

    for j, name in enumerate(PARAM_NAMES):
        lo, hi = PARAM_RANGES[name]
        # LHS: divide [0,1] into n equal intervals, sample one point per interval
        perm = rng.permutation(n_samples)
        u = (perm + rng.uniform(0, 1, n_samples)) / n_samples
        samples[:, j] = lo + u * (hi - lo)

    return samples


def samples_to_params(samples, T=28.0):
    """Convert LHS samples to parameter dicts and compute R0 for each.

    The 'b' and 'eip_days' parameters override the temperature-dependent
    functions at the reference temperature T.
    """
    from model import baseline_params

    n = samples.shape[0]
    r0_values = np.zeros(n)

    # Note: 'b' and 'sigma_v' (via eip_days) are overridden per-sample below,
    # so we do not precompute temperature-dependent rates here.

    for i in range(n):
        p = baseline_params()
        p['nv_ratio'] = samples[i, 0]
        p['beta_hv'] = samples[i, 1]
        p['beta_vh'] = samples[i, 2]
        p['mu_v'] = samples[i, 3]
        p['gamma_h'] = samples[i, 4]

        # For 'b' and 'eip_days', we override the temperature-dependent values
        b_override = samples[i, 5]
        eip_d = samples[i, 6]
        sv_override = 1.0 / eip_d
        sigma_h_override = samples[i, 7]

        # Compute R0 manually with overridden b, sigma_v, and sigma_h
        m = p['nv_ratio']
        beta_vh = p['beta_vh']
        beta_hv = p['beta_hv']
        mu_v = p['mu_v']
        sigma_h = sigma_h_override
        gamma_h = p['gamma_h']
        mu_h = p['mu_h']

        K02 = beta_vh * b_override * sv_override / ((sv_override + mu_v) * mu_v)
        K20 = beta_hv * b_override * m * sigma_h / ((sigma_h + mu_h) * (gamma_h + mu_h))
        r0_values[i] = np.sqrt(K02 * K20)

    return r0_values


def prcc(samples, output):
    """Compute PRCC for each parameter against the output.

    Parameters
    ----------
    samples : ndarray, shape (n, k)
        Input parameter samples.
    output : ndarray, shape (n,)
        Output values (e.g., R0).

    Returns
    -------
    prcc_values : ndarray, shape (k,)
        PRCC values for each parameter.
    p_values : ndarray, shape (k,)
        Two-sided p-values.
    """
    n, k = samples.shape
    # Rank-transform
    rank_samples = np.zeros_like(samples)
    for j in range(k):
        rank_samples[:, j] = stats.rankdata(samples[:, j])
    rank_output = stats.rankdata(output)

    prcc_values = np.zeros(k)
    p_values = np.zeros(k)

    for j in range(k):
        # Regress rank_samples[:, j] on all other ranked inputs
        mask = np.ones(k, dtype=bool)
        mask[j] = False
        X_others = rank_samples[:, mask]
        X_j = rank_samples[:, j]

        # Add intercept
        X_others_aug = np.column_stack([np.ones(n), X_others])

        # Residualize X_j
        beta_j = np.linalg.lstsq(X_others_aug, X_j, rcond=None)[0]
        resid_j = X_j - X_others_aug @ beta_j

        # Residualize rank_output
        beta_y = np.linalg.lstsq(X_others_aug, rank_output, rcond=None)[0]
        resid_y = rank_output - X_others_aug @ beta_y

        # Pearson correlation of residuals
        r, p = stats.pearsonr(resid_j, resid_y)
        prcc_values[j] = r
        p_values[j] = p

    return prcc_values, p_values


def run(n_samples=2000, T=28.0, seed=SEED):
    """Run full LHS + PRCC sensitivity analysis.

    Returns
    -------
    result : dict
        Dictionary with parameter names, PRCC values, p-values, R0 range.
    """
    samples = sample_lhs(n_samples, seed)
    r0_values = samples_to_params(samples, T)
    prcc_vals, p_vals = prcc(samples, r0_values)

    # Sort by absolute PRCC
    order = np.argsort(-np.abs(prcc_vals))

    result = {
        'param_names': [PARAM_NAMES[i] for i in order],
        'param_labels': [PARAM_LABELS[PARAM_NAMES[i]] for i in order],
        'prcc': prcc_vals[order],
        'p_value': p_vals[order],
        'r0_min': float(np.min(r0_values)),
        'r0_max': float(np.max(r0_values)),
        'r0_median': float(np.median(r0_values)),
        'n_samples': n_samples,
        'T': T,
    }
    return result


if __name__ == '__main__':
    result = run()
    print(f"LHS-PRCC sensitivity analysis (n={result['n_samples']}, T={result['T']}°C)")
    print(f"R0 range: {result['r0_min']:.2f} - {result['r0_max']:.2f} (median {result['r0_median']:.2f})")
    print()
    print(f"{'Parameter':<15} {'PRCC':>8} {'p-value':>12}")
    print("-" * 37)
    for i in range(len(result['param_names'])):
        sig = '***' if result['p_value'][i] < 0.001 else ''
        print(f"{result['param_labels'][i]:<15} {result['prcc'][i]:>8.3f} {result['p_value'][i]:>10.2e} {sig}")
