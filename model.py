"""Coupled SEIR/SEI dengue transmission model with temperature forcing.

Implements the ODE system, temperature-dependent biting rate and extrinsic
incubation rate (Li et al., 2016), seasonal temperature, imported-case pulse,
and numerical integration via scipy.integrate.solve_ivp.
"""

import numpy as np
from scipy.integrate import solve_ivp

# ── Constants ────────────────────────────────────────────────────────────────

R_GAS = 1.987  # cal / (K * mol)

# Baseline parameters
NH = 10_000
NV_RATIO = 10          # Nv / Nh
MU_H = 1.0 / (70 * 365)
SIGMA_H = 1.0 / 5.5
GAMMA_H = 1.0 / 5.0
BETA_VH = 0.70
BETA_HV = 0.60
MU_V = 1.0 / 12.0

SEED = 20260816


# ── Temperature-dependent functions ──────────────────────────────────────────

def biting_rate(T):
    """Daily biting rate b(T) for Aedes albopictus (Li et al., 2016).

    Linear in T, clamped to [0.18, 0.23] for 21–32°C.
    """
    b = 0.0043 * T + 0.0943
    return np.clip(b, 0.18, 0.23)


def eip_rate(T):
    """Extrinsic incubation rate sigma_v(T) in 1/day (Li et al., 2016).

    Returns the rate (reciprocal of EIP in days).
    """
    Tk = T + 273.15
    numerator = 24.0 * 0.00333 * (Tk / 298.0) * np.exp(
        (60513.2 / R_GAS) * (1.0 / 298.0 - 1.0 / Tk)
    )
    denominator = 1.0 + np.exp(
        (705550.0 / R_GAS) * (1.0 / 308.352 - 1.0 / Tk)
    )
    return numerator / denominator


def seasonal_temperature(t):
    """Guangzhou seasonal temperature in Celsius.

    T(t) = 22 + 7.8 * cos(2*pi*(t - 200)/365)
    Warmest day around late July (t ~ 200).
    """
    return 22.0 + 7.8 * np.cos(2.0 * np.pi * (t - 200.0) / 365.0)


def import_pulse(t, amp=0.5):
    """Seasonal imported-infectious pulse.

    Gaussian peak at t = 210 (late July), repeated annually.
    """
    result = 0.0
    for k in range(-1, 10):
        center = 210.0 + 365.0 * k
        result += np.exp(-((t - center) / 30.0) ** 2)
    return amp * result


# ── ODE system ───────────────────────────────────────────────────────────────

def rhs(t, y, params):
    """Right-hand side of the coupled SEIR/SEI ODE system.

    State vector: [Sh, Eh, Ih, Rh, Sv, Ev, Iv]
    """
    Sh, Eh, Ih, Rh, Sv, Ev, Iv = y

    T = seasonal_temperature(t) + params.get('temp_offset', 0.0)
    b = biting_rate(T)
    sv = eip_rate(T)

    nh = params['Nh']
    nv_ratio = params['nv_ratio']
    mu_h = params['mu_h']
    sigma_h = params['sigma_h']
    gamma_h = params['gamma_h']
    beta_vh = params['beta_vh']
    beta_hv = params['beta_hv']
    mu_v = params['mu_v']
    control = params.get('control_factor', 1.0)
    import_amp = params.get('import_amp', 0.5)

    nv_eq = nv_ratio * nh * control

    # Human equations
    dSh = mu_h * nh - beta_vh * b * (Iv / nh) * Sh - mu_h * Sh
    dEh = beta_vh * b * (Iv / nh) * Sh - (sigma_h + mu_h) * Eh
    dIh = sigma_h * Eh - (gamma_h + mu_h) * Ih + import_pulse(t, import_amp)
    dRh = gamma_h * Ih - mu_h * Rh

    # Vector equations
    dSv = mu_v * nv_eq - beta_hv * b * (Ih / nh) * Sv - mu_v * Sv
    dEv = beta_hv * b * (Ih / nh) * Sv - (sv + mu_v) * Ev
    dIv = sv * Ev - mu_v * Iv

    return [dSh, dEh, dIh, dRh, dSv, dEv, dIv]


def baseline_params(**overrides):
    """Return baseline parameter dict with optional overrides."""
    p = {
        'Nh': NH,
        'nv_ratio': NV_RATIO,
        'mu_h': MU_H,
        'sigma_h': SIGMA_H,
        'gamma_h': GAMMA_H,
        'beta_vh': BETA_VH,
        'beta_hv': BETA_HV,
        'mu_v': MU_V,
        'control_factor': 1.0,
        'import_amp': 0.5,
        'temp_offset': 0.0,
    }
    p.update(overrides)
    return p


def solve(params, t_span=(0, 1095), n_points=1000):
    """Integrate the ODE system and return t, y arrays.

    Default: 3 years (1095 days), 1000 output points.
    """
    y0 = [
        params['Nh'], 0.0, 0.0, 0.0,           # Sh, Eh, Ih, Rh
        params['nv_ratio'] * params['Nh'], 0.0, 0.0  # Sv, Ev, Iv
    ]
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    sol = solve_ivp(
        rhs, t_span, y0, args=(params,),
        method='RK45', t_eval=t_eval, max_step=1.0, rtol=1e-8, atol=1e-10
    )
    return sol.t, sol.y


if __name__ == '__main__':
    p = baseline_params()
    t, y = solve(p)
    peak_ih = np.max(y[2])
    peak_day = t[np.argmax(y[2])]
    print(f"Baseline peak Ih: {peak_ih:.0f} at day {peak_day:.0f}")
