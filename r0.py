"""Basic reproduction number R0 via the next-generation matrix.

Implements the NGM method (van den Driessche & Watmough, 2002) for the
coupled SEIR/SEI dengue model. Infected compartments: [Eh, Ih, Ev, Iv].
"""

import numpy as np
from model import biting_rate, eip_rate, baseline_params


def next_generation_matrix(params, T=28.0):
    """Construct the next-generation matrix K = F * V^(-1) at temperature T.

    Parameters
    ----------
    params : dict
        Model parameters (must contain beta_vh, beta_hv, mu_v, sigma_h,
        gamma_h, mu_h, nv_ratio).
    T : float
        Temperature in Celsius.

    Returns
    -------
    K : ndarray, shape (4, 4)
        Next-generation matrix.
    """
    b = biting_rate(T)
    sv = eip_rate(T)

    beta_vh = params['beta_vh']
    beta_hv = params['beta_hv']
    mu_v = params['mu_v']
    sigma_h = params['sigma_h']
    gamma_h = params['gamma_h']
    mu_h = params['mu_h']
    m = params['nv_ratio']

    # F: new-infection matrix (4x4)
    # F[0, 3] = beta_vh * b  (Iv -> new Eh)
    # F[2, 1] = beta_hv * b * m  (Ih -> new Ev)
    F = np.zeros((4, 4))
    F[0, 3] = beta_vh * b
    F[2, 1] = beta_hv * b * m

    # V: transition matrix (4x4)
    # V[0,0] = sigma_h + mu_h   (Eh outflow)
    # V[1,0] = -sigma_h          (Eh -> Ih)
    # V[1,1] = gamma_h + mu_h   (Ih outflow)
    # V[2,2] = sigma_v + mu_v   (Ev outflow)
    # V[3,2] = -sigma_v          (Ev -> Iv)
    # V[3,3] = mu_v              (Iv outflow)
    V = np.zeros((4, 4))
    V[0, 0] = sigma_h + mu_h
    V[1, 0] = -sigma_h
    V[1, 1] = gamma_h + mu_h
    V[2, 2] = sv + mu_v
    V[3, 2] = -sv
    V[3, 3] = mu_v

    V_inv = np.linalg.inv(V)
    K = F @ V_inv
    return K


def compute_r0(params, T=28.0):
    """Compute R0 as the spectral radius of the next-generation matrix.

    Parameters
    ----------
    params : dict
        Model parameters.
    T : float
        Temperature in Celsius.

    Returns
    -------
    float
        Basic reproduction number R0.
    """
    K = next_generation_matrix(params, T)
    eigenvalues = np.linalg.eigvals(K)
    return np.max(np.abs(eigenvalues))


def r0_temperature_curve(params, T_range=(15, 35), n_points=100):
    """Compute R0 across a temperature range.

    Returns
    -------
    temps : ndarray
        Temperature values.
    r0s : ndarray
        Corresponding R0 values.
    """
    temps = np.linspace(T_range[0], T_range[1], n_points)
    r0s = np.array([compute_r0(params, T) for T in temps])
    return temps, r0s


if __name__ == '__main__':
    p = baseline_params()
    for T in [15, 18, 19, 20, 22, 25, 28, 30, 32, 35]:
        r0 = compute_r0(p, T)
        print(f"R0({T:.0f} C) = {r0:.3f}")

    # Verify key values
    assert abs(compute_r0(p, 28.0) - 2.928) < 0.01, "R0(28) mismatch"
    assert abs(compute_r0(p, 20.0) - 1.10) < 0.02, "R0(20) mismatch"
    assert abs(compute_r0(p, 25.0) - 2.27) < 0.02, "R0(25) mismatch"
    assert abs(compute_r0(p, 30.0) - 3.26) < 0.02, "R0(30) mismatch"
    assert abs(compute_r0(p, 32.0) - 3.49) < 0.02, "R0(32) mismatch"
    print("\nAll R0 values verified.")
