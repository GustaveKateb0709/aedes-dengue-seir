"""Numerical verification suite for the dengue SEIR/SEI framework.

Every assertion below pins one number that appears in the manuscript, so a
regression in the code (or a silent change of parameters) fails loudly.

Run it either way -- no test runner is required:

    python tests/test_verification.py
    pytest tests/test_verification.py        # if pytest is installed

The file deliberately never imports pytest, so the standalone path works in a
bare environment (numpy/scipy only).
"""

import os
import subprocess
import sys

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
PKG_ROOT = os.path.dirname(TESTS_DIR)
SRC = os.path.join(PKG_ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import numpy as np  # noqa: E402

import config  # noqa: E402
import model  # noqa: E402
from model import biting_rate, eip_rate  # noqa: E402
from model import baseline_params, solve  # noqa: E402
from r0 import compute_r0  # noqa: E402
from scenarios import run_all_scenarios  # noqa: E402

PARAMS = baseline_params()

# Clamp bounds applied inside model.biting_rate (Li et al., 2016).
BITING_RATE_MIN = 0.18
BITING_RATE_MAX = 0.23

# Scenario integration is the slow part (~seconds); solve it once.
_SCENARIOS = None


def scenarios():
    """Return (and cache) the five scenario results."""
    global _SCENARIOS
    if _SCENARIOS is None:
        _SCENARIOS = run_all_scenarios()
    return _SCENARIOS


# ── R0 across temperature ───────────────────────────────────────────────────

def test_r0_at_28c():
    """R0(28 °C) = 2.928 (manuscript headline value)."""
    assert abs(compute_r0(PARAMS, 28.0) - 2.928) < 0.01


def test_r0_at_25c():
    """R0(25 °C) = 2.270."""
    assert abs(compute_r0(PARAMS, 25.0) - 2.270) < 0.01


def test_r0_at_30c():
    """R0(30 °C) = 3.262."""
    assert abs(compute_r0(PARAMS, 30.0) - 3.262) < 0.01


def test_r0_at_32c():
    """R0(32 °C) = 3.492."""
    assert abs(compute_r0(PARAMS, 32.0) - 3.492) < 0.01


def test_r0_below_one_at_19c():
    """R0(19 °C) = 0.940 < 1, i.e. the transmission threshold is near 19-20 °C."""
    r0 = compute_r0(PARAMS, 19.0)
    assert abs(r0 - 0.940) < 0.01
    assert r0 < 1.0


def test_r0_monotonic_over_operational_range():
    """R0 increases monotonically with temperature over 17-33 °C."""
    temps = np.linspace(17.0, 33.0, 17)
    r0s = np.array([compute_r0(PARAMS, T) for T in temps])
    assert np.all(np.diff(r0s) > 0)


# ── Temperature-dependent entomological parameters ──────────────────────────

def test_eip_days_at_28c():
    """Extrinsic incubation period 1/sigma_v(28 °C) = 4.251 days."""
    assert abs(1.0 / eip_rate(28.0) - 4.251) < 1e-3


def test_eip_rate_at_28c():
    """sigma_v(28 °C) = 0.2352 per day."""
    assert abs(eip_rate(28.0) - 0.2352) < 1e-3


def test_biting_rate_lower_clip():
    """b(15 °C) is clamped to the lower bound 0.180 (unclipped value is 0.159)."""
    assert abs(biting_rate(15.0) - 0.180) < 1e-12
    assert abs(0.0043 * 15.0 + 0.0943 - biting_rate(15.0)) > 1e-3
    assert BITING_RATE_MIN <= biting_rate(15.0) <= BITING_RATE_MAX


def test_biting_rate_upper_clip():
    """b(35 °C) is clamped to the upper bound 0.230 (unclipped value is 0.2448)."""
    assert abs(biting_rate(35.0) - 0.230) < 1e-12


def test_biting_rate_interior_value():
    """b(28 °C) = 0.2147 lies inside the clamp and follows the linear relation."""
    assert abs(biting_rate(28.0) - 0.2147) < 1e-4


# ── Scenario simulations ────────────────────────────────────────────────────

def test_scenario_peak_baseline():
    """Baseline peak infectious humans = 1650 (tolerance from scenarios.py)."""
    assert abs(scenarios()[0]["peak_ih"] - 1650) < 15


def test_scenario_peak_vector_control():
    """Vector control peak = 120, i.e. a ~93% reduction."""
    peak = scenarios()[1]["peak_ih"]
    assert abs(peak - 120) < 10
    assert peak < scenarios()[0]["peak_ih"]


def test_scenario_peak_warming_plus1():
    """Warming +1 °C peak = 1796."""
    assert abs(scenarios()[2]["peak_ih"] - 1796) < 15


def test_scenario_peak_warming_plus2():
    """Warming +2 °C peak = 1904."""
    assert abs(scenarios()[3]["peak_ih"] - 1904) < 15


def test_scenario_peak_import_4x():
    """Import 4x peak = 1721."""
    assert abs(scenarios()[4]["peak_ih"] - 1721) < 15


def test_baseline_peak_day():
    """Baseline epidemic peaks at day 233.5 (+/- 2 days)."""
    assert abs(scenarios()[0]["peak_day"] - 233.5) <= 2.0


def test_baseline_peak_from_solver():
    """Independent re-solve of the baseline reproduces the cached peak."""
    t, y = solve(baseline_params(), t_span=(0, 1095), n_points=1000)
    assert abs(float(np.max(y[2])) - scenarios()[0]["peak_ih"]) < 1e-6
    assert abs(float(t[int(np.argmax(y[2]))]) - scenarios()[0]["peak_day"]) < 1e-9


# ── Configuration / seed handling ───────────────────────────────────────────

def test_default_seed():
    """The default seed is 20260816 and is exposed through model.SEED."""
    assert config.DEFAULT_SEED == 20260816
    assert model.SEED == 20260816
    assert config.get_seed() == 20260816


def test_seed_environment_override():
    """DENGUE_SEED is honoured in a fresh interpreter."""
    env = dict(os.environ)
    env["DENGUE_SEED"] = "42"
    out = subprocess.run(
        [sys.executable, "-c", "import model; print(model.SEED)"],
        cwd=SRC,
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    assert out.stdout.strip() == "42"


def test_seed_changes_sensitivity_draws():
    """Different seeds give different LHS draws but leave the model outputs alone."""
    import sensitivity

    a = sensitivity.sample_lhs(n_samples=50, seed=20260816)
    b = sensitivity.sample_lhs(n_samples=50, seed=42)
    assert not np.allclose(a, b)
    # Deterministic quantities must not depend on the seed.
    assert abs(compute_r0(PARAMS, 28.0) - 2.928) < 0.01


# ── Standalone runner ───────────────────────────────────────────────────────

def main():
    """Run every test_* function in definition order and print a summary."""
    tests = [
        (name, obj)
        for name, obj in list(globals().items())
        if name.startswith("test_") and callable(obj)
    ]
    # globals() preserves definition order in CPython; sort only as a fallback.
    failures = []
    print("=" * 72)
    print("Dengue SEIR/SEI framework - numerical verification")
    print(f"Python {sys.version.split()[0]} | numpy {np.__version__} | "
          f"src = {SRC}")
    print("=" * 72)
    for name, func in tests:
        try:
            func()
        except AssertionError as exc:
            failures.append((name, f"AssertionError: {exc}"))
            print(f"FAIL  {name}")
        except Exception as exc:  # noqa: BLE001 - report, do not abort
            failures.append((name, f"{type(exc).__name__}: {exc}"))
            print(f"ERROR {name}: {type(exc).__name__}: {exc}")
        else:
            print(f"ok    {name}")

    print("-" * 72)
    print(f"{len(tests) - len(failures)} passed, {len(failures)} failed "
          f"of {len(tests)} tests")
    if failures:
        print("\nFailed tests:")
        for name, reason in failures:
            print(f"  - {name}: {reason}")
        print("\nRESULT: FAILURE")
        return 1
    print("RESULT: ALL TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
