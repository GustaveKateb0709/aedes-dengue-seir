"""Scenario definitions for dengue transmission modelling.

Five scenarios:
  1. Baseline: seasonal temperature + import pulse (amp=0.5)
  2. Vector control: mosquito abundance reduced to 30% of baseline
  3. Warming +1°C: uniform +1°C temperature offset
  4. Warming +2°C: uniform +2°C temperature offset
  5. Imported cases 4x: import pulse amplitude x4
"""

from model import baseline_params, solve
import numpy as np


SCENARIO_NAMES = [
    'Baseline',
    'Vector control',
    'Warming +1°C',
    'Warming +2°C',
    'Import 4x',
]


def get_scenario(name):
    """Return parameter dict for the named scenario."""
    if name == 'Baseline':
        return baseline_params()
    elif name == 'Vector control':
        return baseline_params(control_factor=0.30)
    elif name == 'Warming +1°C':
        return baseline_params(temp_offset=1.0)
    elif name == 'Warming +2°C':
        return baseline_params(temp_offset=2.0)
    elif name == 'Import 4x':
        return baseline_params(import_amp=2.0)
    else:
        raise ValueError(f"Unknown scenario: {name}")


def run_all_scenarios(t_span=(0, 1095), n_points=1000):
    """Run all 5 scenarios and return results.

    Returns
    -------
    results : list of dict
        Each dict has keys: name, t, y, peak_ih, peak_day.
    """
    results = []
    for name in SCENARIO_NAMES:
        params = get_scenario(name)
        t, y = solve(params, t_span=t_span, n_points=n_points)
        ih = y[2]
        peak_idx = np.argmax(ih)
        results.append({
            'name': name,
            't': t,
            'y': y,
            'peak_ih': float(ih[peak_idx]),
            'peak_day': float(t[peak_idx]),
        })
    return results


def print_summary(results):
    """Print scenario summary table."""
    baseline_peak = results[0]['peak_ih']
    print(f"{'Scenario':<20} {'Peak Ih':>10} {'Day':>8} {'% change':>10}")
    print("-" * 50)
    for r in results:
        pct = (r['peak_ih'] - baseline_peak) / baseline_peak * 100
        print(f"{r['name']:<20} {r['peak_ih']:>10.0f} {r['peak_day']:>8.0f} {pct:>+9.1f}%")


if __name__ == '__main__':
    results = run_all_scenarios()
    print_summary(results)

    # Verify key values (within tolerance for cross-platform reproducibility)
    assert abs(results[0]['peak_ih'] - 1650) < 15, f"Baseline peak mismatch: {results[0]['peak_ih']}"
    assert abs(results[1]['peak_ih'] - 120) < 10, f"Vector control peak mismatch: {results[1]['peak_ih']}"
    assert abs(results[2]['peak_ih'] - 1796) < 15, f"+1°C peak mismatch: {results[2]['peak_ih']}"
    assert abs(results[3]['peak_ih'] - 1904) < 15, f"+2°C peak mismatch: {results[3]['peak_ih']}"
    assert abs(results[4]['peak_ih'] - 1721) < 15, f"Import 4x peak mismatch: {results[4]['peak_ih']}"
    print("\nAll scenario peaks verified.")
