"""Pipeline orchestrator: runs all analyses and generates 4 figures at 300 dpi.

Usage:
    python run_all.py

Output:
    Figure1_r0_temperature.png
    Figure2_prcc_sensitivity.png
    Figure3_scenario_peaks.png
    Figure4_incidence_timeseries.png
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from model import baseline_params, solve, seasonal_temperature
from r0 import compute_r0, r0_temperature_curve
from sensitivity import run as run_sensitivity
from scenarios import run_all_scenarios, SCENARIO_NAMES

OUTPUT_DIR = '.'


def figure1_r0_temperature():
    """Figure 1: R0 as a function of temperature."""
    params = baseline_params()
    temps, r0s = r0_temperature_curve(params, T_range=(15, 35), n_points=200)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(temps, r0s, 'b-', linewidth=2)
    ax.axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.7)
    ax.set_xlabel('Temperature (°C)', fontsize=13)
    ax.set_ylabel('R0', fontsize=13)
    ax.set_title('Basic reproduction number vs temperature', fontsize=14)
    ax.set_xlim(15, 35)
    ax.set_ylim(0, max(r0s) * 1.1)
    ax.grid(True, alpha=0.3)

    # Annotate key temperatures (offset above-left of the marker, white box for legibility)
    for T in [20, 25, 28, 30, 32]:
        r0 = compute_r0(params, T)
        ax.plot(T, r0, 'ro', markersize=5)
        ax.annotate(f'R0={r0:.2f}', xy=(T, r0),
                    xytext=(-13, 9), textcoords='offset points',
                    ha='center', va='bottom', fontsize=9,
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85))

    fig.tight_layout()
    path = f'{OUTPUT_DIR}/Figure1_r0_temperature.png'
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"  -> {path}")
    return path


def figure2_prcc():
    """Figure 2: PRCC sensitivity bar chart."""
    result = run_sensitivity(n_samples=2000)

    labels = result['param_labels']
    prcc_vals = result['prcc']

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ['#d62728' if v < 0 else '#1f77b4' for v in prcc_vals]
    bars = ax.barh(range(len(labels)), prcc_vals, color=colors, edgecolor='black',
                   linewidth=0.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel('PRCC', fontsize=13)
    ax.set_title('Global sensitivity analysis (LHS-PRCC, n=2000)', fontsize=14)
    ax.axvline(x=0, color='black', linewidth=0.8)
    ax.set_xlim(-1.25, 1.25)
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)

    for i, (bar, val) in enumerate(zip(bars, prcc_vals)):
        # place value label just outside the bar end, well inside the extended x-limit
        ax.text(val + (0.03 if val >= 0 else -0.03), i, f'{val:+.3f}',
                va='center', ha='left' if val >= 0 else 'right', fontsize=10)

    fig.tight_layout()
    path = f'{OUTPUT_DIR}/Figure2_prcc_sensitivity.png'
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"  -> {path}")
    return path


def figure3_scenario_peaks():
    """Figure 3: Scenario peak comparison bar chart."""
    results = run_all_scenarios()

    names = [r['name'] for r in results]
    peaks = [r['peak_ih'] for r in results]

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ['#4e79a7', '#59a14f', '#e15759', '#e15759', '#f28e2b']
    bars = ax.bar(range(len(names)), peaks, color=colors, edgecolor='black',
                  linewidth=0.5)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=11, rotation=20, ha='center')
    ax.set_ylabel('Peak infectious humans', fontsize=13)
    ax.set_title('Peak infectious prevalence by scenario', fontsize=14)
    ax.set_ylim(0, max(peaks) * 1.18)
    ax.grid(True, axis='y', alpha=0.3)

    for bar, peak in zip(bars, peaks):
        ax.text(bar.get_x() + bar.get_width() / 2, peak + (max(peaks) * 0.02),
                f'{peak:.0f}', ha='center', va='bottom', fontsize=10)

    fig.tight_layout()
    path = f'{OUTPUT_DIR}/Figure3_scenario_peaks.png'
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"  -> {path}")
    return path


def figure4_timeseries():
    """Figure 4: Infectious human time series (baseline vs vector control)."""
    from scenarios import get_scenario

    fig, ax = plt.subplots(figsize=(10, 5))

    for name, color in [('Baseline', '#4e79a7'), ('Vector control', '#59a14f')]:
        params = get_scenario(name)
        t, y = solve(params, t_span=(0, 1095), n_points=1000)
        ax.plot(t, y[2], label=name, color=color, linewidth=1.5)

    ax.set_xlabel('Time (days)', fontsize=13)
    ax.set_ylabel('Infectious humans', fontsize=13)
    ax.set_title('Dengue outbreak dynamics: baseline vs vector control', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1095)

    fig.tight_layout()
    path = f'{OUTPUT_DIR}/Figure4_incidence_timeseries.png'
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"  -> {path}")
    return path


def main():
    print("=" * 60)
    print("Dengue SEIR/SEI Framework - Full Pipeline")
    print("=" * 60)

    print("\n[1/5] Computing R0-temperature curve...")
    params = baseline_params()
    for T in [20, 25, 28, 30, 32]:
        r0 = compute_r0(params, T)
        print(f"  R0({T}°C) = {r0:.3f}")

    print("\n[2/5] Running LHS-PRCC sensitivity analysis...")
    result = run_sensitivity(n_samples=2000)
    print(f"  R0 range: {result['r0_min']:.2f} - {result['r0_max']:.2f} "
          f"(median {result['r0_median']:.2f})")
    for i in range(len(result['param_names'])):
        print(f"  {result['param_labels'][i]:<15} PRCC = {result['prcc'][i]:+.3f}")

    print("\n[3/5] Running scenario simulations...")
    from scenarios import run_all_scenarios, print_summary
    scenarios = run_all_scenarios()
    print_summary(scenarios)

    print("\n[4/5] Generating figures...")
    figure1_r0_temperature()
    figure2_prcc()
    figure3_scenario_peaks()
    figure4_timeseries()

    print("\n[5/5] Done. All figures generated at 300 dpi.")
    print("=" * 60)


if __name__ == '__main__':
    main()
