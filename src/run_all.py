"""Pipeline orchestrator: runs all analyses and generates 5 figures at 300 dpi.

Usage:
    python run_all.py                 # default seed 20260816
    python run_all.py --seed 42       # alternative seed for the LHS-PRCC step
    DENGUE_SEED=42 python run_all.py  # equivalent environment-variable form

Output:
    Figure1_r0_temperature.png
    Figure2_prcc_sensitivity.png
    Figure3_scenario_peaks.png
    Figure4_incidence_timeseries.png
    Figure5_prcc_designs.png
"""

import argparse
import json
import os


def _parse_args(argv=None):
    """Parse command-line options (seed only)."""
    parser = argparse.ArgumentParser(
        description="Run the full dengue modelling pipeline and write 5 figures."
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help=(
            "Pseudo-random seed for the LHS-PRCC analysis. Overrides the "
            "DENGUE_SEED environment variable; the default (20260816) "
            "reproduces the published results."
        ),
    )
    return parser.parse_args(argv)


# The seed must be fixed before ``model``/``sensitivity`` are imported, because
# those modules resolve it once at import time (see config.py).  Parsing only
# happens for ``python run_all.py``, never when this file is imported.
_ARGS = _parse_args() if __name__ == "__main__" else None
if _ARGS is not None and _ARGS.seed is not None:
    os.environ["DENGUE_SEED"] = str(_ARGS.seed)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config import SEED
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


def figure5_prcc_designs():
    """Figure 5: PRCC under the static, dynamic and equal-width designs."""
    here = os.path.dirname(os.path.abspath(__file__))
    na_dir = os.path.join(here, '..', 'new_analysis')
    files = {'static': 'static_prcc_reference.json',
             'dynamic': 'dynamic_prcc.json',
             'equal': 'equalwidth_prcc.json'}
    for name in files.values():
        path = os.path.join(na_dir, name)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path} not found; Figure 5 reads the archived PRCC tables "
                f"written by the scripts in new_analysis/"
            )
    with open(os.path.join(na_dir, files['static']), encoding='utf-8') as fh:
        static_raw = json.load(fh)
    with open(os.path.join(na_dir, files['dynamic']), encoding='utf-8') as fh:
        dynamic_raw = json.load(fh)
    with open(os.path.join(na_dir, files['equal']), encoding='utf-8') as fh:
        equal_raw = json.load(fh)

    static = dict(zip(static_raw['param_names'], static_raw['prcc']))
    dyn_rows = dynamic_raw['outputs']['peak_ih']['prcc_table']
    dynamic = {r['param']: r['prcc'] for r in dyn_rows}
    equal = {r['param']: r['prcc'] for r in equal_raw['outputs']['peak_ih']['prcc_table']}
    order = [r['param'] for r in sorted(dyn_rows, key=lambda r: r['rank'])]

    labels = {'nv_ratio': r'$N_v/N_h$', 'beta_vh': r'$\beta_{vh}$',
              'beta_hv': r'$\beta_{hv}$', 'gamma_h': r'$\gamma_h$',
              'b': r'$b$', 'mu_v': r'$\mu_v$',
              'eip_days': r'$1/\sigma_v$', 'sigma_h': r'$\sigma_h$'}

    fig, ax = plt.subplots(figsize=(9, 5.2))
    ypos = list(range(len(order)))
    for y, key in zip(ypos, order):
        vals = [static[key], dynamic[key], equal[key]]
        ax.plot([min(vals), max(vals)], [y, y], color='0.6', linewidth=0.9,
                zorder=1)
    ax.axvline(0.0, color='0.4', linewidth=0.8, zorder=0)
    ax.plot([static[k] for k in order], [y + 0.20 for y in ypos], 'o',
            markersize=6.5, markerfacecolor='white', markeredgecolor='#4e79a7',
            markeredgewidth=1.2, linestyle='none', label='Static, $R_0$ at 28 \u00b0C',
            zorder=5)
    ax.plot([dynamic[k] for k in order], ypos, 's', markersize=6.4,
            color='#4e79a7', linestyle='none', label='Dynamic, peak', zorder=4)
    ax.plot([equal[k] for k in order], [y - 0.20 for y in ypos], '^',
            markersize=6.4, color='#E15759', linestyle='none',
            label='Equal-width, peak', zorder=3)

    ax.set_yticks(ypos)
    ax.set_yticklabels([labels.get(k, k) for k in order])
    ax.invert_yaxis()
    ax.set_xlim(-1.06, 1.06)
    ax.set_xlabel('Partial rank correlation coefficient', fontsize=13)
    ax.set_title('Partial rank correlation coefficients under three designs',
                 fontsize=14)
    ax.legend(fontsize=11, frameon=False, loc='upper left')
    ax.grid(True, alpha=0.3, axis='x')

    fig.tight_layout()
    path = f'{OUTPUT_DIR}/Figure5_prcc_designs.png'
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"  -> {path}")
    return path


def main():
    print("=" * 60)
    print("Dengue SEIR/SEI Framework - Full Pipeline")
    print("=" * 60)
    print(f"Random seed: {SEED}")

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
    figure5_prcc_designs()

    print("\n[5/5] Done. All figures generated at 300 dpi.")
    print("=" * 60)


if __name__ == '__main__':
    main()
