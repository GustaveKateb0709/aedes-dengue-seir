"""Publication figure builder for the dengue transmission framework.

Standards: publication preset (style_context "default" + Okabe-Ito),
no in-figure titles, U+2212 minus sign, TIF LZW RGB + PDF at 300 dpi.
Data provenance: model.py / r0.py / sensitivity.py / scenarios.py, seed 20260816,
plus the archived PRCC tables in ../new_analysis for Figure 5.

Portability: this script imports only local modules from this directory plus
matplotlib/numpy/scipy.  No private paths or external style assets are needed:

    python make_figures.py                 # auto-detect the output directory
    python make_figures.py --outdir figs   # explicit output directory

Auto-detection: if this directory sits inside the manuscript package
(``.../04_Code_and_Data/src``), figures are written to ``.../02_Figures``;
otherwise they go to ``./figures`` next to this script.  The environment
variable ``DENGUE_FIG_DIR`` overrides both.
"""
import argparse
import json
import os
import sys

SRC = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SRC)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style_presets import style_context
from figure_export import export_figure

from model import baseline_params
from r0 import compute_r0, r0_temperature_curve
from sensitivity import run as run_sensitivity
from scenarios import run_all_scenarios, get_scenario
from model import solve


def resolve_out_dir(cli_outdir=None):
    """Return the directory the figures are written to."""
    if cli_outdir:
        return os.path.abspath(cli_outdir)
    env = os.environ.get("DENGUE_FIG_DIR")
    if env:
        return os.path.abspath(env)
    if os.path.basename(os.path.dirname(SRC)) == "04_Code_and_Data":
        # Standard manuscript layout: 04_Code_and_Data/src -> 02_Figures
        return os.path.abspath(os.path.join(SRC, "..", "..", "02_Figures"))
    # Extracted archive layout: <root>/src -> <root>/figures
    return os.path.join(os.path.dirname(SRC), "figures")


OUT = resolve_out_dir()

VERIFY = {}

PRCC_LABELS = {"nv_ratio": "$N_v/N_h$", "beta_vh": r"$\beta_{vh}$",
               "beta_hv": r"$\beta_{hv}$", "gamma_h": r"$\gamma_h$",
               "b": "$b$", "mu_v": r"$\mu_v$",
               "eip_days": r"$1/\sigma_v$", "sigma_h": r"$\sigma_h$"}


def load_prcc_designs():
    """Read the archived PRCC tables for the static, dynamic and equal-width designs.

    Returns (static, dynamic, equal, order) where each design maps a parameter
    key to its PRCC and ``order`` lists the keys by dynamic-design rank.
    """
    na_dir = os.path.abspath(os.path.join(SRC, "..", "new_analysis"))
    needed = ["static_prcc_reference.json", "dynamic_prcc.json",
              "equalwidth_prcc.json"]
    for name in needed:
        path = os.path.join(na_dir, name)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path} not found; Figure 5 needs the archived PRCC tables "
                f"written by the scripts in new_analysis/"
            )
    with open(os.path.join(na_dir, "static_prcc_reference.json"), encoding="utf-8") as fh:
        static_raw = json.load(fh)
    with open(os.path.join(na_dir, "dynamic_prcc.json"), encoding="utf-8") as fh:
        dynamic_raw = json.load(fh)
    with open(os.path.join(na_dir, "equalwidth_prcc.json"), encoding="utf-8") as fh:
        equal_raw = json.load(fh)
    static = dict(zip(static_raw["param_names"], static_raw["prcc"]))
    dyn_rows = dynamic_raw["outputs"]["peak_ih"]["prcc_table"]
    eq_rows = equal_raw["outputs"]["peak_ih"]["prcc_table"]
    dynamic = {r["param"]: r["prcc"] for r in dyn_rows}
    equal = {r["param"]: r["prcc"] for r in eq_rows}
    order = [r["param"] for r in sorted(dyn_rows, key=lambda r: r["rank"])]
    return static, dynamic, equal, order


def fig5_prcc_designs(out_dir=None):
    """Figure 5: PRCC under the static, dynamic and equal-width designs."""
    out = resolve_out_dir(out_dir)
    os.makedirs(out, exist_ok=True)
    static, dynamic, equal, order = load_prcc_designs()
    ypos = list(range(len(order)))

    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, ax = plt.subplots(figsize=(6.9, 4.2), layout="constrained")
        # Vertical dodge inside each row (triangle above, square on the row
        # line, circle below) so that coincident estimates stay visible.
        DY = {"tri": -0.20, "sq": 0.0, "circ": 0.20}
        for y, key in zip(ypos, order):
            vals = [static[key], dynamic[key], equal[key]]
            ax.plot([min(vals), max(vals)], [y, y], color="0.60",
                    linewidth=0.9, zorder=1, solid_capstyle="round")
        ax.axvline(0.0, color="#555555", linewidth=0.8, zorder=0)
        ax.plot([static[k] for k in order], [y + DY["circ"] for y in ypos],
                linestyle="none", marker="o", markersize=6.3,
                markerfacecolor="white", markeredgecolor="#0072B2",
                markeredgewidth=1.1,
                label="Static, $R_0$ at 28 \u00b0C", zorder=5)
        ax.plot([dynamic[k] for k in order], [y + DY["sq"] for y in ypos],
                linestyle="none", marker="s", markersize=6.2,
                color="#0072B2", label="Dynamic, peak", zorder=4)
        ax.plot([equal[k] for k in order], [y + DY["tri"] for y in ypos],
                linestyle="none", marker="^", markersize=6.2,
                color="#D55E00", label="Equal-width, peak", zorder=3)
        ax.set_yticks(ypos)
        ax.set_yticklabels([PRCC_LABELS.get(k, k) for k in order])
        ax.invert_yaxis()
        ax.set_xlim(-1.06, 1.06)
        ax.set_xticks([-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels([f"{v:.2f}".replace("-", "\u2212")
                            for v in ax.get_xticks()])
        ax.set_xlabel("Partial rank correlation coefficient")
        ax.legend(frameon=False, fontsize=7.5, loc="upper left",
                  handletextpad=0.4, borderaxespad=0.3)
        export_figure(fig, os.path.join(out, "Fig5_PRCC_designs"),
                      formats=["tiff", "pdf"], dpi=300, bbox_inches=None,
                      tiff_rgb=True,
                      provenance={
                          "source": ("new_analysis/static_prcc_reference.json, "
                                     "dynamic_prcc.json, equalwidth_prcc.json"),
                          "designs": "static R0@28C / dynamic peak / equal-width peak",
                          "n_samples": 2000,
                          "seed": 20260816,
                      },
                      write_manifest=True)
        plt.close(fig)


def build(out_dir=None):
    global OUT
    OUT = resolve_out_dir(out_dir)
    os.makedirs(OUT, exist_ok=True)
    print(f"Writing figures to: {OUT}")
    params = baseline_params()

    # ---- Figure 1: R0 vs temperature -------------------------------------
    temps, r0s = r0_temperature_curve(params, T_range=(15, 35), n_points=200)
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, ax = plt.subplots(figsize=(92/25.4, 73/25.4), layout="constrained")
        ax.plot(temps, r0s, color="#0072B2", linewidth=1.6)
        ax.axhline(y=1.0, color="#D55E00", linestyle="--", linewidth=0.9, alpha=0.85)
        key = {}
        for T in [20, 25, 28, 30, 32]:
            r0 = compute_r0(params, T)
            key[T] = r0
            ax.plot(T, r0, "o", color="#0072B2", markersize=4.5,
                    markeredgecolor="white", markeredgewidth=0.6)
            ax.annotate(f"$R_0$={r0:.2f}", xy=(T, r0),
                        xytext=(-14, 8), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7.5,
                        bbox=dict(boxstyle="round,pad=0.2", fc="white",
                                  ec="none", alpha=0.85))
        ax.annotate("threshold $R_0$ = 1", xy=(34.7, 1.0), xytext=(0, 4),
                    textcoords="offset points", fontsize=7.5, color="#D55E00",
                    ha="right", va="bottom")
        ax.set_xlabel("Temperature (\u00b0C)")
        ax.set_ylabel("$R_0$")
        ax.set_xlim(15, 35)
        ax.set_ylim(0, max(r0s) * 1.12)
        VERIFY["r0_curve"] = key
        VERIFY["r0_28"] = key[28]
        exp = export_figure(fig, os.path.join(OUT, "Fig1_R0_temperature"),
                            formats=["tiff", "pdf"], dpi=300, bbox_inches=None,
                            tiff_rgb=True,
                            provenance={"source": "r0.py r0_temperature_curve; seed 20260816",
                                        "params": "baseline_params()",
                                        "key_values": {str(t): round(v, 3) for t, v in key.items()}},
                            write_manifest=True)
        plt.close(fig)

    # ---- Figure 2: PRCC bars ---------------------------------------------
    result = run_sensitivity(n_samples=2000)
    labels = result["param_labels"]; prcc = result["prcc"]
    VERIFY["prcc"] = dict(zip(labels, [round(float(v), 3) for v in prcc]))
    MATH = {"Nv/Nh": "$N_v/N_h$", "beta_hv": r"$\beta_{hv}$", "beta_vh": r"$\beta_{vh}$",
            "mu_v": r"$\mu_v$", "gamma_h": r"$\gamma_h$", "b": "$b$",
            "1/sigma_v": r"$1/\sigma_v$", "sigma_h": r"$\sigma_h$"}
    pretty = [MATH.get(l, l) for l in labels]
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, ax = plt.subplots(figsize=(92/25.4, 79/25.4), layout="constrained")
        colors = ["#0072B2" if v >= 0 else "#D55E00" for v in prcc]
        bars = ax.barh(range(len(labels)), prcc, color=colors,
                       edgecolor="white", linewidth=0.4, height=0.62)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(pretty)
        ax.axvline(x=0, color="#555555", linewidth=0.8)
        ax.set_xlim(-1.15, 1.15)
        ax.invert_yaxis()
        ax.set_xlabel("PRCC")
        for i, v in enumerate(prcc):
            ax.text(v + (0.025 if v >= 0 else -0.025), i, f"{v:+.2f}",
                    va="center", ha="left" if v >= 0 else "right", fontsize=7)
        exp = export_figure(fig, os.path.join(OUT, "Fig2_PRCC_sensitivity"),
                            formats=["tiff", "pdf"], dpi=300, bbox_inches=None,
                            tiff_rgb=True,
                            provenance={"source": "sensitivity.py run(n=2000); seed 20260816",
                                        "prcc": VERIFY["prcc"]},
                            write_manifest=True)
        plt.close(fig)

    # ---- Figure 3: scenario peaks ----------------------------------------
    results = run_all_scenarios()
    names = [r["name"] for r in results]
    peaks = [r["peak_ih"] for r in results]
    VERIFY["peaks"] = dict(zip(names, [round(float(p), 1) for p in peaks]))
    COL = {"Baseline": "#0072B2", "Vector control": "#009E73",
           "Warming +1\u00b0C": "#E69F00", "Warming +2\u00b0C": "#D55E00",
           "Import 4x": "#CC79A7"}
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, ax = plt.subplots(figsize=(95/25.4, 72/25.4), layout="constrained")
        bars = ax.bar(range(len(names)), peaks,
                      color=[COL[n] for n in names], edgecolor="white",
                      linewidth=0.4, width=0.62)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(["Baseline", "Vector\ncontrol", "Warming\n+1\u00b0C",
                            "Warming\n+2\u00b0C", "Import\n4\u00d7"], fontsize=7.5)
        ax.set_ylabel("Peak infectious humans")
        ax.set_ylim(0, max(peaks) * 1.16)
        for bar, peak in zip(bars, peaks):
            ax.text(bar.get_x() + bar.get_width()/2, peak + max(peaks)*0.02,
                    f"{peak:.0f}", ha="center", va="bottom", fontsize=7.5)
        exp = export_figure(fig, os.path.join(OUT, "Fig3_scenario_peaks"),
                            formats=["tiff", "pdf"], dpi=300, bbox_inches=None,
                            tiff_rgb=True,
                            provenance={"source": "scenarios.py run_all_scenarios(); seed 20260816",
                                        "peaks": VERIFY["peaks"]},
                            write_manifest=True)
        plt.close(fig)

    # ---- Figure 4: timeseries --------------------------------------------
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, ax = plt.subplots(figsize=(140/25.4, 80/25.4), layout="constrained")
        series = {}
        for name, color in [("Baseline", "#0072B2"), ("Vector control", "#009E73")]:
            params_s = get_scenario(name)
            t, y = solve(params_s, t_span=(0, 1095), n_points=1000)
            ax.plot(t, y[2], label=name, color=color, linewidth=1.3)
            series[name] = {"peak": round(float(y[2].max()), 1),
                            "peak_day": round(float(t[y[2].argmax()]), 0)}
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Infectious humans")
        ax.set_xlim(0, 1095)
        ax.legend(frameon=False, fontsize=7.5, loc="upper right")
        VERIFY["timeseries"] = series
        exp = export_figure(fig, os.path.join(OUT, "Fig4_incidence_timeseries"),
                            formats=["tiff", "pdf"], dpi=300, bbox_inches=None,
                            tiff_rgb=True,
                            provenance={"source": "scenarios.get_scenario + model.solve; seed 20260816",
                                        "series": series},
                            write_manifest=True)
        plt.close(fig)

    # ---- Figure 5: PRCC under three designs ------------------------------
    fig5_prcc_designs(out_dir=OUT)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build Figures 1-5 (300 dpi TIF + PDF)."
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="output directory (default: auto-detected; env DENGUE_FIG_DIR wins over the default)",
    )
    args = parser.parse_args()

    build(out_dir=args.outdir)
    print("\n===== VERIFY vs manuscript claims =====")
    print("R0(28C) =", round(VERIFY["r0_28"], 3), "  claim: 2.93")
    for T, v in VERIFY["r0_curve"].items():
        print(f"  R0({T}C) = {v:.3f}")
    print("PRCC:", VERIFY["prcc"], "  claim: Nv/Nh +0.95, beta_vh +0.76~0.78, mu_v -0.72, b +0.57")
    print("Peaks:", VERIFY["peaks"], "  claim: 1650/120/1796/1904/1721")
    print("Timeseries:", VERIFY["timeseries"])
