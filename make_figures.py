"""Restyled figure builder for the Ecological Modelling transfer.

Standards: scientific-visualization skill (style_context default + Okabe-Ito),
figure-qa rules (no in-figure titles, U+2212, TIF LZW RGB + PDF, 300 dpi).
Data provenance: model.py / r0.py / sensitivity.py / scenarios.py, seed 20260816.
"""
import os, sys

SRC = os.path.dirname(os.path.abspath(__file__))
SKILL = "/Users/wangtianyu/.workbuddy/skills/scientific-visualization"
sys.path.insert(0, SRC)
sys.path.insert(0, os.path.join(SKILL, "scripts"))
sys.path.insert(0, os.path.join(SKILL, "assets"))

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

OUT = os.path.abspath(os.path.join(SRC, "..", "..", "02_Figures"))
os.makedirs(OUT, exist_ok=True)

VERIFY = {}

def build():
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

if __name__ == "__main__":
    build()
    print("\n===== VERIFY vs manuscript claims =====")
    print("R0(28C) =", round(VERIFY["r0_28"], 3), "  claim: 2.93")
    for T, v in VERIFY["r0_curve"].items():
        print(f"  R0({T}C) = {v:.3f}")
    print("PRCC:", VERIFY["prcc"], "  claim: Nv/Nh +0.95, beta_vh +0.76~0.78, mu_v -0.72, b +0.57")
    print("Peaks:", VERIFY["peaks"], "  claim: 1650/120/1796/1904/1721")
    print("Timeseries:", VERIFY["timeseries"])
