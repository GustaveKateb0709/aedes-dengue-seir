"""Self-contained Matplotlib style preset used by the figure builder.

This is a trimmed, dependency-free copy of the subset of the
``scientific-visualization`` preset set that ``make_figures.py`` actually uses:
one publication style (``default``) and one colour cycle
(``okabe_ito_on_white``).  It depends only on Matplotlib -- no external asset
files, JSON profiles or private paths -- so the supplementary package is fully
portable.

The rcParams values are byte-for-byte the ones used to produce the figures in
the manuscript, so re-running the pipeline reproduces the same artwork.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Iterator

# ── rcParams of the publication preset ("default") ──────────────────────────

BASE_STYLE: dict[str, Any] = {
    "figure.dpi": 100,
    "figure.facecolor": "white",
    "figure.autolayout": False,
    "figure.constrained_layout.use": False,
    "font.size": 8,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.6,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "axes.labelweight": "normal",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "axes.axisbelow": True,
    "axes.grid": False,
    "xtick.major.size": 3,
    "xtick.minor.size": 2,
    "xtick.major.width": 0.6,
    "xtick.minor.width": 0.5,
    "xtick.labelsize": 7,
    "xtick.direction": "out",
    "ytick.major.size": 3,
    "ytick.minor.size": 2,
    "ytick.major.width": 0.6,
    "ytick.minor.width": 0.5,
    "ytick.labelsize": 7,
    "ytick.direction": "out",
    "lines.linewidth": 1.4,
    "lines.markersize": 4,
    "lines.markeredgewidth": 0.6,
    "patch.linewidth": 0.6,
    "legend.fontsize": 7,
    "legend.frameon": False,
    "savefig.dpi": 300,
    "savefig.format": "pdf",
    # "tight" changes the exported physical page dimensions. Opt in explicitly.
    "savefig.bbox": "standard",
    "savefig.pad_inches": 0.05,
    "savefig.transparent": False,
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "image.cmap": "viridis",
    # Type 42 keeps TrueType fonts editable/embedded in PDF/PS.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    # SVG text remains text but depends on font availability at render time.
    "svg.fonttype": "none",
}

# ── Colour cycle: Okabe-Ito on white (colour-blind safe) ────────────────────

PALETTES: dict[str, list[str]] = {
    "okabe_ito_on_white": ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#000000"],
}


def get_style(
    style_name: str = "default",
    *,
    palette_name: str = "okabe_ito_on_white",
) -> dict[str, Any]:
    """Return a validated plain rcParams dictionary without global mutation."""
    if style_name != "default":
        raise ValueError(
            f"unknown style {style_name!r}; this package ships only 'default'"
        )
    if palette_name not in PALETTES:
        raise ValueError(
            f"unknown palette {palette_name!r}; "
            f"available: {', '.join(sorted(PALETTES))}"
        )
    style = dict(BASE_STYLE)
    # Stored separately because a Matplotlib Cycler is not JSON serializable.
    style["_palette_colors"] = PALETTES[palette_name]
    return style


def _matplotlib_style(style: dict[str, Any]) -> dict[str, Any]:
    import matplotlib as mpl

    prepared = dict(style)
    colors = prepared.pop("_palette_colors")
    prepared["axes.prop_cycle"] = mpl.cycler(color=colors)
    # rcParams validation catches stale or misspelled keys.
    validated = mpl.RcParams()
    validated.update(prepared)
    return dict(validated)


@contextmanager
def style_context(
    style_name: str = "default",
    *,
    palette_name: str = "okabe_ito_on_white",
) -> Iterator[dict[str, Any]]:
    """Temporarily apply the preset without leaking global rcParams changes."""
    import matplotlib as mpl

    style = _matplotlib_style(get_style(style_name, palette_name=palette_name))
    with mpl.rc_context(style):
        yield {
            "style": style_name,
            "palette": palette_name,
            "notice": (
                "General-purpose visual preset; review the rendered figure "
                "in context."
            ),
        }
