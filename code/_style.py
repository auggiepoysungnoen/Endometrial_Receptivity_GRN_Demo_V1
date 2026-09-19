"""
Shared publication styling and validated color palette for all figures.

Palette source: the project's `dataviz` skill reference instance
(references/palette.md), re-validated for our specific slot combinations via
scripts/validate_palette.js this session (all hard checks PASS). Light mode only
— figures are for print/journal, not a themed UI.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---- Categorical palette (fixed hue order, from the validated reference set) ----
COLOR_FERTILE = "#2a78d6"      # blue, slot 1
COLOR_RIF = "#eb6834"          # orange, slot 2
COLOR_EPITHELIAL = "#1baf7a"   # aqua, slot 3
COLOR_TCELL = "#eda100"        # yellow, slot 4
COLOR_UNK = "#e87ba4"          # magenta, slot 5
COLOR_ACTIVATEDNK = "#008300"  # green, slot 6
COLOR_STROMAL = "#4a3aa7"      # violet, slot 7
# slot 8 (red, #e34948) intentionally skipped for compartment identity — it's
# already DIVERGING_LOW ("lost"/"lost from Fertile") in the network-map edge
# coloring, and reusing it as a node-fill color would collide with that.

CONDITION_COLORS = {"Fertile": COLOR_FERTILE, "RIF": COLOR_RIF}
COMPARTMENT_COLORS = {
    "Epithelial": COLOR_EPITHELIAL,
    "Stromal": COLOR_STROMAL,
    "Tcell": COLOR_TCELL,
    "uNK": COLOR_UNK,
    "ActivatedNK": COLOR_ACTIVATEDNK,
}
# Re-validated for the full 7-color combined set via the dataviz skill's
# scripts/validate_palette.js (light mode): all hard checks PASS (worst
# adjacent CVD separation 5.8-9.1 dE, in the "legal with secondary encoding"
# band — satisfied here since every compartment color always appears next to
# a direct text label/legend entry, never color-alone).

# ---- Sequential ramp (magnitude; one hue, light -> dark) ----
SEQUENTIAL_BLUE = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec",
    "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab",
    "#184f95", "#104281", "#0d366b",
]
SEQUENTIAL_CMAP = mpl.colors.LinearSegmentedColormap.from_list("seq_blue", SEQUENTIAL_BLUE)

# ---- Diverging pair (polarity; two hues + neutral midpoint) ----
DIVERGING_LOW = "#e34948"   # red (e.g. "lost")
DIVERGING_HIGH = "#2a78d6"  # blue (e.g. "gained")
DIVERGING_MID = "#f0efec"   # neutral gray (e.g. "retained")
DIVERGING_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "diverging", [DIVERGING_LOW, DIVERGING_MID, DIVERGING_HIGH]
)

# ---- Chart chrome / ink ----
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"

MUTED_NODE = "#c3c2b7"  # for network figures: non-hub nodes
TARGET_NODE_FILL = "#FFA07A"  # light salmon — target-only gene node fill (network maps)
TF_NODE_FILL = COLOR_STROMAL  # "#4a3aa7" — universal TF node fill across every compartment's
# network map (not compartment-specific): user asked for a single fixed TF color rather than
# one that varies by compartment, matching the blue/violet already used for Stromal.


def apply_style():
    """Call once at the top of every figure script."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.titleweight": "bold",
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "figure.titlesize": 20,
        "figure.titleweight": "bold",
        "figure.dpi": 150,
        "savefig.dpi": 800,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": BASELINE,
        "axes.labelcolor": INK_PRIMARY,
        "axes.titlecolor": INK_PRIMARY,
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "text.color": INK_PRIMARY,
        "grid.color": GRIDLINE,
        "grid.linewidth": 0.7,
        "axes.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "legend.frameon": False,
    })


def save_fig(fig, path):
    fig.savefig(path, dpi=800, bbox_inches="tight", facecolor=SURFACE)
    print(f"Saved {path}")
