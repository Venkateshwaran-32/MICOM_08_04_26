from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "Results" / "fba" / "pathway_traces"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PNG = OUT_DIR / "western__Alistipes_shahii_WAL_8301_AGORA1_03__Glycolysis_gluconeogenesis_curated_pipeline.png"


BG = "#fcfbf7"
TEXT = "#15324a"
MET_FILL = "#eef6fb"
MET_EDGE = "#9ab9cf"
RXN_FILL = "#1f5a89"
RXN_TEXT = "white"
ARROW = "#54718a"
BRANCH = "#8a5a2b"


def draw_box(ax, x, y, w, h, text, fc, ec, tc, fontsize=12, weight="normal"):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor=fc,
        edgecolor=ec,
        linewidth=1.5,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize, color=tc, weight=weight)


def draw_arrow(ax, x1, y1, x2, y2, color=ARROW):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=18,
            linewidth=1.8,
            color=color,
        )
    )


def metabolite(ax, x, y, label):
    draw_box(ax, x, y, 0.22, 0.05, label, MET_FILL, MET_EDGE, TEXT, fontsize=13, weight="bold")


def reaction(ax, x, y, label):
    draw_box(ax, x, y, 0.50, 0.06, label, RXN_FILL, RXN_FILL, RXN_TEXT, fontsize=11.5, weight="bold")


def main():
    fig, ax = plt.subplots(figsize=(11, 17))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.05, 0.975, "Alistipes shahii", fontsize=22, weight="bold", color=TEXT, va="top")
    ax.text(0.05, 0.945, "Glycolysis / gluconeogenesis", fontsize=18, weight="bold", color=TEXT, va="top")
    ax.text(0.05, 0.918, "Curated main-chain trace from active community-FBA reactions under western diet", fontsize=11, color="#536879", va="top")

    ax.text(0.08, 0.87, "Start", fontsize=14, weight="bold", color=TEXT)
    ax.text(0.77, 0.87, "End", fontsize=14, weight="bold", color=TEXT)

    left_x = 0.08
    rxn_x = 0.36
    y_positions = [0.82, 0.73, 0.64, 0.55, 0.46, 0.37, 0.28, 0.19]

    metabolite_labels = [
        "f6p[c]",
        "fdp[c]",
        "dhap[c] + g3p[c]",
        "g3p[c]",
        "13dpg[c]",
        "3pg[c]",
        "2pg[c]",
        "pep[c]",
        "pyr[c]",
    ]

    reaction_labels = [
        "PFK / PFK(ppi)\nmake fdp[c]",
        "FBA\nfdp[c] -> dhap[c] + g3p[c]",
        "TPI\ndhap[c] -> g3p[c]",
        "GAPD\ng3p[c] -> 13dpg[c]",
        "PGK\n13dpg[c] -> 3pg[c]",
        "PGM\n3pg[c] -> 2pg[c]",
        "ENO\n2pg[c] -> pep[c]",
        "PYK\npep[c] -> pyr[c]",
    ]

    # First metabolite
    metabolite(ax, left_x, y_positions[0], metabolite_labels[0])
    draw_arrow(ax, left_x + 0.22, y_positions[0] + 0.025, rxn_x, y_positions[0] + 0.03)
    reaction(ax, rxn_x, y_positions[0] - 0.005, reaction_labels[0])

    for i in range(len(reaction_labels)):
        met_y = y_positions[i]
        rxn_y = met_y - 0.005
        next_met_y = met_y - 0.09

        if i < len(reaction_labels) - 1:
            metabolite(ax, left_x, next_met_y, metabolite_labels[i + 1])
            draw_arrow(ax, rxn_x + 0.25, rxn_y, rxn_x + 0.25, next_met_y + 0.05)
            draw_arrow(ax, rxn_x, next_met_y + 0.025, left_x + 0.22, next_met_y + 0.025)
            reaction(ax, rxn_x, next_met_y - 0.005, reaction_labels[i + 1])
        else:
            metabolite(ax, 0.70, next_met_y, metabolite_labels[-1])
            draw_arrow(ax, rxn_x + 0.50, rxn_y + 0.03, 0.70, next_met_y + 0.025)

    # Branch section
    ax.text(0.08, 0.10, "Separate branch observed in the same active pathway set", fontsize=12.5, weight="bold", color=BRANCH)
    metabolite(ax, 0.08, 0.04, "g6p[c]")
    draw_box(ax, 0.36, 0.035, 0.30, 0.06, "PGMT\ng6p[c] -> g1p[c]", BRANCH, BRANCH, "white", fontsize=11.5, weight="bold")
    metabolite(ax, 0.72, 0.04, "g1p[c]")
    draw_arrow(ax, 0.30, 0.065, 0.36, 0.065, color=BRANCH)
    draw_arrow(ax, 0.66, 0.065, 0.72, 0.065, color=BRANCH)

    ax.text(
        0.05,
        0.005,
        "Read each reaction box left-to-right exactly as written. This figure starts at internal f6p[c], not extracellular glucose.",
        fontsize=10.5,
        color=TEXT,
        va="bottom",
    )

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
