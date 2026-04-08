from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_CSV = PROJECT_ROOT / "Results" / "micom" / "lysine_butyrate" / "lysine_related_fluxes_long.csv"
OUT_DIR = PROJECT_ROOT / "Results" / "micom" / "lysine_butyrate"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "lysine_butyrate"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_SUMMARY = OUT_DIR / "lysine_agegroup_metabolism_by_pathway.csv"
OUT_FIG = FIG_DIR / "lysine_agegroup_metabolism_by_pathway.png"

AGE_ORDER = ["21_40", "41_60", "61_70", "71_80", "81_plus"]
DIET_ORDER = ["western", "high_fiber"]
COLORS = {
    "Lysine metabolism": "#1f5a89",
    "Butanoate metabolism": "#9a5d21",
    "Folate metabolism": "#2f7d4a",
    "Cell wall biosynthesis": "#8c4f7d",
}
SUMMARY_TEXT = (
    "Summary\n"
    "1. Lysine metabolism is the dominant pathway signal across age groups.\n"
    "2. Older bins retain strong lysine flux, with 71_80 high_fiber the peak.\n"
    "3. Folate and cell wall pathways are secondary context around the main lysine trend."
)
SOURCE_TEXT = (
    "Sources\n"
    "lysine_related_fluxes_long.csv\n"
    "lysine_agegroup_metabolism_by_pathway.csv"
)


def wrap_title(text: str, width: int = 18) -> str:
    words = str(text).split()
    lines = []
    current = []
    length = 0
    for word in words:
        add_len = len(word) + (1 if current else 0)
        if length + add_len > width and current:
            lines.append(" ".join(current))
            current = [word]
            length = len(word)
        else:
            current.append(word)
            length += add_len
    if current:
        lines.append(" ".join(current))
    return "\n".join(lines)


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    lys = df[df["is_lysine_related"] == True].copy()
    if lys.empty:
        return pd.DataFrame(
            columns=["diet", "age_group", "pathway", "sum_abs_flux", "n_rows"]
        )

    out = (
        lys.groupby(["diet", "age_group", "pathway"], as_index=False)
        .agg(sum_abs_flux=("abs_flux", "sum"), n_rows=("reaction_id", "count"))
    )
    out["age_group"] = pd.Categorical(out["age_group"], categories=AGE_ORDER, ordered=True)

    pathway_order = (
        out.groupby("pathway", as_index=False)["sum_abs_flux"]
        .sum()
        .sort_values("sum_abs_flux", ascending=False)["pathway"]
        .tolist()
    )
    out["pathway"] = pd.Categorical(out["pathway"], categories=pathway_order, ordered=True)
    out = out.sort_values(["diet", "age_group", "pathway"]).reset_index(drop=True)
    return out


def plot_summary(summary: pd.DataFrame, out_path: Path):
    if summary.empty:
        raise ValueError("No lysine-related rows found to plot.")

    pathways = summary["pathway"].cat.categories.tolist() if hasattr(summary["pathway"], "cat") else sorted(summary["pathway"].unique().tolist())
    fig, axes = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
    fig.patch.set_facecolor("#fcfbf7")

    x = np.arange(len(AGE_ORDER))
    n_pathways = max(1, len(pathways))
    total_width = 0.78
    bar_width = total_width / n_pathways

    legend_handles = None
    legend_labels = None

    for ax, diet in zip(axes, DIET_ORDER):
        sub = summary[summary["diet"] == diet].copy()
        pivot = sub.pivot(index="age_group", columns="pathway", values="sum_abs_flux").reindex(AGE_ORDER).fillna(0.0)

        bars = []
        for i, pathway in enumerate(pathways):
            offset = (i - (n_pathways - 1) / 2) * bar_width
            color = COLORS.get(pathway, "#5f7284")
            bar = ax.bar(x + offset, pivot.get(pathway, pd.Series([0.0] * len(AGE_ORDER), index=AGE_ORDER)).values, width=bar_width * 0.92, color=color, label=pathway)
            bars.append(bar)

        if legend_handles is None:
            legend_handles, legend_labels = ax.get_legend_handles_labels()

        ax.set_title(f"{diet.replace('_', ' ').title()} diet", fontsize=12, loc="left", pad=8, color="#173042", weight="bold")
        ax.set_ylabel("Sum abs flux", fontsize=10, color="#173042")
        ax.grid(axis="y", color="#ddd6c8", linewidth=0.8, alpha=0.7)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#bfb6a5")
        ax.spines["bottom"].set_color("#bfb6a5")
        ax.tick_params(axis="y", labelsize=9, colors="#173042")
        ax.tick_params(axis="x", labelsize=10, colors="#173042")

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(AGE_ORDER)

    fig.suptitle("Lysine-Related Metabolism by Age Group", fontsize=16, color="#173042", weight="bold", y=0.98)
    fig.text(0.07, 0.945, "Pathway-level sums from lysine_related_fluxes_long.csv", fontsize=10, color="#5c7280")

    if legend_handles and legend_labels:
        wrapped_labels = [wrap_title(label, width=20) for label in legend_labels]
        fig.legend(
            legend_handles,
            wrapped_labels,
            ncol=min(4, len(wrapped_labels)),
            loc="upper center",
            bbox_to_anchor=(0.5, 0.905),
            frameon=False,
            fontsize=9,
            handlelength=1.8,
            columnspacing=1.6,
        )

    fig.text(
        0.07,
        0.08,
        SUMMARY_TEXT,
        fontsize=8.5,
        color="#173042",
        va="top",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#fff8e8", "edgecolor": "#d8c9a8", "linewidth": 0.8},
    )
    fig.text(
        0.93,
        0.08,
        SOURCE_TEXT,
        fontsize=8,
        color="#5c7280",
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#f7f4ec", "edgecolor": "#d8d2c3", "linewidth": 0.8},
    )

    plt.subplots_adjust(top=0.83, left=0.09, right=0.97, bottom=0.2, hspace=0.28)
    plt.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main():
    if not IN_CSV.exists():
        raise FileNotFoundError(f"Missing input: {IN_CSV}")

    df = pd.read_csv(IN_CSV)
    required = {"age_group", "diet", "pathway", "abs_flux", "reaction_id", "is_lysine_related"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input file must contain columns: {sorted(required)}")

    summary = build_summary(df)
    summary.to_csv(OUT_SUMMARY, index=False)
    plot_summary(summary, OUT_FIG)

    print(f"Saved summary: {OUT_SUMMARY}")
    print(f"Saved figure: {OUT_FIG}")


if __name__ == "__main__":
    main()
