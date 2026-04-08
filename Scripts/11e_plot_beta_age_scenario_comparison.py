from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCENARIO_DIR = PROJECT_ROOT / "Results" / "fba" / "scenarios"
FIG_DIR = PROJECT_ROOT / "Results" / "figures"
OUT_DIR = PROJECT_ROOT / "Results" / "fba"

FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

SCENARIOS = [
    "beta_age_positive_only",
    "beta_age_shifted",
    "beta_age_raw_signed",
]
DIETS = ["western", "high_fiber"]
MIN_GROWTH = 1e-9

OUT_FIG = FIG_DIR / "beta_age_scenario_comparison_biomass_by_diet.png"
OUT_SUMMARY = OUT_DIR / "beta_age_scenario_comparison_summary.csv"

BG = "#fcfbf7"
TEXT = "#173042"
MUTED = "#5c7280"
GRID = "#ddd6c8"

SCENARIO_LABELS = {
    "beta_age_positive_only": "Positive only",
    "beta_age_shifted": "Shifted",
    "beta_age_raw_signed": "Raw signed",
}

SCENARIO_NOTES = {
    "beta_age_positive_only": "Sparse solution dominated by top positive-beta taxa",
    "beta_age_shifted": "Shifted nonnegative weights allow broader participation",
    "beta_age_raw_signed": "Literal negative coefficients produce the strongest winner-take-all behavior",
}

SPECIES_COLORS = {
    "Alistipes shahii": "#1f5a89",
    "Faecalibacterium prausnitzii": "#2f7d4a",
    "Bacteroides xylanisolvens": "#c76b29",
    "Bacteroides dorei": "#6c8b3c",
    "Escherichia coli": "#9a4f7a",
    "Ruminococcus torques": "#b45f04",
    "Parabacteroides merdae": "#657d95",
    "Alistipes onderdonkii": "#7d7d7d",
    "Klebsiella pneumoniae": "#8f6a3b",
    "Bilophila unclassified": "#9aa3ad",
}


def short_species_name(value: str) -> str:
    parts = str(value).split("_")
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return str(value).replace("_", " ")


def fmt_value(value: float) -> str:
    value = float(value)
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}"
    if abs(value) >= 1:
        return f"{value:.2f}"
    if abs(value) > 0:
        return f"{value:.3f}"
    return "0"


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    species_frames = []
    community_frames = []
    for scenario in SCENARIOS:
        species_fp = SCENARIO_DIR / f"{scenario}_species_biomass_by_diet.csv"
        community_fp = SCENARIO_DIR / f"{scenario}_community_summary_by_diet.csv"
        if not species_fp.exists():
            raise FileNotFoundError(f"Missing species biomass file: {species_fp}")
        if not community_fp.exists():
            raise FileNotFoundError(f"Missing community summary file: {community_fp}")
        species_frames.append(pd.read_csv(species_fp))
        community_frames.append(pd.read_csv(community_fp))

    species_df = pd.concat(species_frames, ignore_index=True)
    community_df = pd.concat(community_frames, ignore_index=True)
    species_df["species_label"] = species_df["species"].map(short_species_name)
    return species_df, community_df


def build_species_order(species_df: pd.DataFrame) -> list[str]:
    order = (
        species_df.groupby("species_label", as_index=False)["biomass_flux"]
        .max()
        .sort_values(["biomass_flux", "species_label"], ascending=[False, True])["species_label"]
        .tolist()
    )
    return order


def build_summary_table(species_df: pd.DataFrame, community_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scenario in SCENARIOS:
        for diet in DIETS:
            sub = species_df[(species_df["scenario"] == scenario) & (species_df["diet"] == diet)].copy()
            sub = sub.sort_values("biomass_flux", ascending=False).reset_index(drop=True)
            positive = sub[sub["biomass_flux"] > MIN_GROWTH].copy()
            comm = community_df[(community_df["scenario"] == scenario) & (community_df["diet"] == diet)]
            objective = float(comm["community_objective"].iloc[0]) if not comm.empty else 0.0

            main = positive.iloc[0] if len(positive) >= 1 else None
            second = positive.iloc[1] if len(positive) >= 2 else None
            rows.append(
                {
                    "scenario": scenario,
                    "diet": diet,
                    "community_objective": objective,
                    "n_positive_growth_species": int(len(positive)),
                    "main_grower": main["species_label"] if main is not None else "None",
                    "main_biomass": float(main["biomass_flux"]) if main is not None else 0.0,
                    "secondary_grower": second["species_label"] if second is not None else "None",
                    "secondary_biomass": float(second["biomass_flux"]) if second is not None else 0.0,
                }
            )
    return pd.DataFrame(rows)


def source_box_text() -> str:
    names = []
    for scenario in SCENARIOS:
        names.append(f"{scenario}_community_summary_by_diet.csv")
        names.append(f"{scenario}_species_biomass_by_diet.csv")
    return "Sources\n" + "\n".join(names)


def row_summary_text(summary_df: pd.DataFrame, scenario: str) -> str:
    bits = []
    for diet in DIETS:
        row = summary_df[(summary_df["scenario"] == scenario) & (summary_df["diet"] == diet)].iloc[0]
        bits.append(
            f"{diet}: obj {row['community_objective']:.3f}, "
            f"{int(row['n_positive_growth_species'])} growers, "
            f"winner {row['main_grower']}"
        )
    return " | ".join(bits)


def make_plot(species_df: pd.DataFrame, summary_df: pd.DataFrame):
    species_order = build_species_order(species_df)
    max_flux = float(species_df["biomass_flux"].max())
    x_max = max(1.0, max_flux * 1.18)

    fig, axes = plt.subplots(len(SCENARIOS), len(DIETS), figsize=(16, 12.8), sharex=True)
    fig.patch.set_facecolor(BG)

    for i, scenario in enumerate(SCENARIOS):
        row_text = row_summary_text(summary_df, scenario)
        axes[i, 0].text(
            0.0,
            1.18,
            f"{SCENARIO_LABELS[scenario]}: {SCENARIO_NOTES[scenario]}",
            transform=axes[i, 0].transAxes,
            fontsize=11,
            color=TEXT,
            weight="bold",
            ha="left",
        )
        axes[i, 0].text(
            0.0,
            1.08,
            row_text,
            transform=axes[i, 0].transAxes,
            fontsize=8.8,
            color=MUTED,
            ha="left",
        )

        for j, diet in enumerate(DIETS):
            ax = axes[i, j]
            ax.set_facecolor(BG)
            sub = (
                species_df[(species_df["scenario"] == scenario) & (species_df["diet"] == diet)]
                .set_index("species_label")
                .reindex(species_order)
                .reset_index()
            )
            sub["biomass_flux"] = sub["biomass_flux"].fillna(0.0)

            y = np.arange(len(species_order))
            colors = []
            for _, row in sub.iterrows():
                base = SPECIES_COLORS.get(row["species_label"], "#7f8b96")
                colors.append(base if row["biomass_flux"] > MIN_GROWTH else "#d7d3cc")

            bars = ax.barh(y, sub["biomass_flux"], color=colors, edgecolor="none", height=0.65)
            ax.set_yticks(y)
            if j == 0:
                ax.set_yticklabels(species_order, fontsize=8, color=TEXT)
            else:
                ax.set_yticklabels([])
            ax.invert_yaxis()
            ax.set_xlim(0, x_max)
            ax.set_title(diet.replace("_", " "), fontsize=10.5, color=TEXT, pad=8)
            ax.grid(axis="x", color=GRID, linewidth=0.8, alpha=0.7)
            ax.set_axisbelow(True)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_color("#bfb6a5")
            ax.spines["bottom"].set_color("#bfb6a5")
            ax.tick_params(axis="x", labelsize=8, colors=TEXT)
            ax.tick_params(axis="y", length=0)

            labels = [fmt_value(v) if v > MIN_GROWTH else "" for v in sub["biomass_flux"]]
            ax.bar_label(bars, labels=labels, padding=3, fontsize=7, color=TEXT)

            if i == len(SCENARIOS) - 1:
                ax.set_xlabel("Species biomass flux", fontsize=9, color=TEXT)

    fig.suptitle(
        "Beta-Age COBRApy FBA Scenario Comparison",
        fontsize=17,
        color=TEXT,
        weight="bold",
        y=0.985,
    )
    fig.text(
        0.02,
        0.035,
        "Interpretation\n"
        "1. Positive only remains sparse and favors top positive-beta taxa.\n"
        "2. Shifted weights allow broader participation while retaining a strong winner.\n"
        "3. Raw signed is the most extreme winner-take-all scenario because negative-beta taxa are explicitly penalized.",
        fontsize=9,
        color=TEXT,
        va="bottom",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#fff8e8", "edgecolor": "#d8c9a8", "linewidth": 0.8},
    )
    fig.text(
        0.98,
        0.035,
        source_box_text(),
        fontsize=8,
        color=MUTED,
        ha="right",
        va="bottom",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#f7f4ec", "edgecolor": "#d8d2c3", "linewidth": 0.8},
    )
    plt.tight_layout(rect=[0, 0.14, 1, 0.95], h_pad=2.0, w_pad=1.5)
    plt.savefig(OUT_FIG, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main():
    species_df, community_df = load_data()
    summary_df = build_summary_table(species_df, community_df)
    summary_df.to_csv(OUT_SUMMARY, index=False)
    print(f"Saved: {OUT_SUMMARY}")
    make_plot(species_df, summary_df)
    print(f"Saved: {OUT_FIG}")


if __name__ == "__main__":
    main()
