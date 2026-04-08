from pathlib import Path
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_MICOM = PROJECT_ROOT / "Results" / "micom"
GROWTH_DIR = RESULTS_MICOM / "growth" / "proper_age_bins"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "growth" / "proper_age_bins"
FIG_DIR.mkdir(parents=True, exist_ok=True)

IN_CSV = GROWTH_DIR / "micom_growth_contribution_ranking_by_agegroup_diet.csv"
IN_ABUND = PROJECT_ROOT / "Data" / "Processed" / "micom" / "agegroup_median_taxonomy_for_micom.csv"
OUT_FIG = FIG_DIR / "micom_all_growing_contributors_by_agegroup_diet_proper_age_bins.png"
OUT_GROWTH_FIG = FIG_DIR / "micom_all_growing_species_growth_rates_by_agegroup_diet_proper_age_bins.png"
MIN_GROWTH = 1e-6
TRADEOFF_FRACTION_LABEL = "TF = 0.5"
NOTE_BOX_X = 0.9825
NOTE_BOX_Y = 0.05


def display_name(row: pd.Series) -> str:
    for col in ["table2_taxon", "profile_taxon", "id"]:
        val = row.get(col)
        if pd.notna(val) and str(val).strip():
            return str(val)
    return str(row["id"])


def pretty_age_group(value: str) -> str:
    mapping = {
        "21_40": "Ages 21-40",
        "41_60": "Ages 41-60",
        "61_70": "Ages 61-70",
        "71_80": "Ages 71-80",
        "81_plus": "Ages 81+",
    }
    return mapping.get(str(value), str(value))


def main():
    if not IN_CSV.exists():
        raise FileNotFoundError(f"Missing ranking file: {IN_CSV}. Run Script 14b first.")
    if not IN_ABUND.exists():
        raise FileNotFoundError(f"Missing MICOM abundance file: {IN_ABUND}. Run Script 13 first.")

    df = pd.read_csv(IN_CSV)
    abund = pd.read_csv(IN_ABUND).rename(columns={"sample_id": "age_group"})
    required = {"age_group", "diet", "abundance_weighted_contribution"}
    if not required.issubset(df.columns):
        raise ValueError(f"Ranking file must contain columns: {sorted(required)}")

    top = df[df["growth_rate"] > MIN_GROWTH].copy()
    top["species_label"] = top.apply(display_name, axis=1)
    top["panel"] = top["age_group"].astype(str) + " | " + top["diet"].astype(str)
    species_order = (
        pd.concat(
            [
                top["species_label"],
                df.apply(display_name, axis=1),
            ],
            ignore_index=True,
        )
        .drop_duplicates()
        .tolist()
    )
    palette = [
        "#2f6c8f",
        "#c76b29",
        "#4f8a3c",
        "#a23b72",
        "#8b5e3c",
        "#5b7c99",
        "#b84a3a",
        "#6f9d3c",
        "#7a4fa3",
        "#3c8c84",
    ]
    species_colors = {sp: palette[i % len(palette)] for i, sp in enumerate(species_order)}

    age_groups = sorted(top["age_group"].astype(str).drop_duplicates().tolist())
    panels = [f"{age_group} | {diet}" for age_group in age_groups for diet in ["high_fiber", "western"]]
    fig, axes = plt.subplots(len(age_groups), 2, figsize=(14, max(12, 2.4 * len(age_groups))))
    axes = axes.flatten()

    for ax, panel in zip(axes, panels):
        sub = top[top["panel"] == panel].copy()
        sub = sub.sort_values("abundance_weighted_contribution", ascending=True)
        if sub.empty:
            ax.axis("off")
            continue
        bar_colors = [species_colors[s] for s in sub["species_label"]]
        bars = ax.barh(sub["species_label"], sub["abundance_weighted_contribution"], color=bar_colors)
        age_group, diet = panel.split(" | ", 1)
        ax.set_title(f"{pretty_age_group(age_group)} | {diet} | {TRADEOFF_FRACTION_LABEL}")
        ax.set_xlabel("Abundance-weighted contribution")
        ax.set_ylabel("")
        ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=9)

    for ax in axes[len(panels):]:
        ax.axis("off")

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    fig.text(
        0.985,
        0.025,
        "abundance_weighted_contribution =\nabundance * growth_rate",
        ha="right",
        va="bottom",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.9, "edgecolor": "#888888"},
    )
    fig.text(
        0.015,
        0.025,
        "Key points:\n"
        "1. Faecalibacterium prausnitzii is the dominant contributor across panels.\n"
        "2. Older-age bins show the strongest overall contribution magnitudes.\n"
        "3. Rankings are shaped by both species abundance and MICOM growth rate.",
        ha="left",
        va="bottom",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.9, "edgecolor": "#888888"},
    )
    plt.savefig(OUT_FIG, dpi=220)
    plt.close()
    print(f"Saved: {OUT_FIG}")

    growing = df[df["growth_rate"] > MIN_GROWTH].copy()
    growing["species_label"] = growing.apply(display_name, axis=1)
    growing["panel"] = growing["age_group"].astype(str) + " | " + growing["diet"].astype(str)

    growth_panels = panels
    fig, axes = plt.subplots(len(age_groups), 2, figsize=(14, max(12, 2.4 * len(age_groups))))
    axes = axes.flatten()

    for ax, panel in zip(axes, growth_panels):
        sub = growing[growing["panel"] == panel].copy()
        if sub.empty:
            ax.axis("off")
            continue
        sub = sub.sort_values("growth_rate", ascending=True)
        bar_colors = [species_colors[s] for s in sub["species_label"]]
        bars = ax.barh(sub["species_label"], sub["growth_rate"], color=bar_colors)
        age_group, diet = panel.split(" | ", 1)
        ax.set_title(f"{pretty_age_group(age_group)} | {diet} | {TRADEOFF_FRACTION_LABEL}")
        ax.set_xlabel("Growth rate")
        ax.set_ylabel("")
        ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=9)

        all_species = abund[abund["age_group"].astype(str) == age_group].copy()
        all_species["species_label"] = all_species.apply(display_name, axis=1)
        growing_species = set(sub["species_label"].tolist())
        nongrowing_species = sorted(
            label for label in all_species["species_label"].tolist() if label not in growing_species
        )
        if nongrowing_species:
            note = "Non-growing:\n" + "\n".join(f"{label}: 0.0" for label in nongrowing_species)
        else:
            note = "Non-growing:\nNone"
        ax.text(
            NOTE_BOX_X,
            NOTE_BOX_Y,
            note,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "#888888"},
        )

    for ax in axes[len(growth_panels):]:
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(OUT_GROWTH_FIG, dpi=220)
    plt.close()
    print(f"Saved: {OUT_GROWTH_FIG}")


if __name__ == "__main__":
    main()
