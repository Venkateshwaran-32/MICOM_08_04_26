from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TRADEOFF_DIR = PROJECT_ROOT / "Results" / "micom" / "tradeoff_sensitivity" / "proper_age_bins"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "growth" / "proper_age_bins"
FIG_DIR.mkdir(parents=True, exist_ok=True)

IN_SPECIES = TRADEOFF_DIR / "tradeoff_sensitivity_species_growth.csv"
IN_SUMMARY = TRADEOFF_DIR / "tradeoff_sensitivity_summary.csv"
IN_ABUND = PROJECT_ROOT / "Data" / "Processed" / "micom" / "agegroup_median_taxonomy_for_micom.csv"

OUT_SPECIES = FIG_DIR / "micom_tradeoff_species_growth_curves_proper_age_bins.png"
OUT_COMMUNITY = FIG_DIR / "micom_tradeoff_community_growth_curves_proper_age_bins.png"

MIN_GROWTH = 1e-6
COMMUNITY_RIGHT_MARGIN = 0.78
COMMUNITY_NOTE_X = 0.95
COMMUNITY_NOTE_LABEL_WIDTH = 24
COMMUNITY_NOTE_Y_GAP = 0.008


def pretty_species_name(value: str) -> str:
    text = str(value).strip()
    if not text:
        return text
    parts = text.split("_AGORA1_")[0].split("_")
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return text


def display_name_from_row(row: pd.Series) -> str:
    for col in ["table2_taxon", "profile_taxon", "id"]:
        val = row.get(col)
        if pd.notna(val) and str(val).strip():
            return str(val)
    return str(row.get("id", ""))


def pretty_age_group(value: str) -> str:
    mapping = {
        "21_40": "Ages 21-40",
        "41_60": "Ages 41-60",
        "61_70": "Ages 61-70",
        "71_80": "Ages 71-80",
        "81_plus": "Ages 81+",
    }
    return mapping.get(str(value), str(value))


def pretty_diet(value: str) -> str:
    return str(value).replace("_", " ")


def abbreviate_species_label(name: str, max_len: int) -> str:
    text = str(name).strip()
    if len(text) <= max_len:
        return text

    parts = text.split()
    if len(parts) >= 2:
        genus = parts[0]
        species = parts[1]
        reserve = len(genus) + 1 + 2  # space + ".."
        species_room = max(1, max_len - reserve)
        return f"{genus} {species[:species_room]}.."

    return text[: max(1, max_len - 2)] + ".."


def fixed_width_species_note(species_names: list[str]) -> str:
    lines = ["Growing species:"]
    for idx, name in enumerate(species_names, start=1):
        clipped = abbreviate_species_label(name, COMMUNITY_NOTE_LABEL_WIDTH)
        lines.append(f"{idx:>2}. {clipped:<{COMMUNITY_NOTE_LABEL_WIDTH}}")
    return "\n".join(lines)


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not IN_SPECIES.exists():
        raise FileNotFoundError(f"Missing species tradeoff table: {IN_SPECIES}. Run Script 14d first.")
    if not IN_SUMMARY.exists():
        raise FileNotFoundError(f"Missing tradeoff summary table: {IN_SUMMARY}. Run Script 14d first.")
    if not IN_ABUND.exists():
        raise FileNotFoundError(f"Missing MICOM abundance table: {IN_ABUND}. Run Script 13 first.")

    species = pd.read_csv(IN_SPECIES)
    summary = pd.read_csv(IN_SUMMARY)
    abund = pd.read_csv(IN_ABUND)

    species_need = {"tradeoff_fraction", "age_group", "diet", "id", "growth_rate"}
    summary_need = {"tradeoff_fraction", "age_group", "diet", "community_growth_rate", "n_growing_in_diet"}
    if not species_need.issubset(species.columns):
        raise ValueError(f"Species tradeoff table missing columns: {sorted(species_need)}")
    if not summary_need.issubset(summary.columns):
        raise ValueError(f"Tradeoff summary table missing columns: {sorted(summary_need)}")

    species["tradeoff_fraction"] = species["tradeoff_fraction"].astype(float)
    summary["tradeoff_fraction"] = summary["tradeoff_fraction"].astype(float)
    species["growth_rate"] = species["growth_rate"].astype(float)
    summary["community_growth_rate"] = summary["community_growth_rate"].astype(float)
    summary["n_growing_in_diet"] = summary["n_growing_in_diet"].astype(int)
    species["species_label"] = species["id"].map(pretty_species_name)
    species["grows_in_diet"] = species["growth_rate"] > MIN_GROWTH
    abund["age_group"] = abund["sample_id"].astype(str)
    abund["species_label"] = abund.apply(display_name_from_row, axis=1)
    return species, summary, abund


def build_species_palette(species: pd.DataFrame) -> dict[str, tuple]:
    labels = sorted(species["species_label"].dropna().unique().tolist())
    cmap = plt.get_cmap("tab10")
    return {label: cmap(i % 10) for i, label in enumerate(labels)}


def plot_species_curves(species: pd.DataFrame, summary: pd.DataFrame, abund: pd.DataFrame):
    age_groups = sorted(species["age_group"].astype(str).unique().tolist())
    diets = ["high_fiber", "western"]
    palette = build_species_palette(species)

    fig, axes = plt.subplots(len(age_groups), len(diets), figsize=(16, 12.5), sharex=True, sharey=False)
    if len(age_groups) == 1:
        axes = [axes]

    for row_idx, age_group in enumerate(age_groups):
        for col_idx, diet in enumerate(diets):
            ax = axes[row_idx][col_idx]
            sub = species[(species["age_group"] == age_group) & (species["diet"] == diet)].copy()
            if sub.empty:
                ax.axis("off")
                continue

            ordered_species = (
                sub.groupby("species_label", as_index=False)["growth_rate"]
                .max()
                .sort_values("growth_rate", ascending=False)["species_label"]
                .tolist()
            )
            top3_species = ordered_species[:3]

            for label in ordered_species:
                cur = sub[sub["species_label"] == label].sort_values("tradeoff_fraction")
                ax.plot(
                    cur["tradeoff_fraction"],
                    cur["growth_rate"],
                    marker="o",
                    linewidth=1.8,
                    markersize=4,
                    color=palette[label],
                    alpha=0.9,
                )
                if label in top3_species and not cur.empty:
                    end_x = float(cur["tradeoff_fraction"].iloc[-1])
                    end_y = float(cur["growth_rate"].iloc[-1])
                    ax.annotate(
                        label,
                        (end_x, end_y),
                        textcoords="offset points",
                        xytext=(6, 0),
                        ha="left",
                        va="center",
                        fontsize=7.5,
                        color=palette[label],
                    )

            n_growing = int(sub.groupby("tradeoff_fraction")["grows_in_diet"].sum().max())
            panel_summary = summary[(summary["age_group"] == age_group) & (summary["diet"] == diet)]
            if panel_summary.empty:
                n_species = int(abund[abund["age_group"].astype(str) == age_group]["species_label"].nunique())
            else:
                n_species = int(panel_summary["n_organisms"].max())
            ax.set_title(f"{pretty_age_group(age_group)} | {pretty_diet(diet)} | {n_growing}/{n_species} growing")
            ax.set_xlim(0.0, 1.08)
            ax.set_xticks([round(x, 1) for x in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]])
            ax.grid(axis="y", alpha=0.25)
            ax.set_xlabel("Tradeoff fraction")
            ax.tick_params(axis="x", labelbottom=True)
            if col_idx == 0:
                ax.set_ylabel("Individual growth rate")

    handles = [
        plt.Line2D([0], [0], color=color, marker="o", linewidth=2, markersize=5, label=label)
        for label, color in palette.items()
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        title="Species",
        frameon=True,
        fontsize=8,
        title_fontsize=9,
        ncol=3,
        columnspacing=1.2,
        handletextpad=0.6,
        borderaxespad=0.2,
    )
    fig.suptitle(
        "MICOM Tradeoff Sensitivity: Individual Species Growth by Age Group and Diet\n"
        "Panels show species with nonzero growth across the tradeoff sweep.",
        y=0.995,
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.82])
    fig.savefig(OUT_SPECIES, dpi=220)
    plt.close(fig)
    print(f"Saved: {OUT_SPECIES}")


def plot_community_curves(summary: pd.DataFrame, species: pd.DataFrame):
    age_groups = sorted(summary["age_group"].astype(str).unique().tolist())
    diets = ["high_fiber", "western"]
    diet_colors = {"high_fiber": "#2f6c8f", "western": "#c76b29"}
    label_offsets = {"high_fiber": (0, -8, "top"), "western": (0, 6, "bottom")}

    fig, axes = plt.subplots(len(age_groups), 1, figsize=(13, 10), sharex=True)
    if len(age_groups) == 1:
        axes = [axes]

    for ax, age_group in zip(axes, age_groups):
        sub = summary[summary["age_group"] == age_group].copy()
        for diet in diets:
            cur = sub[sub["diet"] == diet].sort_values("tradeoff_fraction")
            if cur.empty:
                continue
            label = f"{diet} | growing {int(cur['n_growing_in_diet'].max())}/{int(cur['n_organisms'].max())}"
            ax.plot(
                cur["tradeoff_fraction"],
                cur["community_growth_rate"],
                marker="o",
                linewidth=2.2,
                markersize=5,
                color=diet_colors[diet],
                label=label,
            )
            dx, dy, va = label_offsets[diet]
            for x, y in zip(cur["tradeoff_fraction"], cur["community_growth_rate"]):
                ax.annotate(
                    f"{y:.2f}",
                    (x, y),
                    textcoords="offset points",
                    xytext=(dx, dy),
                    ha="center",
                    va=va,
                    fontsize=7,
                    color=diet_colors[diet],
                )
        ax.set_title(pretty_age_group(age_group))
        ax.set_ylabel("Community growth rate")
        ax.set_xlabel("Tradeoff fraction")
        ax.tick_params(axis="x", labelbottom=True)
        ax.grid(alpha=0.25)
        ax.legend(frameon=False)
        ref_species = species[
            (species["age_group"] == age_group)
            & (species["diet"] == "western")
            & (species["grows_in_diet"])
        ].copy()
        if ref_species.empty:
            ref_species = species[
                (species["age_group"] == age_group)
                & (species["grows_in_diet"])
            ].copy()
        ordered_species = (
            ref_species.groupby("species_label", as_index=False)["growth_rate"]
            .max()
            .sort_values("growth_rate", ascending=False)["species_label"]
            .tolist()
        )
    axes[-1].set_xlim(0.0, 1.0)
    axes[-1].set_xticks([round(x, 1) for x in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]])
    fig.suptitle("MICOM Tradeoff Sensitivity: Community Growth by Age Group and Diet", y=0.995)

    # Reserve right-side figure space for the growing-species note boxes.
    fig.tight_layout(rect=[0, 0, COMMUNITY_RIGHT_MARGIN, 0.97])

    for ax, age_group in zip(axes, age_groups):
        ref_species = species[
            (species["age_group"] == age_group)
            & (species["diet"] == "western")
            & (species["grows_in_diet"])
        ].copy()
        if ref_species.empty:
            ref_species = species[
                (species["age_group"] == age_group)
                & (species["grows_in_diet"])
            ].copy()
        ordered_species = (
            ref_species.groupby("species_label", as_index=False)["growth_rate"]
            .max()
            .sort_values("growth_rate", ascending=False)["species_label"]
            .tolist()
        )
        if not ordered_species:
            continue

        pos = ax.get_position()
        species_note = fixed_width_species_note(ordered_species)
        fig.text(
            COMMUNITY_NOTE_X,
            pos.y0 + COMMUNITY_NOTE_Y_GAP,
            species_note,
            ha="right",
            va="bottom",
            fontsize=8,
            family="monospace",
            multialignment="left",
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "#888888"},
        )

    fig.savefig(OUT_COMMUNITY, dpi=220)
    plt.close(fig)
    print(f"Saved: {OUT_COMMUNITY}")


def main():
    species, summary, abund = load_inputs()
    plot_species_curves(species, summary, abund)
    plot_community_curves(summary, species)


if __name__ == "__main__":
    main()
