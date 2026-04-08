from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GROWTH_DIR = PROJECT_ROOT / "Results" / "micom" / "growth" / "proper_age_bins"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "growth" / "proper_age_bins"
FIG_DIR.mkdir(parents=True, exist_ok=True)

IN_SPECIES = GROWTH_DIR / "scaled_uptake_species_growth_by_agegroup_diet.csv"
IN_COMMUNITY = GROWTH_DIR / "scaled_uptake_community_growth_by_agegroup_diet.csv"
IN_ABUND = PROJECT_ROOT / "Data" / "Processed" / "micom" / "agegroup_median_taxonomy_for_micom.csv"

MIN_GROWTH = 1e-6
NOTE_BOX_X = 0.9825
NOTE_BOX_Y = 0.05


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


def display_name(row: pd.Series) -> str:
    for col in ["table2_taxon", "profile_taxon", "id"]:
        val = row.get(col)
        if pd.notna(val) and str(val).strip():
            return str(val)
    return str(row.get("id", ""))


def build_palette(labels: list[str]) -> dict[str, tuple]:
    cmap = plt.get_cmap("tab10")
    return {label: cmap(i % 10) for i, label in enumerate(sorted(labels))}


def scale_description(scale_label: str) -> str:
    if scale_label == "scaled5x":
        return "All diet max_uptake bounds multiplied by 5 relative to the original media files"
    if scale_label == "scaled10x":
        return "All diet max_uptake bounds multiplied by 10 relative to the original media files"
    return scale_label


def plot_species(df: pd.DataFrame, abund: pd.DataFrame, scale_label: str):
    sub = df[df["uptake_scale"] == scale_label].copy()
    sub["species_label"] = sub["species_label"].astype(str)
    palette = build_palette(sub["species_label"].drop_duplicates().tolist())
    age_groups = sorted(sub["age_group"].drop_duplicates().tolist())
    diets = ["high_fiber", "western"]

    fig, axes = plt.subplots(len(age_groups), len(diets), figsize=(15, 12))
    axes = axes.flatten()

    for idx, (age_group, diet) in enumerate([(a, d) for a in age_groups for d in diets]):
        ax = axes[idx]
        panel = sub[(sub["age_group"] == age_group) & (sub["diet"] == diet) & (sub["growth_rate"] > MIN_GROWTH)].copy()
        all_species = abund[abund["age_group"].astype(str) == age_group]["species_label"].drop_duplicates().tolist()

        if panel.empty:
            ax.axis("off")
            continue

        panel = panel.sort_values("growth_rate", ascending=True)
        colors = [palette[label] for label in panel["species_label"]]
        bars = ax.barh(panel["species_label"], panel["growth_rate"], color=colors)
        ax.set_title(f"{pretty_age_group(age_group)} | {pretty_diet(diet)} | {scale_label}")
        ax.set_xlabel("Growth rate")
        ax.set_ylabel("")
        ax.grid(axis="x", alpha=0.25)
        ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=8)

        growing_species = set(panel["species_label"].tolist())
        nongrowing = [name for name in all_species if name not in growing_species]
        note = "Non-growing:\n" + ("\n".join(f"{name}: 0.0" for name in nongrowing) if nongrowing else "None")
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

    handles = [
        plt.Line2D([0], [0], color=palette[label], linewidth=6, label=label)
        for label in sorted(palette)
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.95),
        title="Species",
        frameon=True,
        fontsize=8,
        title_fontsize=9,
        ncol=3,
        columnspacing=1.2,
        handlelength=1.8,
    )
    out = FIG_DIR / f"micom_{scale_label}_species_growth_by_agegroup_diet_proper_age_bins.png"
    fig.suptitle(
        f"MICOM Species Growth by Age Group and Diet | {scale_label}\n{scale_description(scale_label)}",
        y=0.99,
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.83])
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"Saved: {out}")


def plot_community(df: pd.DataFrame, scale_label: str):
    sub = df[df["uptake_scale"] == scale_label].copy()
    age_groups = sorted(sub["age_group"].drop_duplicates().tolist())
    diets = ["high_fiber", "western"]
    colors = {"high_fiber": "#2f6c8f", "western": "#c76b29"}

    fig, axes = plt.subplots(len(age_groups), 1, figsize=(10, 10))
    if len(age_groups) == 1:
        axes = [axes]

    for ax, age_group in zip(axes, age_groups):
        panel = sub[sub["age_group"] == age_group].copy()
        if panel.empty:
            ax.axis("off")
            continue
        for diet in diets:
            cur = panel[panel["diet"] == diet].copy()
            if cur.empty:
                continue
            x = [diet]
            y = cur["community_growth_rate"].tolist()
            bars = ax.bar(x, y, color=colors[diet], width=0.6, label=diet)
            ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=8)
        ax.set_title(f"{pretty_age_group(age_group)} | {scale_label}")
        ax.set_ylabel("Community growth rate")
        ax.grid(axis="y", alpha=0.25)
        ax.legend(frameon=False)

    out = FIG_DIR / f"micom_{scale_label}_community_growth_by_agegroup_diet_proper_age_bins.png"
    fig.suptitle(
        f"MICOM Community Growth by Age Group and Diet | {scale_label}\n{scale_description(scale_label)}",
        y=0.99,
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"Saved: {out}")


def main():
    if not IN_SPECIES.exists():
        raise FileNotFoundError(f"Missing scaled species growth table: {IN_SPECIES}. Run Script 14f first.")
    if not IN_COMMUNITY.exists():
        raise FileNotFoundError(f"Missing scaled community growth table: {IN_COMMUNITY}. Run Script 14f first.")
    if not IN_ABUND.exists():
        raise FileNotFoundError(f"Missing MICOM abundance table: {IN_ABUND}. Run Script 13 first.")

    species = pd.read_csv(IN_SPECIES)
    community = pd.read_csv(IN_COMMUNITY)
    abund = pd.read_csv(IN_ABUND).rename(columns={"sample_id": "age_group"})

    species["growth_rate"] = species["growth_rate"].astype(float)
    species["species_label"] = species["id"].astype(str)
    species = species.merge(
        abund[["age_group", "id", "table2_taxon", "profile_taxon"]],
        on=["age_group", "id"],
        how="left",
    )
    species["species_label"] = species.apply(display_name, axis=1)
    abund["species_label"] = abund.apply(display_name, axis=1)

    for scale_label in ["scaled5x", "scaled10x"]:
        plot_species(species, abund, scale_label)
        plot_community(community, scale_label)


if __name__ == "__main__":
    main()
