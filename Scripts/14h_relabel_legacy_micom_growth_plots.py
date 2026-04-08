from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_MICOM = PROJECT_ROOT / "Results" / "micom"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "growth"
FIG_DIR.mkdir(parents=True, exist_ok=True)

IN_RANKING = RESULTS_MICOM / "growth" / "micom_growth_contribution_ranking_by_agegroup_diet.csv"
IN_GROWTH = RESULTS_MICOM / "growth" / "organism_growth_rates_by_agegroup_diet.csv"
IN_TRADEOFF_SPECIES = RESULTS_MICOM / "tradeoff_sensitivity" / "tradeoff_sensitivity_species_growth.csv"
IN_TRADEOFF_SUMMARY = RESULTS_MICOM / "tradeoff_sensitivity" / "tradeoff_sensitivity_summary.csv"
IN_SCALED_SPECIES = RESULTS_MICOM / "growth" / "scaled_uptake_species_growth_by_agegroup_diet.csv"
IN_SCALED_COMMUNITY = RESULTS_MICOM / "growth" / "scaled_uptake_community_growth_by_agegroup_diet.csv"

OUT_CONTRIB = FIG_DIR / "micom_top5_growth_contributors_by_agegroup_diet.png"
OUT_GROWTH = FIG_DIR / "micom_top5_growth_rates_by_agegroup_diet.png"
OUT_ALL_GROWING = FIG_DIR / "micom_all_growing_species_growth_rates_by_agegroup_diet.png"
OUT_TRADEOFF_SPECIES = FIG_DIR / "micom_tradeoff_species_growth_curves.png"
OUT_TRADEOFF_COMMUNITY = FIG_DIR / "micom_tradeoff_community_growth_curves.png"
OUT_SCALED5_SPECIES = FIG_DIR / "micom_scaled5x_species_growth_by_agegroup_diet.png"
OUT_SCALED5_COMMUNITY = FIG_DIR / "micom_scaled5x_community_growth_by_agegroup_diet.png"
OUT_SCALED10_SPECIES = FIG_DIR / "micom_scaled10x_species_growth_by_agegroup_diet.png"
OUT_SCALED10_COMMUNITY = FIG_DIR / "micom_scaled10x_community_growth_by_agegroup_diet.png"

MIN_GROWTH = 1e-6


def pretty_age_group(value: str) -> str:
    mapping = {
        "20_40": "Ages 20-39",
        "40_60": "Ages 40-59",
        "60_plus": "Ages 60+",
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


def pretty_species_name(value: str) -> str:
    text = str(value).strip()
    if not text:
        return text
    parts = text.split("_AGORA1_")[0].split("_")
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return text


def scale_description(scale_label: str) -> str:
    if scale_label == "scaled5x":
        return "All diet max_uptake bounds multiplied by 5 relative to the original media files"
    if scale_label == "scaled10x":
        return "All diet max_uptake bounds multiplied by 10 relative to the original media files"
    return scale_label


def build_palette(labels: list[str]) -> dict[str, tuple]:
    cmap = plt.get_cmap("tab10")
    return {label: cmap(i % 10) for i, label in enumerate(sorted(labels))}


def load_abundance_like_table() -> pd.DataFrame:
    df = pd.read_csv(IN_RANKING)
    df["species_label"] = df.apply(display_name, axis=1)
    cols = [c for c in ["age_group", "id", "species_label", "abundance"] if c in df.columns]
    return df[cols].drop_duplicates()


def plot_top_contributors_and_growth():
    df = pd.read_csv(IN_RANKING)
    df["species_label"] = df.apply(display_name, axis=1)
    abundance_like = load_abundance_like_table()

    top = df[df["growth_rate"] > MIN_GROWTH].copy()
    top["panel"] = top["age_group"].astype(str) + " | " + top["diet"].astype(str)
    species_order = pd.concat([top["species_label"], df["species_label"]], ignore_index=True).drop_duplicates().tolist()
    palette = [
        "#2f6c8f", "#c76b29", "#4f8a3c", "#a23b72", "#8b5e3c",
        "#5b7c99", "#b84a3a", "#6f9d3c", "#7a4fa3", "#3c8c84",
    ]
    species_colors = {sp: palette[i % len(palette)] for i, sp in enumerate(species_order)}
    age_groups = sorted(top["age_group"].astype(str).drop_duplicates().tolist())
    panels = [f"{age_group} | {diet}" for age_group in age_groups for diet in ["high_fiber", "western"]]

    fig, axes = plt.subplots(len(age_groups), 2, figsize=(14, max(12, 2.4 * len(age_groups))))
    axes = axes.flatten()
    for ax, panel in zip(axes, panels):
        sub = top[top["panel"] == panel].copy().sort_values("abundance_weighted_contribution", ascending=True).tail(5)
        if sub.empty:
            ax.axis("off")
            continue
        age_group, diet = panel.split(" | ", 1)
        bars = ax.barh(sub["species_label"], sub["abundance_weighted_contribution"], color=[species_colors[s] for s in sub["species_label"]])
        ax.set_title(f"{pretty_age_group(age_group)} | {pretty_diet(diet)} | TF = 0.5")
        ax.set_xlabel("Abundance-weighted contribution")
        ax.set_ylabel("")
        ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=9)
    for ax in axes[len(panels):]:
        ax.axis("off")
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.savefig(OUT_CONTRIB, dpi=220)
    plt.close()

    fig, axes = plt.subplots(len(age_groups), 2, figsize=(14, max(12, 2.4 * len(age_groups))))
    axes = axes.flatten()
    growing = df[df["growth_rate"] > MIN_GROWTH].copy()
    growing["panel"] = growing["age_group"].astype(str) + " | " + growing["diet"].astype(str)
    for ax, panel in zip(axes, panels):
        sub = growing[growing["panel"] == panel].copy().sort_values("growth_rate", ascending=True).tail(5)
        if sub.empty:
            ax.axis("off")
            continue
        age_group, diet = panel.split(" | ", 1)
        bars = ax.barh(sub["species_label"], sub["growth_rate"], color=[species_colors[s] for s in sub["species_label"]])
        ax.set_title(f"{pretty_age_group(age_group)} | {pretty_diet(diet)} | TF = 0.5")
        ax.set_xlabel("Growth rate")
        ax.set_ylabel("")
        ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=9)

        all_species = abundance_like[abundance_like["age_group"].astype(str) == age_group]["species_label"].drop_duplicates().tolist()
        growing_species = set(sub["species_label"].tolist())
        nongrowing = sorted(label for label in all_species if label not in growing_species)
        note = "Non-growing:\n" + ("\n".join(f"{label}: 0.0" for label in nongrowing) if nongrowing else "None")
        ax.text(
            0.98, 0.02, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "#888888"},
        )
    for ax in axes[len(panels):]:
        ax.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_GROWTH, dpi=220)
    plt.close()


def plot_all_growing():
    df = pd.read_csv(IN_RANKING)
    abund = load_abundance_like_table()
    df = df[df["growth_rate"] > MIN_GROWTH].copy()
    df["species_label"] = df.apply(display_name, axis=1)
    species_order = pd.concat([df["species_label"], abund["species_label"]], ignore_index=True).drop_duplicates().tolist()
    palette = build_palette(species_order)
    age_groups = sorted(df["age_group"].astype(str).drop_duplicates().tolist())
    panels = [f"{age_group} | {diet}" for age_group in age_groups for diet in ["high_fiber", "western"]]

    fig, axes = plt.subplots(len(age_groups), 2, figsize=(14, max(12, 2.4 * len(age_groups))))
    axes = axes.flatten()
    for ax, panel in zip(axes, panels):
        age_group, diet = panel.split(" | ", 1)
        sub = df[(df["age_group"].astype(str) == age_group) & (df["diet"].astype(str) == diet)].copy()
        if sub.empty:
            ax.axis("off")
            continue
        sub = sub.sort_values("growth_rate", ascending=True)
        bars = ax.barh(sub["species_label"], sub["growth_rate"], color=[palette[s] for s in sub["species_label"]])
        ax.set_title(f"{pretty_age_group(age_group)} | {pretty_diet(diet)} | TF = 0.5")
        ax.set_xlabel("Growth rate")
        ax.set_ylabel("")
        ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=8)

        all_species = abund[abund["age_group"].astype(str) == age_group]["species_label"].drop_duplicates().tolist()
        growing_species = set(sub["species_label"].tolist())
        nongrowing = [name for name in all_species if name not in growing_species]
        note = "Non-growing:\n" + ("\n".join(f"{name}: 0.0" for name in nongrowing) if nongrowing else "None")
        ax.text(
            0.98, 0.02, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "#888888"},
        )
    for ax in axes[len(panels):]:
        ax.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_ALL_GROWING, dpi=220)
    plt.close()


def plot_tradeoff_curves():
    species = pd.read_csv(IN_TRADEOFF_SPECIES)
    summary = pd.read_csv(IN_TRADEOFF_SUMMARY)
    species["tradeoff_fraction"] = species["tradeoff_fraction"].astype(float)
    summary["tradeoff_fraction"] = summary["tradeoff_fraction"].astype(float)
    species["growth_rate"] = species["growth_rate"].astype(float)
    summary["community_growth_rate"] = summary["community_growth_rate"].astype(float)
    species["species_label"] = species["id"].map(pretty_species_name)
    species["grows_in_diet"] = species["growth_rate"] > MIN_GROWTH

    age_groups = sorted(species["age_group"].astype(str).unique().tolist())
    diets = ["high_fiber", "western"]
    palette = build_palette(species["species_label"].dropna().unique().tolist())

    fig, axes = plt.subplots(len(age_groups), len(diets), figsize=(16, 10.5), sharex=True, sharey=False)
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
                .max().sort_values("growth_rate", ascending=False)["species_label"].tolist()
            )
            for label in ordered_species:
                cur = sub[sub["species_label"] == label].sort_values("tradeoff_fraction")
                ax.plot(cur["tradeoff_fraction"], cur["growth_rate"], marker="o", linewidth=1.8, markersize=4, color=palette[label], alpha=0.9)
            panel_summary = summary[(summary["age_group"] == age_group) & (summary["diet"] == diet)]
            n_growing = int(sub.groupby("tradeoff_fraction")["grows_in_diet"].sum().max())
            n_species = int(panel_summary["n_organisms"].max()) if not panel_summary.empty else int(sub["id"].nunique())
            ax.set_title(f"{pretty_age_group(age_group)} | {pretty_diet(diet)} | {n_growing}/{n_species} growing")
            ax.set_xlim(0.0, 1.08)
            ax.grid(axis="y", alpha=0.25)
            ax.set_xlabel("Tradeoff fraction")
            if col_idx == 0:
                ax.set_ylabel("Individual growth rate")
    fig.tight_layout()
    fig.savefig(OUT_TRADEOFF_SPECIES, dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(len(age_groups), 1, figsize=(10, 10), sharex=True)
    if len(age_groups) == 1:
        axes = [axes]
    colors = {"high_fiber": "#2f6c8f", "western": "#c76b29"}
    for ax, age_group in zip(axes, age_groups):
        sub = summary[summary["age_group"] == age_group].copy()
        for diet in diets:
            cur = sub[sub["diet"] == diet].sort_values("tradeoff_fraction")
            if cur.empty:
                continue
            ax.plot(cur["tradeoff_fraction"], cur["community_growth_rate"], marker="o", linewidth=2.2, markersize=5, color=colors[diet], label=diet)
        ax.set_title(pretty_age_group(age_group))
        ax.set_ylabel("Community growth rate")
        ax.set_xlabel("Tradeoff fraction")
        ax.grid(alpha=0.25)
        ax.legend(frameon=False)
        ref_species = species[(species["age_group"] == age_group) & (species["diet"] == "western") & (species["grows_in_diet"])].copy()
        ordered_species = (
            ref_species.groupby("species_label", as_index=False)["growth_rate"]
            .max().sort_values("growth_rate", ascending=False)["species_label"].tolist()
        )
        if ordered_species:
            anchored = AnchoredText(
                "Growing species:\n" + "\n".join(f"{idx}. {name}" for idx, name in enumerate(ordered_species, start=1)),
                loc="lower right", prop={"size": 8}, frameon=True, borderpad=0.25,
            )
            anchored.patch.set_facecolor("white")
            anchored.patch.set_alpha(0.9)
            ax.add_artist(anchored)
    fig.tight_layout()
    fig.savefig(OUT_TRADEOFF_COMMUNITY, dpi=220)
    plt.close(fig)


def plot_scaled_uptake():
    species = pd.read_csv(IN_SCALED_SPECIES)
    community = pd.read_csv(IN_SCALED_COMMUNITY)
    abund = load_abundance_like_table()
    species = species.merge(abund[["age_group", "id", "species_label"]], on=["age_group", "id"], how="left")
    age_groups = sorted(species["age_group"].drop_duplicates().tolist())
    diets = ["high_fiber", "western"]

    for scale_label, species_out, community_out in [
        ("scaled5x", OUT_SCALED5_SPECIES, OUT_SCALED5_COMMUNITY),
        ("scaled10x", OUT_SCALED10_SPECIES, OUT_SCALED10_COMMUNITY),
    ]:
        sub = species[species["uptake_scale"] == scale_label].copy()
        palette = build_palette(sub["species_label"].dropna().unique().tolist())
        fig, axes = plt.subplots(len(age_groups), len(diets), figsize=(15, 10.5))
        axes = axes.flatten()
        for idx, (age_group, diet) in enumerate([(a, d) for a in age_groups for d in diets]):
            ax = axes[idx]
            panel = sub[(sub["age_group"] == age_group) & (sub["diet"] == diet) & (sub["growth_rate"] > MIN_GROWTH)].copy()
            all_species = abund[abund["age_group"].astype(str) == str(age_group)]["species_label"].drop_duplicates().tolist()
            if panel.empty:
                ax.axis("off")
                continue
            panel = panel.sort_values("growth_rate", ascending=True)
            bars = ax.barh(panel["species_label"], panel["growth_rate"], color=[palette[label] for label in panel["species_label"]])
            ax.set_title(f"{pretty_age_group(age_group)} | {pretty_diet(diet)} | {scale_label}")
            ax.set_xlabel("Growth rate")
            ax.set_ylabel("")
            ax.grid(axis="x", alpha=0.25)
            ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=8)
            nongrowing = [name for name in all_species if name not in set(panel["species_label"].tolist())]
            note = "Non-growing:\n" + ("\n".join(f"{name}: 0.0" for name in nongrowing) if nongrowing else "None")
            ax.text(
                0.98, 0.02, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "#888888"},
            )
        fig.suptitle(f"MICOM Species Growth by Age Group and Diet | {scale_label}\n{scale_description(scale_label)}", y=0.99, fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        fig.savefig(species_out, dpi=220)
        plt.close(fig)

        sub_comm = community[community["uptake_scale"] == scale_label].copy()
        fig, axes = plt.subplots(len(age_groups), 1, figsize=(10, 8.5))
        if len(age_groups) == 1:
            axes = [axes]
        colors = {"high_fiber": "#2f6c8f", "western": "#c76b29"}
        for ax, age_group in zip(axes, age_groups):
            panel = sub_comm[sub_comm["age_group"] == age_group].copy()
            if panel.empty:
                ax.axis("off")
                continue
            for diet in diets:
                cur = panel[panel["diet"] == diet].copy()
                if cur.empty:
                    continue
                bars = ax.bar([diet], cur["community_growth_rate"].tolist(), color=colors[diet], width=0.6, label=diet)
                ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=8)
            ax.set_title(f"{pretty_age_group(age_group)} | {scale_label}")
            ax.set_ylabel("Community growth rate")
            ax.grid(axis="y", alpha=0.25)
            ax.legend(frameon=False)
        fig.suptitle(f"MICOM Community Growth by Age Group and Diet | {scale_label}\n{scale_description(scale_label)}", y=0.99, fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fig.savefig(community_out, dpi=220)
        plt.close(fig)


def main():
    required = [
        IN_RANKING, IN_TRADEOFF_SPECIES, IN_TRADEOFF_SUMMARY, IN_SCALED_SPECIES, IN_SCALED_COMMUNITY,
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required legacy MICOM input files:\n" + "\n".join(missing))

    plot_top_contributors_and_growth()
    plot_all_growing()
    plot_tradeoff_curves()
    plot_scaled_uptake()
    print("Rebuilt legacy MICOM growth figures with explicit age labels.")


if __name__ == "__main__":
    main()
