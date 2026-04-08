from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_CSV = PROJECT_ROOT / "Results" / "fba" / "community_reaction_fluxes_by_diet.csv"
OUT_DIR = PROJECT_ROOT / "Results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_BAR = OUT_DIR / "active_reaction_counts_by_species_diet.png"
OUT_HEAT = OUT_DIR / "active_reaction_counts_main_species_heatmap.png"


def prettify_species_name(text: str) -> str:
    # Turn model-style IDs into labels that are easier to read on plots.
    return (
        str(text)
        .replace("_AGORA1_03", "")
        .replace("__", "_")
        .replace("_", " ")
    )


def main():
    df = pd.read_csv(IN_CSV)
    required = {"diet", "species", "reaction_id"}
    if not required.issubset(df.columns):
        raise ValueError("Input CSV must contain diet, species, and reaction_id columns.")

    # Keep only rows where we can identify both the species and the reaction.
    df = df.dropna(subset=["species", "reaction_id"]).copy()
    if df.empty:
        raise ValueError("No active reaction rows found in the input CSV.")

    # Count how many unique active reactions each species has under each diet.
    counts = (
        df.groupby(["diet", "species"], as_index=False)["reaction_id"]
        .nunique()
        .rename(columns={"reaction_id": "n_active_reactions"})
    )
    counts["species_label"] = counts["species"].map(prettify_species_name)

    # Plot 1: bar chart showing active reaction counts for all species by diet.
    plot_df = counts.sort_values(["diet", "n_active_reactions"], ascending=[True, False]).copy()
    species_order = (
        counts.groupby("species_label", as_index=False)["n_active_reactions"]
        .mean()
        .sort_values("n_active_reactions", ascending=False)["species_label"]
        .tolist()
    )

    # Reshape the table so rows are species and columns are diets.
    pivot_bar = (
        counts.pivot(index="species_label", columns="diet", values="n_active_reactions")
        .fillna(0)
        .reindex(species_order)
    )

    ax = pivot_bar.plot(kind="bar", figsize=(13, 6), color=["#3b7a57", "#c97b63"])
    ax.set_title("Active reaction counts by species and diet")
    ax.set_xlabel("Species")
    ax.set_ylabel("Number of active reactions")
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()
    plt.savefig(OUT_BAR, dpi=300, bbox_inches="tight")
    plt.close()

    # Plot 2: heatmap for the main organisms (top 4 by mean active reaction count).
    top_species = (
        counts.groupby("species", as_index=False)["n_active_reactions"]
        .mean()
        .sort_values("n_active_reactions", ascending=False)
        .head(4)
    )
    main_species = top_species["species"].tolist()

    heat = (
        counts[counts["species"].isin(main_species)]
        .pivot(index="species_label", columns="diet", values="n_active_reactions")
        .fillna(0)
    )
    heat = heat.reindex([prettify_species_name(s) for s in main_species])

    fig, ax = plt.subplots(figsize=(6, 4.5))
    im = ax.imshow(heat.values, aspect="auto", cmap="YlGnBu")
    ax.set_title("Active reaction counts for main organisms")
    ax.set_xlabel("Diet")
    ax.set_ylabel("Species")
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index)

    for i in range(len(heat.index)):
        for j in range(len(heat.columns)):
            ax.text(j, i, int(heat.iloc[i, j]), ha="center", va="center", color="black")

    fig.colorbar(im, ax=ax, label="Active reactions")
    fig.tight_layout()
    fig.savefig(OUT_HEAT, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {OUT_BAR}")
    print(f"Saved: {OUT_HEAT}")


if __name__ == "__main__":
    main()
