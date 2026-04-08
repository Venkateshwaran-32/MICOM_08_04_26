from pathlib import Path
import pandas as pd


# -----------------------------
# Project paths
# -----------------------------
# Keep paths in one place so you can edit easily later.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
FBA_DIR = PROJECT_ROOT / "Results" / "fba"
FIG_DIR = PROJECT_ROOT / "Results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
PROC_DIR = PROJECT_ROOT / "Data" / "Processed"


# -----------------------------
# Input CSV files (from Script 06)
# -----------------------------
PATHWAY_ACTIVITY_CSV = FBA_DIR / "community_pathway_activity_by_diet.csv"
PATHWAY_COMPARE_CSV = FBA_DIR / "community_pathway_diet_comparison.csv"
BIOMASS_PATHWAY_CSV = FBA_DIR / "community_biomass_associated_pathways_by_diet.csv"
SG90_MEDIAN_CSV = PROC_DIR / "sg90_median_abundance_by_taxon.csv"


# -----------------------------
# Figure output files
# -----------------------------
OUT_HEATMAP_TOP20 = FIG_DIR / "pathway_activity_heatmap_top20.png"
OUT_DELTA_TOP20 = FIG_DIR / "pathway_diet_delta_top20.png"
OUT_DELTA_ALISTIPES = FIG_DIR / "pathway_diet_delta_top15_alistipes_shahii.png"
OUT_DELTA_FAECALI = FIG_DIR / "pathway_diet_delta_top15_faecalibacterium_prausnitzii.png"
OUT_BIOMASS_TOP15 = FIG_DIR / "biomass_associated_pathways_by_diet_top15.png"
OUT_SG90_MEDIAN = FIG_DIR / "sg90_median_taxa_abundance.png"
OUT_SG90_WEIGHT = FIG_DIR / "sg90_median_taxa_normalized_weight.png"


def load_required_csv(path: Path) -> pd.DataFrame:
    """
    Read a CSV and stop with a clear error if missing.
    This helps beginners quickly understand what to run first.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Missing required input: {path}\n"
            "Run Script 06 first so this figure script has data to plot."
        )
    return pd.read_csv(path)


def make_pathway_activity_heatmap(df_pathway: pd.DataFrame):
    """
    Plot top pathways by total flux activity across diets as a heatmap.
    """
    import matplotlib.pyplot as plt

    top = (
        df_pathway.groupby("pathway", as_index=False)["sum_abs_flux"]
        .sum()
        .sort_values("sum_abs_flux", ascending=False)
        .head(20)["pathway"]
        .tolist()
    )
    heat = (
        df_pathway[df_pathway["pathway"].isin(top)]
        .groupby(["pathway", "diet"], as_index=False)["sum_abs_flux"]
        .sum()
        .pivot(index="pathway", columns="diet", values="sum_abs_flux")
        .fillna(0.0)
    )
    if heat.empty:
        print("Skipping heatmap: no pathway rows found.")
        return

    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(heat.values, aspect="auto")
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=45, ha="right")
    ax.set_title("Top Pathway Activity by Diet")
    fig.colorbar(im, ax=ax, label="sum_abs_flux")
    fig.tight_layout()
    fig.savefig(OUT_HEATMAP_TOP20, dpi=200)
    plt.close(fig)
    print(f"Saved: {OUT_HEATMAP_TOP20}")


def make_pathway_delta_top20(df_compare: pd.DataFrame):
    """
    Plot top pathway changes between diets (all species combined rows).
    """
    import matplotlib.pyplot as plt

    if df_compare.empty:
        print("Skipping top20 delta plot: comparison table is empty.")
        return

    top_delta = (
        df_compare.assign(abs_delta=lambda d: d["delta_high_fiber_minus_western"].abs())
        .sort_values("abs_delta", ascending=False)
        .head(20)
        .copy()
    )
    if top_delta.empty:
        print("Skipping top20 delta plot: no rows after filtering.")
        return

    top_delta["label"] = top_delta["species"] + " | " + top_delta["pathway"]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(top_delta["label"], top_delta["delta_high_fiber_minus_western"])
    ax.set_title("Top Pathway Flux Changes (high_fiber - western)")
    ax.set_xlabel("delta sum_abs_flux")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUT_DELTA_TOP20, dpi=200)
    plt.close(fig)
    print(f"Saved: {OUT_DELTA_TOP20}")


def make_species_delta_plot(df_compare: pd.DataFrame, species: str, out_file: Path):
    """
    Plot top 15 pathway deltas for one species.
    Green bars = increased in high_fiber.
    Red bars = decreased in high_fiber.
    """
    import matplotlib.pyplot as plt

    sub = df_compare[df_compare["species"] == species].copy()
    if sub.empty:
        print(f"Skipping species plot ({species}): no rows found.")
        return

    sub["abs_delta"] = sub["delta_high_fiber_minus_western"].abs()
    top = sub.sort_values("abs_delta", ascending=False).head(15).copy()
    top = top.sort_values("delta_high_fiber_minus_western", ascending=True)

    fig, ax = plt.subplots(figsize=(10, 6.5))
    colors = ["#2E7D32" if x > 0 else "#C62828" for x in top["delta_high_fiber_minus_western"]]
    ax.barh(top["pathway"], top["delta_high_fiber_minus_western"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_title(f"Top Pathway Flux Shifts: {species} (high_fiber - western)")
    ax.set_xlabel("delta sum_abs_flux")
    fig.tight_layout()
    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    print(f"Saved: {out_file}")


def make_biomass_pathway_plot(df_biomass: pd.DataFrame):
    """
    Plot pathways most associated with biomass production, split by diet.
    Score is summed across species to give a community-level view.
    """
    import matplotlib.pyplot as plt

    if df_biomass.empty:
        print("Skipping biomass pathway plot: table is empty.")
        return

    agg = (
        df_biomass.groupby(["diet", "pathway"], as_index=False)["biomass_association_score"]
        .sum()
    )

    top_pathways = (
        agg.groupby("pathway", as_index=False)["biomass_association_score"]
        .sum()
        .sort_values("biomass_association_score", ascending=False)
        .head(15)["pathway"]
    )

    plot_df = agg[agg["pathway"].isin(top_pathways)].copy()
    pivot = (
        plot_df.pivot(index="pathway", columns="diet", values="biomass_association_score")
        .fillna(0.0)
    )
    pivot = pivot.loc[top_pathways]

    fig, ax = plt.subplots(figsize=(11, 8))
    idx = range(len(pivot.index))
    width = 0.38
    west = pivot["western"] if "western" in pivot.columns else [0.0] * len(pivot.index)
    hfib = pivot["high_fiber"] if "high_fiber" in pivot.columns else [0.0] * len(pivot.index)

    ax.barh([i - width / 2 for i in idx], west, height=width, label="western", color="#C62828")
    ax.barh([i + width / 2 for i in idx], hfib, height=width, label="high_fiber", color="#2E7D32")
    ax.set_yticks(list(idx))
    ax.set_yticklabels(pivot.index)
    ax.invert_yaxis()
    ax.set_xlabel("Summed biomass_association_score (across species)")
    ax.set_title("Pathways Most Associated With Biomass by Diet")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_BIOMASS_TOP15, dpi=220)
    plt.close(fig)
    print(f"Saved: {OUT_BIOMASS_TOP15}")


def make_sg90_taxa_plots(df_sg90: pd.DataFrame):
    """
    Plot SG90 median abundance and normalized objective weights for modeled taxa.
    """
    import matplotlib.pyplot as plt

    required = {"table2_taxon", "median_abundance", "normalized_weight"}
    if not required.issubset(df_sg90.columns):
        print(f"Skipping SG90 taxa plots: missing columns {sorted(required)}")
        return

    d = df_sg90.copy().sort_values("median_abundance", ascending=False)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(d["table2_taxon"], d["median_abundance"], color="#1565C0")
    ax.set_title("SG90 Median Taxa Abundance (Modeled Species)")
    ax.set_ylabel("median_abundance")
    ax.set_xlabel("taxon")
    ax.set_yscale("log")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(OUT_SG90_MEDIAN, dpi=220)
    plt.close(fig)
    print(f"Saved: {OUT_SG90_MEDIAN}")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(d["table2_taxon"], d["normalized_weight"], color="#2E7D32")
    ax.set_title("SG90 Normalized Weights Used in Community FBA")
    ax.set_ylabel("normalized_weight")
    ax.set_xlabel("taxon")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(OUT_SG90_WEIGHT, dpi=220)
    plt.close(fig)
    print(f"Saved: {OUT_SG90_WEIGHT}")


def main():
    """
    Build all currently used project figures from Script 06 outputs.
    """
    try:
        import matplotlib  # noqa: F401
    except Exception as e:
        raise RuntimeError(
            "matplotlib is required for figure generation. "
            "Install with: .venv\\Scripts\\python -m pip install matplotlib"
        ) from e

    df_pathway = load_required_csv(PATHWAY_ACTIVITY_CSV)
    df_compare = load_required_csv(PATHWAY_COMPARE_CSV)
    df_biomass = load_required_csv(BIOMASS_PATHWAY_CSV)
    df_sg90 = load_required_csv(SG90_MEDIAN_CSV)

    make_pathway_activity_heatmap(df_pathway)
    make_pathway_delta_top20(df_compare)
    make_species_delta_plot(
        df_compare,
        "Alistipes_shahii_WAL_8301_AGORA1_03",
        OUT_DELTA_ALISTIPES,
    )
    make_species_delta_plot(
        df_compare,
        "Faecalibacterium_prausnitzii_M21_2__AGORA1_03",
        OUT_DELTA_FAECALI,
    )
    make_biomass_pathway_plot(df_biomass)
    make_sg90_taxa_plots(df_sg90)

    print(f"\nAll figures saved to: {FIG_DIR}")


if __name__ == "__main__":
    main()
