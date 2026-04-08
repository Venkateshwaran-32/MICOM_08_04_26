from pathlib import Path
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# Script 16: find and plot pathways different between age groups
# Inputs:
# - Results/micom/pathway_flux/reaction_fluxes_long_by_agegroup_diet.csv
# Outputs:
# - Results/micom/pathway_flux/pathway_flux_by_agegroup_diet.csv
# - Results/micom/pathway_flux/differential_pathways_by_agegroup_diet.csv
# - Results/micom/pathway_flux/pairwise_pathway_flux_differences.csv
# - Results/figures/micom/pathway_flux/differential_pathways_heatmap_top25.png
# - Results/figures/micom/pathway_flux/differential_pathways_top20_range.png
# Run:
# - .venv\Scripts\python Scripts\16_find_differential_pathways.py
# Expected runtime:
# - ~2 to 15 seconds
# -------------------------------------------------------------------

# This script answers:
# "Which pathways differ most between age groups?"
# It computes per-pathway flux totals, then derives differential scores
# (range/std across age groups) for each diet.

PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_FLUX = PROJECT_ROOT / "Results" / "micom" / "pathway_flux" / "reaction_fluxes_long_by_agegroup_diet.csv"
OUT_DIR = PROJECT_ROOT / "Results" / "micom" / "pathway_flux"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "pathway_flux"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_PATHWAY = OUT_DIR / "pathway_flux_by_agegroup_diet.csv"
OUT_DIFF = OUT_DIR / "differential_pathways_by_agegroup_diet.csv"
OUT_PAIRWISE = OUT_DIR / "pairwise_pathway_flux_differences.csv"
FIG_HEAT = FIG_DIR / "differential_pathways_heatmap_top25.png"
FIG_BAR = FIG_DIR / "differential_pathways_top20_range.png"

FLUX_EPS = 1e-9


def summarize_by_pathway(df: pd.DataFrame) -> pd.DataFrame:
    # Collapse reaction-level fluxes to pathway-level totals.
    d = df[df["abs_flux"] > FLUX_EPS].copy()
    out = (
        d.groupby(["diet", "age_group", "pathway"], as_index=False)
        .agg(sum_abs_flux=("abs_flux", "sum"), n_reactions=("reaction_id", "count"))
    )
    return out


def make_diff_scores(pathway_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for diet, sub in pathway_df.groupby("diet"):
        # Wide table per diet: rows=pathways, cols=age groups.
        piv = sub.pivot(index="pathway", columns="age_group", values="sum_abs_flux").fillna(0.0)
        if piv.shape[1] < 2:
            continue

        for pathway, vals in piv.iterrows():
            vals_arr = vals.values
            # range_sum_abs_flux is the main differential score used in plots.
            rows.append(
                {
                    "diet": diet,
                    "pathway": pathway,
                    "max_sum_abs_flux": float(vals_arr.max()),
                    "min_sum_abs_flux": float(vals_arr.min()),
                    "range_sum_abs_flux": float(vals_arr.max() - vals_arr.min()),
                    "std_sum_abs_flux": float(vals_arr.std()),
                }
            )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["diet", "range_sum_abs_flux"], ascending=[True, False])
    return out


def make_pairwise(pathway_df: pd.DataFrame) -> pd.DataFrame:
    # Explicit pairwise differences (age_group_b - age_group_a) for each pathway.
    rows = []
    for (diet, pathway), sub in pathway_df.groupby(["diet", "pathway"]):
        grp_vals = {str(r["age_group"]): float(r["sum_abs_flux"]) for _, r in sub.iterrows()}
        groups = sorted(grp_vals.keys())
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                g1, g2 = groups[i], groups[j]
                delta = grp_vals[g2] - grp_vals[g1]
                rows.append(
                    {
                        "diet": diet,
                        "pathway": pathway,
                        "age_group_a": g1,
                        "age_group_b": g2,
                        "sum_abs_flux_a": grp_vals[g1],
                        "sum_abs_flux_b": grp_vals[g2],
                        "delta_b_minus_a": delta,
                        "abs_delta_b_minus_a": abs(delta),
                    }
                )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["diet", "abs_delta_b_minus_a"], ascending=[True, False])
    return out


def plot_heatmap(pathway_df: pd.DataFrame):
    top = (
        # Select globally active pathways first, then compare their age/diet pattern.
        pathway_df.groupby("pathway", as_index=False)["sum_abs_flux"]
        .sum()
        .sort_values("sum_abs_flux", ascending=False)
        .head(25)["pathway"]
        .tolist()
    )

    heat = (
        pathway_df[pathway_df["pathway"].isin(top)]
        .assign(age_diet=lambda d: d["age_group"].astype(str) + "|" + d["diet"].astype(str))
        .pivot(index="pathway", columns="age_diet", values="sum_abs_flux")
        .fillna(0.0)
    )
    if heat.empty:
        print("Skipping heatmap: no rows.")
        return

    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(heat.values, aspect="auto")
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=45, ha="right")
    ax.set_title("Top Pathway Fluxes Across Age Groups and Diets")
    fig.colorbar(im, ax=ax, label="sum abs flux")
    fig.tight_layout()
    fig.savefig(FIG_HEAT, dpi=220)
    plt.close(fig)
    print(f"Saved: {FIG_HEAT}")


def plot_top_ranges(diff_df: pd.DataFrame):
    if diff_df.empty:
        print("Skipping top-range bar plot: no differential pathways.")
        return

    top = diff_df.sort_values("range_sum_abs_flux", ascending=False).head(20).copy()
    top["label"] = top["diet"].astype(str) + " | " + top["pathway"].astype(str)
    top = top.sort_values("range_sum_abs_flux", ascending=True)

    fig, ax = plt.subplots(figsize=(11, 8))
    ax.barh(top["label"], top["range_sum_abs_flux"])
    ax.set_title("Top Differential Pathways by Age Group (Range of Sum Abs Flux)")
    ax.set_xlabel("range(sum_abs_flux) across age groups")
    fig.tight_layout()
    fig.savefig(FIG_BAR, dpi=220)
    plt.close(fig)
    print(f"Saved: {FIG_BAR}")


def main():
    if not IN_FLUX.exists():
        raise FileNotFoundError(f"Missing input: {IN_FLUX}. Run Script 14 first.")

    df = pd.read_csv(IN_FLUX)
    required = {"age_group", "diet", "pathway", "reaction_id", "abs_flux"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input flux table must contain columns: {sorted(required)}")

    by_path = summarize_by_pathway(df)
    diff = make_diff_scores(by_path)
    pairwise = make_pairwise(by_path)

    # Output files:
    # - OUT_PATHWAY: raw pathway totals by age group and diet
    # - OUT_DIFF: ranked differential scores
    # - OUT_PAIRWISE: pairwise deltas between age groups
    by_path.to_csv(OUT_PATHWAY, index=False)
    diff.to_csv(OUT_DIFF, index=False)
    pairwise.to_csv(OUT_PAIRWISE, index=False)

    if by_path["age_group"].nunique() < 2:
        print("Only one age group present in inputs. Differential analysis across age groups is not possible yet.")
    else:
        plot_heatmap(by_path)
        plot_top_ranges(diff)

    print(f"Saved: {OUT_PATHWAY}")
    print(f"Saved: {OUT_DIFF}")
    print(f"Saved: {OUT_PAIRWISE}")


if __name__ == "__main__":
    main()
