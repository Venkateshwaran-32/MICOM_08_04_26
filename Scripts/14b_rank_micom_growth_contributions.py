from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROC_MICOM = PROJECT_ROOT / "Data" / "Processed" / "micom"
RESULTS_MICOM = PROJECT_ROOT / "Results" / "micom"
GROWTH_DIR = RESULTS_MICOM / "growth" / "proper_age_bins"

IN_GROWTH = GROWTH_DIR / "organism_growth_rates_by_agegroup_diet.csv"
IN_ABUND = PROC_MICOM / "agegroup_median_taxonomy_for_micom.csv"
OUT_CSV = GROWTH_DIR / "micom_growth_contribution_ranking_by_agegroup_diet.csv"


def main():
    if not IN_GROWTH.exists():
        raise FileNotFoundError(f"Missing MICOM organism growth file: {IN_GROWTH}. Run Script 14 first.")
    if not IN_ABUND.exists():
        raise FileNotFoundError(f"Missing MICOM abundance file: {IN_ABUND}. Run Script 13 first.")

    growth = pd.read_csv(IN_GROWTH)
    abund = pd.read_csv(IN_ABUND)

    required_growth = {"age_group", "diet", "id", "growth_rate"}
    required_abund = {"sample_id", "id", "abundance"}
    if not required_growth.issubset(growth.columns):
        raise ValueError(f"Growth table must contain columns: {sorted(required_growth)}")
    if not required_abund.issubset(abund.columns):
        raise ValueError(f"Abundance table must contain columns: {sorted(required_abund)}")

    abund = abund.rename(columns={"sample_id": "age_group"})
    meta_cols = [c for c in ["table2_taxon", "profile_taxon", "model_species_id"] if c in abund.columns]

    merged = growth.merge(
        abund[["age_group", "id", "abundance"] + meta_cols],
        on=["age_group", "id"],
        how="left",
    )
    merged["abundance"] = merged["abundance"].fillna(0.0)
    merged["abundance_weighted_contribution"] = merged["growth_rate"] * merged["abundance"]

    merged["rank_by_growth_rate"] = (
        merged.groupby(["age_group", "diet"])["growth_rate"]
        .rank(method="dense", ascending=False)
        .astype(int)
    )
    merged["rank_by_contribution"] = (
        merged.groupby(["age_group", "diet"])["abundance_weighted_contribution"]
        .rank(method="dense", ascending=False)
        .astype(int)
    )

    out_cols = [
        "age_group",
        "diet",
        "id",
        "table2_taxon",
        "profile_taxon",
        "model_species_id",
        "abundance",
        "growth_rate",
        "abundance_weighted_contribution",
        "rank_by_growth_rate",
        "rank_by_contribution",
        "grows_in_diet",
        "growth_rate_full_access",
        "grows_full_access",
        "nongrowing_due_to_competition_or_diet",
    ]
    out_cols = [c for c in out_cols if c in merged.columns]

    merged = merged.sort_values(
        ["age_group", "diet", "rank_by_contribution", "rank_by_growth_rate", "id"]
    )
    merged[out_cols].to_csv(OUT_CSV, index=False)

    print(f"Saved: {OUT_CSV} ({len(merged)} rows)")


if __name__ == "__main__":
    main()
