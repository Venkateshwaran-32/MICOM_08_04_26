from pathlib import Path
import sys
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]


class Validator:
    def __init__(self):
        self.errors = []
        self.warnings = []

    def ok(self, msg: str):
        print(f"[OK] {msg}")

    def warn(self, msg: str):
        self.warnings.append(msg)
        print(f"[WARN] {msg}")

    def fail(self, msg: str):
        self.errors.append(msg)
        print(f"[FAIL] {msg}")


def read_csv_checked(v: Validator, path: Path, required_cols: set[str] | None = None) -> pd.DataFrame | None:
    if not path.exists():
        v.fail(f"Missing file: {path}")
        return None
    try:
        df = pd.read_csv(path)
    except Exception as e:
        v.fail(f"Could not read CSV: {path} ({type(e).__name__})")
        return None

    if required_cols is not None:
        missing = required_cols - set(df.columns)
        if missing:
            v.fail(f"Missing columns in {path.name}: {sorted(missing)}")
            return None
    v.ok(f"Loaded {path.name} ({len(df)} rows)")
    return df


def check_allcohort_branch(v: Validator):
    proc = PROJECT_ROOT / "Data" / "Processed"
    inp = PROJECT_ROOT / "Results" / "inputs"
    fba = PROJECT_ROOT / "Results" / "fba"
    sc = fba / "scenarios"
    fig = PROJECT_ROOT / "Results" / "figures"

    # Inputs
    p_age = proc / "allcohort_agebin_median_abundance_by_taxon.csv"
    p_input = inp / "allcohort_agebin_input_for_fba.csv"
    p_sum = proc / "allcohort_agebin_input_summary.csv"

    df_age = read_csv_checked(v, p_age, {"age_group", "model_species_id", "median_abundance", "normalized_weight"})
    df_in = read_csv_checked(v, p_input, {"age_group", "model_species_id", "normalized_weight"})
    df_sum = read_csv_checked(v, p_sum, {"age_group", "n_subjects"})

    expected_bins = {"21_40", "41_60", "61_70", "71_80", "81_90"}

    if df_sum is not None:
        bins = set(df_sum["age_group"].astype(str).unique())
        missing = expected_bins - bins
        if missing:
            v.fail(f"Missing age bins in summary: {sorted(missing)}")
        else:
            v.ok("All expected age bins are present in allcohort summary")

    if df_in is not None:
        sums = df_in.groupby("age_group", as_index=False)["normalized_weight"].sum()
        bad = sums[(sums["normalized_weight"] - 1.0).abs() > 1e-6]
        if len(bad) > 0:
            v.fail(f"Normalized weights do not sum to ~1 for bins: {bad['age_group'].astype(str).tolist()}")
        else:
            v.ok("Normalized weights sum to ~1 for every allcohort age bin")

    # Scenario outputs
    base_name = "allcohort_equal_abundance"
    req_base = [
        sc / f"{base_name}_community_summary_by_diet.csv",
        sc / f"{base_name}_species_biomass_by_diet.csv",
        sc / f"{base_name}_community_exchange_fluxes_by_diet.csv",
    ]
    for p in req_base:
        if not p.exists():
            v.fail(f"Missing baseline scenario file: {p}")
        else:
            v.ok(f"Found {p.name}")

    # Check each age-bin scenario files
    for b in sorted(expected_bins):
        tag = f"allcohort_agebin_{b}"
        files = [
            sc / f"{tag}_community_summary_by_diet.csv",
            sc / f"{tag}_species_biomass_by_diet.csv",
            sc / f"{tag}_community_exchange_fluxes_by_diet.csv",
            sc / f"{tag}_vs_allcohort_equal_exchange_flux_comparison_by_diet.csv",
        ]
        for p in files:
            if not p.exists():
                v.fail(f"Missing scenario file: {p}")
            else:
                v.ok(f"Found {p.name}")

        # community summary sanity
        s = read_csv_checked(v, sc / f"{tag}_community_summary_by_diet.csv", {"diet", "status", "community_objective"})
        if s is not None:
            diets = set(s["diet"].astype(str).unique())
            if diets != {"western", "high_fiber"}:
                v.fail(f"{tag} summary diets mismatch: {sorted(diets)}")
            if not (s["status"].astype(str) == "optimal").all():
                v.warn(f"{tag} has non-optimal statuses")

    # Master comparison
    p_master = fba / "allcohort_agebin_vs_equal_exchange_comparisons.csv"
    m = read_csv_checked(v, p_master, {"scenario", "diet", "exchange_id", "delta_flux_vs_equal"})
    if m is not None:
        scens = sorted(set(m["scenario"].astype(str).unique()))
        if len(scens) != 5:
            v.warn(f"Expected 5 age-bin scenarios in master comparison, found {len(scens)}")
        if len(m) == 0:
            v.fail("Master comparison table is empty")

    # Figures
    fig_required = [
        fig / "allcohort_agebin_vs_equal_flux_delta_top30.png",
        fig / "allcohort_agebin_vs_equal_flux_delta_heatmap_top40.png",
        fig / "allcohort_agebin_top3_species_by_diet.png",
        fig / "allcohort_agebin_nonzero_species_count_by_diet.png",
        fig / "allcohort_agebin_cohort_sample_counts_used.png",
    ]
    for p in fig_required:
        if not p.exists():
            v.fail(f"Missing figure: {p}")
        else:
            v.ok(f"Found figure {p.name}")


def check_micom_branch(v: Validator):
    mic = PROJECT_ROOT / "Results" / "micom"
    req = [
        mic / "growth" / "organism_growth_rates_by_agegroup_diet.csv",
        mic / "growth" / "community_growth_summary_by_agegroup_diet.csv",
        mic / "pathway_flux" / "reaction_fluxes_long_by_agegroup_diet.csv",
        mic / "growth" / "organism_growth_capability_check.csv",
    ]
    for p in req:
        if not p.exists():
            v.warn(f"MICOM output missing (branch may not have been rerun): {p}")
        else:
            v.ok(f"Found MICOM output {p}")


def main():
    v = Validator()
    print("== Validating All-Cohort Age-Bin Branch ==")
    check_allcohort_branch(v)
    print("\n== Checking MICOM Branch Presence ==")
    check_micom_branch(v)

    print("\n== Validation Summary ==")
    print(f"Errors: {len(v.errors)}")
    print(f"Warnings: {len(v.warnings)}")

    if v.errors:
        print("Validation FAILED")
        sys.exit(1)
    print("Validation PASSED")


if __name__ == "__main__":
    main()
