from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SUPP_DIR = PROJECT_ROOT / "Data" / "Supplementary"
INPUT_DIR = PROJECT_ROOT / "Results" / "inputs"
PROC_DIR = PROJECT_ROOT / "Data" / "Processed"
INPUT_DIR.mkdir(parents=True, exist_ok=True)
PROC_DIR.mkdir(parents=True, exist_ok=True)

BETA_CSV = SUPP_DIR / "beta_coeffecient_by_species_from_gut_metagenomes_article.csv"
SG90_INPUT = INPUT_DIR / "sg90_median_input_for_fba.csv"

OUT_MAP = PROC_DIR / "beta_age_taxon_mapping.csv"
OUT_INPUT = INPUT_DIR / "beta_age_input_for_fba.csv"


def normalize_nonnegative(series: pd.Series) -> pd.Series:
    total = float(series.sum())
    if total <= 0:
        return pd.Series([0.0] * len(series), index=series.index)
    return series / total


def main():
    beta = pd.read_csv(BETA_CSV)
    sg90 = pd.read_csv(SG90_INPUT)

    required_beta = {"taxon", "beta_age"}
    required_sg90 = {"table2_taxon", "profile_taxon", "model_file", "model_species_id"}
    if not required_beta.issubset(beta.columns):
        raise ValueError("Beta CSV must contain taxon and beta_age columns.")
    if not required_sg90.issubset(sg90.columns):
        raise ValueError("SG90 input must contain table2_taxon, profile_taxon, model_file, and model_species_id.")

    beta["taxon"] = beta["taxon"].astype(str).str.strip()
    sg90["table2_taxon"] = sg90["table2_taxon"].astype(str).str.strip()

    # Use the existing SG90 input as the mapping reference from paper taxon names to model species IDs.
    mapped = beta.merge(
        sg90[["table2_taxon", "profile_taxon", "model_file", "model_species_id"]],
        left_on="taxon",
        right_on="table2_taxon",
        how="left",
    )
    mapped["match_status"] = mapped["model_species_id"].notna().map({True: "matched", False: "unmatched"})
    mapped.to_csv(OUT_MAP, index=False)

    matched = mapped[mapped["match_status"] == "matched"].copy()
    if matched.empty:
        raise ValueError("No beta taxa matched the modeled species list.")

    # Keep the raw beta, but also create two nonnegative weight schemes for downstream FBA use.
    matched["beta_positive_only"] = matched["beta_age"].clip(lower=0.0)
    min_beta = float(matched["beta_age"].min())
    matched["beta_shifted_nonnegative"] = matched["beta_age"] - min_beta

    matched["normalized_weight_positive_only"] = normalize_nonnegative(matched["beta_positive_only"])
    matched["normalized_weight_shifted"] = normalize_nonnegative(matched["beta_shifted_nonnegative"])

    out = matched[
        [
            "taxon",
            "profile_taxon",
            "model_file",
            "model_species_id",
            "beta_age",
            "beta_positive_only",
            "beta_shifted_nonnegative",
            "normalized_weight_positive_only",
            "normalized_weight_shifted",
        ]
    ].sort_values("normalized_weight_positive_only", ascending=False)

    out.to_csv(OUT_INPUT, index=False)

    print(f"Saved: {OUT_MAP}")
    print(f"Saved: {OUT_INPUT}")
    print(f"Matched taxa: {len(matched)}")
    print(f"Unmatched taxa: {int((mapped['match_status'] == 'unmatched').sum())}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
