from pathlib import Path
import re
import pandas as pd

# -------------------------------------------------------------------
# Script 13: prepare MICOM input tables from supplementary data
# Inputs:
# - Data/Supplementary/supplementary_data_file_1.xlsx (Metadata sheet)
# - Data/Supplementary/taxonomic_profiles_filtered.csv
# - Metadata/models_manifest.csv/models_manifest.csv
# Outputs:
# - Data/Processed/micom/taxon_to_model_mapping_for_micom.csv
# - Data/Processed/micom/subjects_with_age_groups.csv
# - Data/Processed/micom/agegroup_median_taxonomy_for_micom.csv
# - Data/Processed/micom/agegroup_input_summary.csv
# Run:
# - .venv\Scripts\python Scripts\13_prepare_micom_inputs.py
# Expected runtime:
# - ~2 to 10 seconds
# -------------------------------------------------------------------

# This script prepares MICOM-ready abundance inputs per age group.
# Output rows follow MICOM's expected columns: sample_id, id, abundance, file.

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SUPP_DIR = PROJECT_ROOT / "Data" / "Supplementary"
PROC_DIR = PROJECT_ROOT / "Data" / "Processed"
MICOM_DIR = PROC_DIR / "micom"
META_MANIFEST = PROJECT_ROOT / "Metadata" / "models_manifest.csv" / "models_manifest.csv"

META_XLSX = SUPP_DIR / "supplementary_data_file_1.xlsx"
TAXON_CSV = SUPP_DIR / "taxonomic_profiles_filtered.csv"

OUT_MAPPING = MICOM_DIR / "taxon_to_model_mapping_for_micom.csv"
OUT_AGE_SUBJECTS = MICOM_DIR / "subjects_with_age_groups.csv"
OUT_AGE_MEDIAN = MICOM_DIR / "agegroup_median_taxonomy_for_micom.csv"
OUT_AGE_SUMMARY = MICOM_DIR / "agegroup_input_summary.csv"


def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", text)


def to_profile_label(taxon: str) -> str:
    return "s__" + taxon.replace(" ", "_")


def age_bin(age: float):
    # Age bins aligned to your SG90 analysis plan.
    if pd.isna(age):
        return None
    if 21 <= age <= 40:
        return "21_40"
    if 41 <= age <= 60:
        return "41_60"
    if 61 <= age <= 70:
        return "61_70"
    if 71 <= age <= 80:
        return "71_80"
    if age >= 81:
        return "81_plus"
    return None


def normalize(values: pd.Series) -> pd.Series:
    # MICOM abundance weights should sum to 1 inside each community/sample.
    total = float(values.sum())
    if total <= 0:
        return pd.Series([1.0 / len(values)] * len(values), index=values.index)
    return values / total


def first_existing_column(df: pd.DataFrame, candidates):
    # Accept slight column-name variation across supplementary files.
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"None of the expected columns found: {candidates}")


def build_mapping(manifest: pd.DataFrame, taxa_labels: set[str]) -> pd.DataFrame:
    # Connect taxonomy names to local model files and stable model IDs.
    rows = []
    for _, r in manifest.iterrows():
        table_taxon = str(r["table2_taxon"]).strip()
        model_file = str(r["local_file"]).strip()
        abs_model_file = (PROJECT_ROOT / model_file).resolve()
        model_stem = Path(model_file).stem

        profile_taxon = to_profile_label(table_taxon)
        # One manual proxy used in your project mapping rules.
        if table_taxon.lower() == "bilophila unclassified":
            profile_taxon = "s__Bilophila_wadsworthia"

        if profile_taxon not in taxa_labels:
            vmh_taxon = str(r.get("vmh_taxon_page", "")).strip()
            fallback = to_profile_label(vmh_taxon.split(" AGORA")[0].split(" DSM")[0].split(" ATCC")[0])
            if fallback in taxa_labels:
                profile_taxon = fallback

        rows.append(
            {
                "table2_taxon": table_taxon,
                "profile_taxon": profile_taxon,
                "model_file": model_file,
                "model_file_abs": str(abs_model_file),
                "model_species_id": sanitize(model_stem),
            }
        )
    return pd.DataFrame(rows)


def main():
    # Store all intermediate MICOM inputs together.
    MICOM_DIR.mkdir(parents=True, exist_ok=True)

    for p in [META_XLSX, TAXON_CSV, META_MANIFEST]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required input: {p}")

    meta = pd.read_excel(META_XLSX, sheet_name="Metadata")
    tax = pd.read_csv(TAXON_CSV)
    manifest = pd.read_csv(META_MANIFEST)

    subject_col = first_existing_column(meta, ["Subject ID", "subject_id", "SampleID", "sample_id"])
    age_col = first_existing_column(meta, ["Age (in years)", "Age", "age"])

    taxon_col = tax.columns[0]
    taxa_labels = set(tax[taxon_col].astype(str))

    mapping = build_mapping(manifest, taxa_labels)
    # Mapping table is useful for debugging any missing taxa/models.
    mapping.to_csv(OUT_MAPPING, index=False)

    meta = meta.copy()
    meta[subject_col] = meta[subject_col].astype(str)
    meta["age_group"] = meta[age_col].apply(age_bin)
    meta = meta[meta["age_group"].notna()].copy()
    meta.to_csv(OUT_AGE_SUBJECTS, index=False)

    sample_cols = set(map(str, tax.columns[1:]))
    matched_subjects = [sid for sid in meta[subject_col].tolist() if sid in sample_cols]
    if not matched_subjects:
        raise RuntimeError("No metadata subject IDs matched taxonomy profile sample columns.")

    prof_taxa = set(mapping["profile_taxon"])
    sub = tax[tax[taxon_col].astype(str).isin(prof_taxa)].copy()
    sub = sub[[taxon_col] + matched_subjects]

    rows = []
    summary_rows = []

    for grp in ["21_40", "41_60", "61_70", "71_80", "81_plus"]:
        grp_subjects = [
            sid
            for sid in meta.loc[meta["age_group"] == grp, subject_col].astype(str).tolist()
            if sid in sample_cols
        ]
        if not grp_subjects:
            continue

        med = (
            # Median abundance by taxon within each age group.
            sub.set_index(taxon_col)[grp_subjects]
            .median(axis=1)
            .rename("abundance")
            .reset_index()
            .rename(columns={taxon_col: "profile_taxon"})
        )
        out = mapping.merge(med, on="profile_taxon", how="left").fillna({"abundance": 0.0})
        out["abundance"] = normalize(out["abundance"])
        # MICOM "sample_id" is community label; here we use age-group label.
        out["sample_id"] = grp
        # MICOM expects "id" (organism ID) and "file" (model path).
        out["id"] = out["model_species_id"]
        out["file"] = out["model_file_abs"]

        keep = [
            "sample_id",
            "id",
            "abundance",
            "file",
            "table2_taxon",
            "profile_taxon",
            "model_file",
            "model_species_id",
        ]
        rows.append(out[keep])

        summary_rows.append(
            {
                "age_group": grp,
                "n_subjects": len(grp_subjects),
                "n_taxa_mapped": int((out["abundance"] > 0).sum()),
                "total_abundance_before_norm": float(med["abundance"].sum()),
            }
        )

    if not rows:
        raise RuntimeError("No age groups with matched subjects were found.")

    age_tax = pd.concat(rows, ignore_index=True)
    # Main MICOM input table used by Script 14.
    age_tax.to_csv(OUT_AGE_MEDIAN, index=False)
    pd.DataFrame(summary_rows).to_csv(OUT_AGE_SUMMARY, index=False)

    print(f"Saved: {OUT_MAPPING}")
    print(f"Saved: {OUT_AGE_SUBJECTS} (rows={len(meta)})")
    print(f"Saved: {OUT_AGE_MEDIAN} (rows={len(age_tax)})")
    print(f"Saved: {OUT_AGE_SUMMARY}")


if __name__ == "__main__":
    main()
