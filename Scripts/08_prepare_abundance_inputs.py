from pathlib import Path
import re
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SUPP_DIR = PROJECT_ROOT / "Data" / "Supplementary"
PROC_DIR = PROJECT_ROOT / "Data" / "Processed"
INPUT_DIR = PROJECT_ROOT / "Results" / "inputs"
META_MANIFEST = PROJECT_ROOT / "Metadata" / "models_manifest.csv" / "models_manifest.csv"

PROC_DIR.mkdir(parents=True, exist_ok=True)
INPUT_DIR.mkdir(parents=True, exist_ok=True)

META_XLSX = SUPP_DIR / "supplementary_data_file_1.xlsx"
TAXON_CSV = SUPP_DIR / "taxonomic_profiles_filtered.csv"

OUT_SG90_SUBJECTS = PROC_DIR / "sg90_subjects.csv"
OUT_MAPPING = PROC_DIR / "taxon_to_sbml_mapping.csv"
OUT_SG90_MEDIAN = PROC_DIR / "sg90_median_abundance_by_taxon.csv"
OUT_AGE_MEDIAN = PROC_DIR / "sg90_median_abundance_by_agegroup.csv"
OUT_SG90_INPUT = INPUT_DIR / "sg90_median_input_for_fba.csv"
OUT_AGE_INPUT = INPUT_DIR / "agegroup_median_input_for_fba.csv"


def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", text)


def to_profile_label(taxon: str) -> str:
    return "s__" + taxon.replace(" ", "_")


def normalize_weights(df: pd.DataFrame, value_col: str) -> pd.Series:
    total = float(df[value_col].sum())
    if total <= 0:
        return pd.Series([1.0 / len(df)] * len(df), index=df.index)
    return df[value_col] / total


def age_bin(age: float) -> str | None:
    if pd.isna(age):
        return None
    if 20 <= age < 40:
        return "20_40"
    if 40 <= age < 60:
        return "40_60"
    if age >= 60:
        return "60_plus"
    return None


def build_mapping(manifest: pd.DataFrame, taxa_labels: set[str]) -> pd.DataFrame:
    rows = []
    for _, r in manifest.iterrows():
        table_taxon = str(r["table2_taxon"]).strip()
        model_file = str(r["local_file"]).strip()
        model_stem = Path(model_file).stem
        species_id = sanitize(model_stem)

        # Proxy handling: Bilophila unclassified -> Bilophila wadsworthia
        profile_taxon = to_profile_label(table_taxon)
        if table_taxon.lower() == "bilophila unclassified":
            profile_taxon = "s__Bilophila_wadsworthia"

        if profile_taxon not in taxa_labels:
            # Fallback for edge cases: try with taxon page if present
            vmh_taxon = str(r.get("vmh_taxon_page", "")).strip()
            fallback = to_profile_label(vmh_taxon.split(" AGORA")[0].split(" DSM")[0].split(" ATCC")[0])
            if fallback in taxa_labels:
                profile_taxon = fallback

        rows.append(
            {
                "table2_taxon": table_taxon,
                "profile_taxon": profile_taxon,
                "model_file": model_file,
                "model_species_id": species_id,
            }
        )
    return pd.DataFrame(rows)


def main():
    if not META_XLSX.exists():
        raise FileNotFoundError(f"Missing metadata file: {META_XLSX}")
    if not TAXON_CSV.exists():
        raise FileNotFoundError(f"Missing taxonomy file: {TAXON_CSV}")
    if not META_MANIFEST.exists():
        raise FileNotFoundError(f"Missing model manifest: {META_MANIFEST}")

    meta = pd.read_excel(META_XLSX, sheet_name="Metadata")
    tax = pd.read_csv(TAXON_CSV)
    manifest = pd.read_csv(META_MANIFEST)

    required_meta_cols = {"Subject ID", "Age (in years)", "Cohort"}
    if not required_meta_cols.issubset(meta.columns):
        raise ValueError(f"Metadata file must include columns: {sorted(required_meta_cols)}")

    taxon_col = tax.columns[0]
    taxa_labels = set(tax[taxon_col].astype(str))

    mapping = build_mapping(manifest, taxa_labels)
    mapping.to_csv(OUT_MAPPING, index=False)

    sg90_meta = meta[meta["Cohort"].astype(str) == "SG90"].copy()
    sg90_meta["Subject ID"] = sg90_meta["Subject ID"].astype(str)
    sg90_meta["age_group"] = sg90_meta["Age (in years)"].apply(age_bin)
    sg90_meta.to_csv(OUT_SG90_SUBJECTS, index=False)

    sample_cols = set(tax.columns[1:].astype(str))
    sg90_subjects = [sid for sid in sg90_meta["Subject ID"].tolist() if sid in sample_cols]
    if not sg90_subjects:
        raise RuntimeError("No SG90 subject IDs matched taxonomy table columns.")

    # Keep only mapped taxa rows and SG90 sample columns.
    prof_taxa = set(mapping["profile_taxon"])
    sub = tax[tax[taxon_col].astype(str).isin(prof_taxa)].copy()
    sub = sub[[taxon_col] + sg90_subjects]

    # Scenario A: median per taxon across all SG90 subjects.
    med = sub.set_index(taxon_col)[sg90_subjects].median(axis=1).rename("median_abundance").reset_index()
    med = med.rename(columns={taxon_col: "profile_taxon"})
    sg90_out = mapping.merge(med, on="profile_taxon", how="left").fillna({"median_abundance": 0.0})
    sg90_out["normalized_weight"] = normalize_weights(sg90_out, "median_abundance")
    sg90_out.to_csv(OUT_SG90_MEDIAN, index=False)
    sg90_out.to_csv(OUT_SG90_INPUT, index=False)

    # Scenario B: median per taxon by age group within SG90.
    age_rows = []
    for grp in ["20_40", "40_60", "60_plus"]:
        grp_subjects = [sid for sid in sg90_meta.loc[sg90_meta["age_group"] == grp, "Subject ID"].tolist() if sid in sample_cols]
        if not grp_subjects:
            continue
        gmed = sub.set_index(taxon_col)[grp_subjects].median(axis=1).rename("median_abundance").reset_index()
        gmed = gmed.rename(columns={taxon_col: "profile_taxon"})
        gout = mapping.merge(gmed, on="profile_taxon", how="left").fillna({"median_abundance": 0.0})
        gout["age_group"] = grp
        gout["n_subjects"] = len(grp_subjects)
        gout["normalized_weight"] = normalize_weights(gout, "median_abundance")
        age_rows.append(gout)

    age_out = pd.concat(age_rows, ignore_index=True) if age_rows else pd.DataFrame()
    age_out.to_csv(OUT_AGE_MEDIAN, index=False)
    age_out.to_csv(OUT_AGE_INPUT, index=False)

    print(f"Saved: {OUT_MAPPING}")
    print(f"Saved: {OUT_SG90_SUBJECTS} (rows={len(sg90_meta)}, matched_tax_samples={len(sg90_subjects)})")
    print(f"Saved: {OUT_SG90_MEDIAN} (rows={len(sg90_out)})")
    print(f"Saved: {OUT_AGE_MEDIAN} (rows={len(age_out)})")
    print(f"Saved: {OUT_SG90_INPUT}")
    print(f"Saved: {OUT_AGE_INPUT}")


if __name__ == "__main__":
    main()
