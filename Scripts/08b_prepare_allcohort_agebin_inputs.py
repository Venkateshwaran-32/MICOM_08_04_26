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

OUT_MAPPING = PROC_DIR / "taxon_to_sbml_mapping_allcohort.csv"
OUT_SUBJECTS = PROC_DIR / "allcohort_subjects_age_21_90.csv"
OUT_MEDIAN = PROC_DIR / "allcohort_agebin_median_abundance_by_taxon.csv"
OUT_SUMMARY = PROC_DIR / "allcohort_agebin_input_summary.csv"
OUT_INPUT = INPUT_DIR / "allcohort_agebin_input_for_fba.csv"


def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", text)


def to_profile_label(taxon: str) -> str:
    return "s__" + taxon.replace(" ", "_")


def normalize_weights(df: pd.DataFrame, value_col: str) -> pd.Series:
    total = float(df[value_col].sum())
    if total <= 0:
        return pd.Series([1.0 / len(df)] * len(df), index=df.index)
    return df[value_col] / total


def age_bin_21_90(age: float) -> str | None:
    if pd.isna(age):
        return None
    a = float(age)
    if a < 21 or a > 90:
        return None
    if a <= 40:
        return "21_40"
    if a <= 60:
        return "41_60"
    if a <= 70:
        return "61_70"
    if a <= 80:
        return "71_80"
    return "81_90"


def first_existing_column(df: pd.DataFrame, candidates: list[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"None of the expected columns found: {candidates}")


def build_mapping(manifest: pd.DataFrame, taxa_labels: set[str]) -> pd.DataFrame:
    rows = []
    for _, r in manifest.iterrows():
        table_taxon = str(r["table2_taxon"]).strip()
        model_file = str(r["local_file"]).strip()
        model_stem = Path(model_file).stem

        profile_taxon = to_profile_label(table_taxon)
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
                "model_species_id": sanitize(model_stem),
            }
        )
    return pd.DataFrame(rows)


def main():
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
    mapping.to_csv(OUT_MAPPING, index=False)

    meta = meta.copy()
    meta[subject_col] = meta[subject_col].astype(str)
    meta["age_group"] = meta[age_col].apply(age_bin_21_90)
    meta = meta[meta["age_group"].notna()].copy()
    meta.to_csv(OUT_SUBJECTS, index=False)

    sample_cols = set(tax.columns[1:].astype(str))
    subjects = [sid for sid in meta[subject_col].tolist() if sid in sample_cols]
    if not subjects:
        raise RuntimeError("No subjects in metadata matched taxonomy profile columns.")

    prof_taxa = set(mapping["profile_taxon"])
    sub = tax[tax[taxon_col].astype(str).isin(prof_taxa)].copy()
    sub = sub[[taxon_col] + subjects]

    age_rows = []
    summary_rows = []
    for grp in ["21_40", "41_60", "61_70", "71_80", "81_90"]:
        grp_subjects = [
            sid
            for sid in meta.loc[meta["age_group"] == grp, subject_col].tolist()
            if sid in sample_cols
        ]
        if not grp_subjects:
            continue

        gmed = sub.set_index(taxon_col)[grp_subjects].median(axis=1).rename("median_abundance").reset_index()
        gmed = gmed.rename(columns={taxon_col: "profile_taxon"})
        gout = mapping.merge(gmed, on="profile_taxon", how="left").fillna({"median_abundance": 0.0})
        gout["age_group"] = grp
        gout["n_subjects"] = len(grp_subjects)
        gout["normalized_weight"] = normalize_weights(gout, "median_abundance")
        age_rows.append(gout)

        summary_rows.append(
            {
                "age_group": grp,
                "n_subjects": len(grp_subjects),
                "n_taxa_mapped_nonzero": int((gout["median_abundance"] > 0).sum()),
                "total_median_abundance_pre_norm": float(gout["median_abundance"].sum()),
            }
        )

    if not age_rows:
        raise RuntimeError("No valid age bins (21-90) found with matching taxonomy subjects.")

    age_out = pd.concat(age_rows, ignore_index=True)
    age_out.to_csv(OUT_MEDIAN, index=False)
    age_out.to_csv(OUT_INPUT, index=False)
    pd.DataFrame(summary_rows).to_csv(OUT_SUMMARY, index=False)

    print(f"Saved: {OUT_MAPPING}")
    print(f"Saved: {OUT_SUBJECTS} (rows={len(meta)})")
    print(f"Saved: {OUT_MEDIAN} (rows={len(age_out)})")
    print(f"Saved: {OUT_INPUT}")
    print(f"Saved: {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
