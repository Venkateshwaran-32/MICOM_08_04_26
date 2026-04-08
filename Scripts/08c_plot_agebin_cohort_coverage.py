from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SUPP_XLSX = PROJECT_ROOT / "Data" / "Supplementary" / "supplementary_data_file_1.xlsx"
TAX_CSV = PROJECT_ROOT / "Data" / "Supplementary" / "taxonomic_profiles_filtered.csv"

OUT_PROC = PROJECT_ROOT / "Data" / "Processed"
OUT_FIG = PROJECT_ROOT / "Results" / "figures"
OUT_PROC.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

OUT_COUNTS = OUT_PROC / "allcohort_agebin_cohort_sample_counts_used.csv"
OUT_PCT = OUT_PROC / "allcohort_agebin_cohort_sample_percent_used.csv"
OUT_PLOT = OUT_FIG / "allcohort_agebin_cohort_sample_counts_used.png"

AGE_ORDER = ["21_40", "41_60", "61_70", "71_80", "81_90"]


def age_bin_21_90(age: float):
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


def first_existing_column(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"Missing expected columns. Tried: {candidates}")


def main():
    if not SUPP_XLSX.exists() or not TAX_CSV.exists():
        raise FileNotFoundError("Missing supplementary metadata or taxonomy file.")

    meta = pd.read_excel(SUPP_XLSX, sheet_name="Metadata")
    tax = pd.read_csv(TAX_CSV)

    subject_col = first_existing_column(meta, ["Subject ID", "subject_id", "SampleID", "sample_id"])
    cohort_col = first_existing_column(meta, ["Cohort", "cohort"])
    age_col = first_existing_column(meta, ["Age (in years)", "Age", "age"])

    meta = meta.copy()
    meta[subject_col] = meta[subject_col].astype(str)
    meta["age_group"] = meta[age_col].apply(age_bin_21_90)
    meta = meta[meta["age_group"].notna()].copy()

    # "Used" samples = subjects present in taxonomy abundance table columns.
    used_samples = set(map(str, tax.columns[1:]))
    meta_used = meta[meta[subject_col].isin(used_samples)].copy()

    counts = (
        meta_used.groupby(["age_group", cohort_col], as_index=False)
        .size()
        .rename(columns={"size": "n_samples_used", cohort_col: "cohort"})
    )
    counts["age_group"] = pd.Categorical(counts["age_group"], categories=AGE_ORDER, ordered=True)
    counts = counts.sort_values(["age_group", "cohort"])
    counts.to_csv(OUT_COUNTS, index=False)

    totals = counts.groupby("age_group", as_index=False)["n_samples_used"].sum().rename(columns={"n_samples_used": "age_total"})
    pct = counts.merge(totals, on="age_group", how="left")
    pct["pct_within_age_group"] = 100.0 * pct["n_samples_used"] / pct["age_total"]
    pct.to_csv(OUT_PCT, index=False)

    pivot = counts.pivot(index="age_group", columns="cohort", values="n_samples_used").fillna(0).reindex(AGE_ORDER).fillna(0)

    ax = pivot.plot(kind="bar", stacked=True, figsize=(10, 6), width=0.75)
    ax.set_title("Samples Used per Age Bin by Cohort")
    ax.set_xlabel("age bin")
    ax.set_ylabel("number of samples used")
    ax.tick_params(axis="x", rotation=0)
    ax.legend(title="cohort", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")

    # Add per-segment labels like "T2D:54" so each cohort contribution is explicit.
    for i, age in enumerate(pivot.index):
        running = 0.0
        for cohort in pivot.columns:
            val = float(pivot.loc[age, cohort])
            if val <= 0:
                continue
            y = running + (val / 2.0)
            ax.text(
                i,
                y,
                f"{cohort}:{int(val)}",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                bbox=dict(facecolor="black", alpha=0.25, edgecolor="none", pad=1.5),
            )
            running += val

    # Add total n above each bar.
    totals_by_age = pivot.sum(axis=1)
    for i, total in enumerate(totals_by_age):
        ax.text(i, total + max(1, pivot.values.max() * 0.02), f"n={int(total)}", ha="center", va="bottom", fontsize=9)

    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=220)
    plt.close()

    print(f"Saved: {OUT_COUNTS}")
    print(f"Saved: {OUT_PCT}")
    print(f"Saved: {OUT_PLOT}")


if __name__ == "__main__":
    main()
