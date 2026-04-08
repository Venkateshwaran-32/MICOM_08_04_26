from pathlib import Path
import pandas as pd

try:
    from micom import Community
except Exception as e:
    raise RuntimeError(
        "MICOM is required. Install it in your project environment first (pip install micom)."
    ) from e


# -------------------------------------------------------------------
# Script 14b: tradeoff sensitivity analysis for MICOM communities
# Inputs:
# - Data/Processed/micom/agegroup_median_taxonomy_for_micom.csv
# - Media/western.csv
# - Media/high_fiber.csv
# Outputs:
# - Results/micom/tradeoff_sensitivity/tradeoff_sensitivity_summary.csv
# - Results/micom/tradeoff_sensitivity/tradeoff_sensitivity_species_growth.csv
# Run:
# - HOME="$PWD" python Scripts/14b_tradeoff_sensitivity_micom.py
# Expected runtime:
# - ~10 to 60 minutes depending on machine / solver
# -------------------------------------------------------------------

# This script evaluates how MICOM results change across cooperative-tradeoff
# fractions. It is useful for checking whether the default 0.5 gives a
# reasonable balance between community growth and the number of taxa growing.

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MICOM_INPUT = PROJECT_ROOT / "Data" / "Processed" / "micom" / "agegroup_median_taxonomy_for_micom.csv"
MEDIA_DIR = PROJECT_ROOT / "Media"

WESTERN_CSV = MEDIA_DIR / "western.csv"
FIBER_CSV = MEDIA_DIR / "high_fiber.csv"

OUT_DIR = PROJECT_ROOT / "Results" / "micom" / "tradeoff_sensitivity" / "proper_age_bins"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_SUMMARY = OUT_DIR / "tradeoff_sensitivity_summary.csv"
OUT_SPECIES = OUT_DIR / "tradeoff_sensitivity_species_growth.csv"

MIN_GROWTH = 1e-6
TRADEOFF_FRACTIONS = [round(x, 1) for x in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]]


def load_medium(csv_path: Path) -> dict[str, float]:
    df = pd.read_csv(csv_path)
    if not {"exchange_id", "max_uptake"}.issubset(df.columns):
        raise ValueError(f"{csv_path.name} must contain columns exchange_id,max_uptake")
    return dict(zip(df["exchange_id"].astype(str), df["max_uptake"].astype(float)))


def to_micom_exchange_id(exchange_id: str) -> str:
    rid = str(exchange_id)
    rid = rid.replace("(e)", "_m")
    rid = rid.replace("[e]", "_m")
    if rid.endswith("_e"):
        rid = rid[:-2] + "_m"
    return rid


def build_community(tax_sub: pd.DataFrame, label: str) -> Community:
    taxonomy = tax_sub[["id", "abundance", "file"]].copy()
    return Community(taxonomy=taxonomy, id=label)


def apply_medium_to_community(com: Community, medium: dict[str, float]) -> tuple[dict[str, float], int]:
    available = {r.id for r in com.exchanges}
    mapped = {}
    for rid, val in medium.items():
        for cand in (rid, to_micom_exchange_id(rid)):
            if cand in available:
                mapped[cand] = float(val)
                break

    if not mapped:
        example = sorted(list(available))[:10]
        raise ValueError(
            "No medium IDs matched MICOM exchange reactions. "
            f"Example MICOM exchange IDs: {example}"
        )

    com.medium = mapped
    return mapped, len(medium) - len(mapped)


def safe_members_df(sol) -> pd.DataFrame:
    mem = sol.members.copy().reset_index()
    if "compartments" in mem.columns:
        mem = mem.rename(columns={"compartments": "id"})
    if "taxon" in mem.columns and "id" not in mem.columns:
        mem = mem.rename(columns={"taxon": "id"})

    growth_col = None
    for c in ["growth_rate", "growth", "mu"]:
        if c in mem.columns:
            growth_col = c
            break
    if growth_col is None:
        raise ValueError("Could not find growth-rate column in MICOM members table.")
    mem = mem.rename(columns={growth_col: "growth_rate"})

    if "id" not in mem.columns:
        raise ValueError("Could not find organism ID column in MICOM members table.")

    return mem


def run_tradeoff(com: Community, fraction: float, min_growth: float):
    return com.cooperative_tradeoff(
        fraction=fraction,
        min_growth=min_growth,
        fluxes=False,
        pfba=False,
    )


def main():
    if not MICOM_INPUT.exists():
        raise FileNotFoundError(f"Missing MICOM input table: {MICOM_INPUT}. Run Script 13 first.")

    tax = pd.read_csv(MICOM_INPUT)
    required = {"sample_id", "id", "abundance", "file"}
    if not required.issubset(tax.columns):
        raise ValueError(f"MICOM input must contain columns: {sorted(required)}")

    western = load_medium(WESTERN_CSV)
    high_fiber = load_medium(FIBER_CSV)

    summary_rows = []
    species_rows = []

    for age_group in sorted(tax["sample_id"].astype(str).unique()):
        tax_sub = tax[tax["sample_id"] == age_group].copy()
        org_ids = sorted(tax_sub["id"].astype(str).unique())

        for diet_name, medium in [("western", western), ("high_fiber", high_fiber)]:
            for fraction in TRADEOFF_FRACTIONS:
                m = build_community(tax_sub, label=f"{age_group}_{diet_name}_tradeoff_{fraction}")
                applied, missing = apply_medium_to_community(m, medium)

                status = "optimal"
                note = ""
                try:
                    sol = run_tradeoff(m, fraction=fraction, min_growth=MIN_GROWTH)
                except Exception:
                    sol = run_tradeoff(m, fraction=fraction, min_growth=0.0)
                    note = "strict_min_growth_infeasible_reran_with_min_growth_0"

                mem = safe_members_df(sol)
                mem = mem[mem["id"].isin(org_ids)].copy()
                mem["age_group"] = age_group
                mem["diet"] = diet_name
                mem["tradeoff_fraction"] = fraction
                mem["grows_in_diet"] = mem["growth_rate"] > MIN_GROWTH

                species_rows.append(
                    mem[["tradeoff_fraction", "age_group", "diet", "id", "growth_rate", "grows_in_diet"]]
                )

                summary_rows.append(
                    {
                        "tradeoff_fraction": fraction,
                        "age_group": age_group,
                        "diet": diet_name,
                        "status": status,
                        "community_growth_rate": float(sol.growth_rate),
                        "n_organisms": len(org_ids),
                        "n_growing_in_diet": int(mem["grows_in_diet"].sum()),
                        "fraction_growing_in_diet": float(mem["grows_in_diet"].sum()) / float(len(org_ids)),
                        "applied_exchanges": len(applied),
                        "missing_exchanges": missing,
                        "note": note,
                    }
                )
                print(
                    f"{age_group} | {diet_name} | tradeoff={fraction:.1f} | "
                    f"growth={float(sol.growth_rate):.6f} | "
                    f"growing={int(mem['grows_in_diet'].sum())}/{len(org_ids)}"
                )

    df_summary = pd.DataFrame(summary_rows).sort_values(
        ["age_group", "diet", "tradeoff_fraction"], ascending=[True, True, True]
    )
    df_species = pd.concat(species_rows, ignore_index=True).sort_values(
        ["tradeoff_fraction", "age_group", "diet", "growth_rate"], ascending=[True, True, True, False]
    )

    df_summary.to_csv(OUT_SUMMARY, index=False)
    df_species.to_csv(OUT_SPECIES, index=False)

    print(f"Saved: {OUT_SUMMARY}")
    print(f"Saved: {OUT_SPECIES}")


if __name__ == "__main__":
    main()
