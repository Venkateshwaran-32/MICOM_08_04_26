from pathlib import Path

import pandas as pd

try:
    from micom import Community
except Exception as e:
    raise RuntimeError(
        "MICOM is required. Install it in your project environment first (pip install micom)."
    ) from e


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MICOM_INPUT = PROJECT_ROOT / "Data" / "Processed" / "micom" / "agegroup_median_taxonomy_for_micom.csv"
MEDIA_DIR = PROJECT_ROOT / "Media"

WESTERN_CSV = MEDIA_DIR / "western.csv"
FIBER_CSV = MEDIA_DIR / "high_fiber.csv"

OUT_DIR = PROJECT_ROOT / "Results" / "micom" / "growth" / "proper_age_bins"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_SPECIES = OUT_DIR / "scaled_uptake_species_growth_by_agegroup_diet.csv"
OUT_COMMUNITY = OUT_DIR / "scaled_uptake_community_growth_by_agegroup_diet.csv"

MIN_GROWTH = 1e-6
TRADEOFF_FRACTION = 0.5
SCALES = [5, 10]


def load_medium(csv_path: Path) -> dict[str, float]:
    df = pd.read_csv(csv_path)
    if not {"exchange_id", "max_uptake"}.issubset(df.columns):
        raise ValueError(f"{csv_path.name} must contain columns exchange_id,max_uptake")
    return dict(zip(df["exchange_id"].astype(str), df["max_uptake"].astype(float)))


def scale_medium(medium: dict[str, float], factor: int) -> dict[str, float]:
    return {rid: float(val) * float(factor) for rid, val in medium.items()}


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


def run_tradeoff(com: Community, min_growth: float):
    return com.cooperative_tradeoff(
        fraction=TRADEOFF_FRACTION,
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

    species_rows = []
    community_rows = []

    for scale in SCALES:
        scaled_media = {
            "western": scale_medium(western, scale),
            "high_fiber": scale_medium(high_fiber, scale),
        }
        scale_label = f"scaled{scale}x"

        for age_group in sorted(tax["sample_id"].astype(str).unique()):
            tax_sub = tax[tax["sample_id"] == age_group].copy()
            org_ids = sorted(tax_sub["id"].astype(str).unique())

            for diet_name, medium in scaled_media.items():
                m = build_community(tax_sub, label=f"{age_group}_{diet_name}_{scale_label}")
                applied, missing = apply_medium_to_community(m, medium)

                note = ""
                try:
                    sol = run_tradeoff(m, min_growth=MIN_GROWTH)
                except Exception:
                    sol = run_tradeoff(m, min_growth=0.0)
                    note = "strict_min_growth_infeasible_reran_with_min_growth_0"

                mem = safe_members_df(sol)
                mem = mem[mem["id"].isin(org_ids)].copy()
                mem["age_group"] = age_group
                mem["diet"] = diet_name
                mem["uptake_scale"] = scale_label
                mem["grows_in_diet"] = mem["growth_rate"] > MIN_GROWTH
                species_rows.append(mem[["uptake_scale", "age_group", "diet", "id", "growth_rate", "grows_in_diet"]])

                community_rows.append(
                    {
                        "uptake_scale": scale_label,
                        "age_group": age_group,
                        "diet": diet_name,
                        "community_growth_rate": float(sol.growth_rate),
                        "n_organisms": len(org_ids),
                        "n_growing_in_diet": int(mem["grows_in_diet"].sum()),
                        "applied_exchanges": len(applied),
                        "missing_exchanges": missing,
                        "note": note,
                    }
                )
                print(
                    f"{scale_label} | {age_group} | {diet_name} | "
                    f"community_growth={float(sol.growth_rate):.6f} | "
                    f"growing={int(mem['grows_in_diet'].sum())}/{len(org_ids)}"
                )

    df_species = pd.concat(species_rows, ignore_index=True).sort_values(
        ["uptake_scale", "age_group", "diet", "growth_rate"], ascending=[True, True, True, False]
    )
    df_community = pd.DataFrame(community_rows).sort_values(
        ["uptake_scale", "age_group", "diet"], ascending=[True, True, True]
    )

    df_species.to_csv(OUT_SPECIES, index=False)
    df_community.to_csv(OUT_COMMUNITY, index=False)
    print(f"Saved: {OUT_SPECIES}")
    print(f"Saved: {OUT_COMMUNITY}")


if __name__ == "__main__":
    main()
