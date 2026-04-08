from pathlib import Path
import pandas as pd
import cobra

# -------------------------------------------------------------------
# Script 14: run MICOM communities by age group and diet
# Inputs:
# - Data/Processed/micom/agegroup_median_taxonomy_for_micom.csv
# - Media/western.csv
# -"Media"
#  Outputs:
# - Results/micom/growth/organism_growth_rates_by_agegroup_diet.csv
# - Results/micom/growth/community_growth_summary_by_agegroup_diet.csv
# - Results/micom/pathway_flux/reaction_fluxes_long_by_agegroup_diet.csv
# - Results/micom/growth/organism_growth_capability_check.csv
# - Results/micom/pathway_flux/reaction_annotations_by_organism.csv
# Run:
# - .venv\Scripts\python Scripts\14_run_micom_by_agegroup.py
# Expected runtime:
# - ~5 to 20 minutes (depends on solver and machine)
# -------------------------------------------------------------------

# Core MICOM runner:
# 1) builds one community per age group,
# 2) applies diet media,
# 3) runs cooperative tradeoff,
# 4) exports organism growth + reaction fluxes,
# 5) checks "all can grow" under full-access medium.

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

OUT_DIR = PROJECT_ROOT / "Results" / "micom"
GROWTH_DIR = OUT_DIR / "growth" / "proper_age_bins"
PATHWAY_DIR = OUT_DIR / "pathway_flux" / "proper_age_bins"
GROWTH_DIR.mkdir(parents=True, exist_ok=True)
PATHWAY_DIR.mkdir(parents=True, exist_ok=True)

OUT_GROWTH = GROWTH_DIR / "organism_growth_rates_by_agegroup_diet.csv"
OUT_COMM = GROWTH_DIR / "community_growth_summary_by_agegroup_diet.csv"
OUT_FLUX = PATHWAY_DIR / "reaction_fluxes_long_by_agegroup_diet.csv"
OUT_CAP = GROWTH_DIR / "organism_growth_capability_check.csv"
OUT_ANNOT = PATHWAY_DIR / "reaction_annotations_by_organism.csv"

MIN_GROWTH = 1e-6
TRADEOFF_FRACTION = 0.5
FULL_ACCESS_UPTAKE = 1000.0


def load_medium(csv_path: Path) -> dict[str, float]:
    df = pd.read_csv(csv_path)
    if not {"exchange_id", "max_uptake"}.issubset(df.columns):
        raise ValueError(f"{csv_path.name} must contain columns exchange_id,max_uptake")
    return dict(zip(df["exchange_id"].astype(str), df["max_uptake"].astype(float)))


def to_micom_exchange_id(exchange_id: str) -> str:
    # Your media files use EX_xxx(e), while MICOM community exchanges use EX_xxx_m.
    # This helper makes that translation so media constraints apply correctly.
    rid = str(exchange_id)
    rid = rid.replace("(e)", "_m")
    rid = rid.replace("[e]", "_m")
    if rid.endswith("_e"):
        rid = rid[:-2] + "_m"
    return rid


def first_nonempty(value):
    if value is None:
        return None
    if isinstance(value, (list, tuple, set)):
        vals = [str(v).strip() for v in value if str(v).strip()]
        return " | ".join(vals) if vals else None
    s = str(value).strip()
    return s if s else None


def reaction_pathway(rxn: cobra.Reaction) -> str:
    # Pull pathway/subsystem labels from reaction fields/annotations for downstream plots.
    cand = first_nonempty(getattr(rxn, "subsystem", None))
    if cand:
        return cand
    ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
    for key in ["subsystem", "pathway", "kegg.pathway", "metacyc.reaction", "sbo"]:
        cand = first_nonempty(ann.get(key))
        if cand:
            return cand
    return "Unassigned"


def build_reaction_annotation(tax_df: pd.DataFrame) -> pd.DataFrame:
    # Build lookup table: (organism, reaction) -> pathway name.
    # Needed because MICOM flux output itself does not always carry pathway labels.
    rows = []
    unique_models = tax_df[["id", "file"]].drop_duplicates()
    for _, r in unique_models.iterrows():
        org = str(r["id"])
        model_fp = Path(str(r["file"]))
        if not model_fp.exists():
            raise FileNotFoundError(f"Missing model file for {org}: {model_fp}")

        model = cobra.io.read_sbml_model(str(model_fp))
        for rxn in model.reactions:
            rows.append(
                {
                    "id": org,
                    "reaction_id": rxn.id,
                    "reaction_name": rxn.name,
                    "pathway": reaction_pathway(rxn),
                }
            )
    ann = pd.DataFrame(rows)
    ann.to_csv(OUT_ANNOT, index=False)
    return ann


def build_community(tax_sub: pd.DataFrame, label: str) -> Community:
    # MICOM Community is built directly from taxonomy table columns: id, abundance, file.
    taxonomy = tax_sub[["id", "abundance", "file"]].copy()
    return Community(taxonomy=taxonomy, id=label)


def apply_medium_to_community(com: Community, medium: dict[str, float]) -> tuple[dict[str, float], int]:
    available = {r.id for r in com.exchanges}

    mapped = {}
    for rid, val in medium.items():
        # Try exact ID first, then translated MICOM ID.
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
    # Normalize MICOM solution schema across versions.
    # Different versions can name organism or growth columns differently.
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
    # MICOM cooperative tradeoff:
    # maximize community growth subject to a fraction, then distribute individual growth.
    return com.cooperative_tradeoff(
        fraction=TRADEOFF_FRACTION,
        min_growth=min_growth,
        fluxes=True,
        pfba=False,
    )


def run_full_access_capability(tax_sub: pd.DataFrame, organism_ids: list[str]) -> pd.DataFrame:
    # Diagnostic run: open all exchanges to a high uptake bound.
    # If a taxon grows here but not in diet run, limitation is diet/competition, not model incapability.
    m = build_community(tax_sub, label="full_access_check")
    m.medium = {r.id: FULL_ACCESS_UPTAKE for r in m.exchanges}
    sol = run_tradeoff(m, min_growth=0.0)
    mem = safe_members_df(sol)
    mem = mem[mem["id"].isin(organism_ids)].copy()
    mem["grows_full_access"] = mem["growth_rate"] > MIN_GROWTH
    return mem[["id", "growth_rate", "grows_full_access"]].rename(
        columns={"growth_rate": "growth_rate_full_access"}
    )


def solution_flux_long(sol, organism_ids: list[str], age_group: str, diet: str) -> pd.DataFrame:
    # Convert wide MICOM flux matrix into long table for easier filtering/plotting.
    fl = sol.fluxes.copy()
    fl.index.name = "id"
    long = fl.reset_index().melt(id_vars="id", var_name="reaction_id", value_name="flux")
    long = long[long["id"].isin(organism_ids)].copy()
    long["age_group"] = age_group
    long["diet"] = diet
    long["abs_flux"] = long["flux"].abs()
    return long


def main():
    if not MICOM_INPUT.exists():
        raise FileNotFoundError(f"Missing MICOM input table: {MICOM_INPUT}. Run Script 13 first.")

    tax = pd.read_csv(MICOM_INPUT)
    required = {"sample_id", "id", "abundance", "file"}
    if not required.issubset(tax.columns):
        raise ValueError(f"MICOM input must contain columns: {sorted(required)}")

    western = load_medium(WESTERN_CSV)
    high_fiber = load_medium(FIBER_CSV)

    ann = build_reaction_annotation(tax)

    growth_rows = []
    community_rows = []
    flux_rows = []
    capability_rows = []

    for age_group in sorted(tax["sample_id"].astype(str).unique()):
        # One MICOM community per age group.
        tax_sub = tax[tax["sample_id"] == age_group].copy()
        org_ids = sorted(tax_sub["id"].astype(str).unique())

        full_access = run_full_access_capability(tax_sub, org_ids)

        for diet_name, medium in [("western", western), ("high_fiber", high_fiber)]:
            # Rebuild per-run community to avoid carry-over state between diets.
            m = build_community(tax_sub, label=f"{age_group}_{diet_name}")
            applied, missing = apply_medium_to_community(m, medium)

            status = "optimal"
            note = ""
            try:
                # Preferred run with minimum per-organism growth threshold.
                sol = run_tradeoff(m, min_growth=MIN_GROWTH)
            except Exception:
                # Fallback keeps workflow robust if strict min-growth becomes infeasible.
                sol = run_tradeoff(m, min_growth=0.0)
                note = "strict_min_growth_infeasible_reran_with_min_growth_0"

            mem = safe_members_df(sol)
            mem = mem[mem["id"].isin(org_ids)].copy()
            mem["age_group"] = age_group
            mem["diet"] = diet_name
            mem["grows_in_diet"] = mem["growth_rate"] > MIN_GROWTH
            mem = mem.merge(full_access, on="id", how="left")
            mem["nongrowing_due_to_competition_or_diet"] = (~mem["grows_in_diet"]) & (mem["grows_full_access"])

            growth_rows.append(
                mem[
                    [
                        "age_group",
                        "diet",
                        "id",
                        "growth_rate",
                        "grows_in_diet",
                        "growth_rate_full_access",
                        "grows_full_access",
                        "nongrowing_due_to_competition_or_diet",
                    ]
                ]
            )

            n_total = len(mem)
            n_grow_diet = int(mem["grows_in_diet"].sum())
            n_grow_full = int(mem["grows_full_access"].sum())
            community_rows.append(
                {
                    "age_group": age_group,
                    "diet": diet_name,
                    "status": status,
                    "community_growth_rate": float(sol.growth_rate),
                    "n_organisms": n_total,
                    "n_growing_in_diet": n_grow_diet,
                    "n_not_growing_in_diet": n_total - n_grow_diet,
                    "n_growing_full_access": n_grow_full,
                    "applied_exchanges": len(applied),
                    "missing_exchanges": missing,
                    "note": note,
                }
            )

            fl = solution_flux_long(sol, org_ids, age_group, diet_name)
            fl = fl.merge(ann, on=["id", "reaction_id"], how="left")
            flux_rows.append(fl)

        cap = full_access.copy()
        cap["age_group"] = age_group
        capability_rows.append(cap[["age_group", "id", "growth_rate_full_access", "grows_full_access"]])

    df_growth = pd.concat(growth_rows, ignore_index=True).sort_values(
        ["age_group", "diet", "growth_rate"], ascending=[True, True, False]
    )
    df_comm = pd.DataFrame(community_rows).sort_values(["age_group", "diet"])
    df_flux = pd.concat(flux_rows, ignore_index=True).sort_values(
        ["age_group", "diet", "id", "abs_flux"], ascending=[True, True, True, False]
    )
    df_cap = pd.concat(capability_rows, ignore_index=True).drop_duplicates().sort_values(["age_group", "id"])

    df_growth.to_csv(OUT_GROWTH, index=False)
    df_comm.to_csv(OUT_COMM, index=False)
    df_flux.to_csv(OUT_FLUX, index=False)
    df_cap.to_csv(OUT_CAP, index=False)
    # OUT_ANNOT is written inside build_reaction_annotation().

    print(f"Saved: {OUT_GROWTH} ({len(df_growth)} rows)")
    print(f"Saved: {OUT_COMM} ({len(df_comm)} rows)")
    print(f"Saved: {OUT_FLUX} ({len(df_flux)} rows)")
    print(f"Saved: {OUT_CAP} ({len(df_cap)} rows)")
    print(f"Saved: {OUT_ANNOT}")


if __name__ == "__main__":
    main()
