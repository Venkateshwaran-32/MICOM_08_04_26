from pathlib import Path
import pandas as pd
import cobra

# --------- PATHS ----------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODELS_DIR   = PROJECT_ROOT / "Models" / "vmh_agora_sbml"
MEDIA_DIR    = PROJECT_ROOT / "Media"

WESTERN_CSV = MEDIA_DIR / "western.csv"
FIBER_CSV   = MEDIA_DIR / "high_fiber.csv"

OUT_DIR = PROJECT_ROOT / "Results" / "fba"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CSV = OUT_DIR / "individual_growth_full_access_medium.csv"

# --------- SETTINGS ----------
# "Full access" bound for EVERY candidate exchange (diagnostic: set high).
UPTAKE_BOUND_PER_EXCHANGE = 1.0

# Numerical tolerance for "nonzero growth"
GROWTH_EPS = 1e-9

# --------- HELPERS ----------
def load_exchange_ids(csv_path: Path) -> set[str]:
    """
    Read medium CSV with columns: exchange_id,max_uptake
    Return a set of exchange reaction IDs.
    """
    df = pd.read_csv(csv_path)
    if "exchange_id" not in df.columns:
        raise ValueError(f"{csv_path.name} must contain column: exchange_id")
    return set(df["exchange_id"].astype(str))

def block_oxygen_uptake(model: cobra.Model):
    """
    Gut context is mostly anaerobic: block oxygen uptake if an oxygen exchange exists.
    """
    for rid in ["EX_o2(e)", "EX_o2_e", "EX_o2"]:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = 0.0

def apply_full_access_medium(model: cobra.Model, candidate_ex_ids: set[str], bound: float):
    """
    Apply a standardized medium where all candidate exchanges are opened to the SAME uptake bound.
    Only apply exchanges that actually exist AND are true exchanges in this model.
    Returns (n_applied, n_missing).
    """
    model_exchange_ids = {rxn.id for rxn in model.exchanges}
    available = {rid: float(bound) for rid in candidate_ex_ids if rid in model_exchange_ids}

    # This sets all other exchange uptakes to 0 unless included in 'available'
    model.medium = available

    n_applied = len(available)
    n_missing = len(candidate_ex_ids) - n_applied
    return n_applied, n_missing

def run_growth(model: cobra.Model):
    sol = model.optimize()
    growth = sol.objective_value if sol.status == "optimal" else None
    return sol.status, growth

# --------- MAIN ----------
def main():
    # candidate exchanges = union of both diet lists (your "food/gut nutrient candidates")
    western_ids = load_exchange_ids(WESTERN_CSV)
    fiber_ids   = load_exchange_ids(FIBER_CSV)
    candidate_ex_ids = western_ids | fiber_ids

    print("Candidate exchange IDs (union diets):", len(candidate_ex_ids))

    rows = []
    model_files = sorted(MODELS_DIR.glob("*.xml"))
    print("Models found:", len(model_files))

    for fp in model_files:
        species = fp.stem  # keep it simple + human-readable

        print("Running full-access FBA for:", fp.name)
        model = cobra.io.read_sbml_model(str(fp))

        # keep your gut assumption consistent with the rest of your pipeline
        block_oxygen_uptake(model)

        # apply standardized full-access medium
        applied, missing = apply_full_access_medium(model, candidate_ex_ids, UPTAKE_BOUND_PER_EXCHANGE)

        # run FBA
        try:
            status, growth = run_growth(model)
            notes = ""
        except Exception as e:
            status, growth = f"error: {type(e).__name__}", None
            notes = repr(e)

        # extra interpretability note
        if status == "optimal" and (growth is None or growth <= GROWTH_EPS):
            # optimal but zero growth = feasible, but cannot make biomass under these constraints
            if notes:
                notes += " | "
            notes += "optimal_but_zero_growth"

        rows.append({
            # requested columns
            "species": species,
            "status": status,
            "growth_full_access": growth,
            "applied_exchanges": applied,
            "missing_exchanges": missing,
            "uptake_bound_per_exchange": UPTAKE_BOUND_PER_EXCHANGE,
            "notes": notes,

            # helpful extras (keep or delete if PI wants strict schema)
            "file": fp.name,
            "objective": str(model.objective.expression),
        })

    df = pd.DataFrame(rows).sort_values(by="growth_full_access", ascending=False, na_position="last")
    df.to_csv(OUT_CSV, index=False)

    print("\nSaved:", OUT_CSV)
    print(df[["species", "status", "growth_full_access", "applied_exchanges", "missing_exchanges"]])

if __name__ == "__main__":
    main()