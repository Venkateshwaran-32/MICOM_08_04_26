from pathlib import Path
import pandas as pd
import cobra

# --------- PATHS ----------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

MODELS_DIR = PROJECT_ROOT / "Models" / "vmh_agora_sbml"
MEDIA_DIR = PROJECT_ROOT / "Media"

WESTERN_CSV = MEDIA_DIR / "western.csv"
FIBER_CSV = MEDIA_DIR / "high_fiber.csv"

OUT_DIR = PROJECT_ROOT / "Results" / "fba"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CSV = OUT_DIR / "single_species_growth_by_diet.csv"


def load_medium(csv_path: Path) -> dict:
    """
    Read a medium CSV with columns: exchange_id,max_uptake.
    Return a dict suitable for model.medium.
    """
    df = pd.read_csv(csv_path)
    if not {"exchange_id", "max_uptake"}.issubset(df.columns):
        raise ValueError(f"{csv_path.name} must have columns: exchange_id,max_uptake")
    return dict(zip(df["exchange_id"].astype(str), df["max_uptake"].astype(float)))


def apply_medium(model: cobra.Model, medium: dict):
    """
    Apply a medium dict to the model, but only for exchange reactions that exist in this model.
    """
    available = {rid: val for rid, val in medium.items() if rid in model.reactions}
    model.medium = available
    return len(available), (len(medium) - len(available))


def block_oxygen_uptake(model: cobra.Model):
    """
    Gut context is mostly anaerobic: block oxygen uptake if an oxygen exchange exists.
    """
    for rid in ["EX_o2(e)", "EX_o2_e", "EX_o2"]:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = 0.0


def run_growth(model: cobra.Model):
    """
    Run FBA and return (status, growth_rate).
    """
    sol = model.optimize()
    return sol.status, (sol.objective_value if sol.status == "optimal" else None)


def main():
    western_medium = load_medium(WESTERN_CSV)
    fiber_medium = load_medium(FIBER_CSV)

    print("Loaded media:")
    print(" - Western rows:", len(western_medium))
    print(" - High-fiber rows:", len(fiber_medium))

    rows = []
    model_files = sorted(MODELS_DIR.glob("*.xml"))

    print("\nModels found:", len(model_files))
    for fp in model_files:
        print("\n=== Model:", fp.name, "===")
        model = cobra.io.read_sbml_model(str(fp))

        # Keep gut-like assumption consistent
        block_oxygen_uptake(model)

        # --- Western ---
        n_applied_w, n_missing_w = apply_medium(model, western_medium)
        status_w, growth_w = run_growth(model)
        print(f"Western: applied {n_applied_w} exchanges, missing {n_missing_w}, status={status_w}, growth={growth_w}")

        # Reload model before applying second medium to avoid carry-over constraints
        model = cobra.io.read_sbml_model(str(fp))
        block_oxygen_uptake(model)

        # --- High fiber ---
        n_applied_h, n_missing_h = apply_medium(model, fiber_medium)
        status_h, growth_h = run_growth(model)
        print(f"HighFiber: applied {n_applied_h} exchanges, missing {n_missing_h}, status={status_h}, growth={growth_h}")

        rows.append(
            {
                "file": fp.name,
                "status_western": status_w,
                "growth_western": growth_w,
                "applied_exchanges_western": n_applied_w,
                "missing_exchanges_western": n_missing_w,
                "status_high_fiber": status_h,
                "growth_high_fiber": growth_h,
                "applied_exchanges_high_fiber": n_applied_h,
                "missing_exchanges_high_fiber": n_missing_h,
                "delta_western_minus_fiber": (growth_w - growth_h) if (growth_w is not None and growth_h is not None) else None,
            }
        )

    df = pd.DataFrame(rows)
    df = df.sort_values(by="growth_western", ascending=False, na_position="last")
    df.to_csv(OUT_CSV, index=False)
    print("\nSaved:", OUT_CSV)
    print(df[["file", "growth_western", "growth_high_fiber", "delta_western_minus_fiber"]])


if __name__ == "__main__":
    main()
