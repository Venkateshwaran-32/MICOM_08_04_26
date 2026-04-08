from pathlib import Path
import pandas as pd
import cobra

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Folder containing your SBML files
MODELS_DIR = PROJECT_ROOT / "Models" / "vmh_agora_sbml"

# Where to save growth results
OUT_CSV = PROJECT_ROOT / "Results" / "fba" / "single_species_growth_default.csv"
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)


def block_oxygen_uptake(model):
    """
    Gut context is mostly anaerobic.
    This function tries common oxygen exchange IDs and blocks oxygen uptake (sets lower_bound = 0).
    If the model uses a different oxygen exchange ID, nothing breaks, it just will not change anything.
    """
    for rid in ["EX_o2(e)", "EX_o2_e", "EX_o2"]:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = 0.0


def main():
    rows = []

    # Loop through every SBML file in your folder
    for fp in sorted(MODELS_DIR.glob("*.xml")):
        print("Running FBA for:", fp.name)

        # Load the model
        model = cobra.io.read_sbml_model(str(fp))

        # (Optional but recommended) block oxygen uptake for gut microbes
        block_oxygen_uptake(model)

        # FAST: only compute objective value (growth), not all fluxes
        # If you later need fluxes, use: sol = model.optimize()
        try:
            growth = model.slim_optimize()
            status = "optimal"
        except Exception as e:
            growth = None
            status = f"error: {type(e).__name__}"

        rows.append(
            {
                "file": fp.name,
                "status": status,
                "growth_rate": growth,
                "objective": str(model.objective.expression),
            }
        )

    # Save to CSV
    df = pd.DataFrame(rows).sort_values(by="growth_rate", ascending=False, na_position="last")
    df.to_csv(OUT_CSV, index=False)

    print("\nSaved growth table to:", OUT_CSV)
    print(df)


if __name__ == "__main__":
    main()
