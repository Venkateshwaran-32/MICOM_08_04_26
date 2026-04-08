from pathlib import Path
import importlib.util
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_06 = PROJECT_ROOT / "Scripts" / "06_flux_diagnostics.py"
OUT_DIR = PROJECT_ROOT / "Results" / "fba"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT_DIR / "community_model_summary.csv"


def load_script_06_module():
    spec = importlib.util.spec_from_file_location("flux_diagnostics", SCRIPT_06)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    mod = load_script_06_module()

    western = mod.load_medium(mod.WESTERN_CSV)
    fiber = mod.load_medium(mod.FIBER_CSV)
    model_files = sorted(mod.MODELS_DIR.glob("*.xml"))
    medium_ids = set(western.keys()) | set(fiber.keys())

    community, connector_meta, species_ids, biomass_rxn_ids = mod.build_community_model(model_files, medium_ids)

    community_exchange_count = sum(
        1 for rxn in community.reactions if rxn.id in medium_ids
    )
    species_exchange_count = 0
    species_internal_count = 0
    connector_count = 0

    for rxn in community.reactions:
        if rxn.id.startswith("COMM__"):
            connector_count += 1
        elif rxn.id in medium_ids:
            continue
        elif "__EX_" in rxn.id:
            species_exchange_count += 1
        else:
            species_internal_count += 1

    summary = pd.DataFrame(
        [
            {
                "n_species_models": len(model_files),
                "n_species_ids": len(species_ids),
                "n_reactions_total": len(community.reactions),
                "n_metabolites_total": len(community.metabolites),
                "n_genes_total": len(community.genes),
                "n_biomass_reactions": len(biomass_rxn_ids),
                "n_connector_reactions": connector_count,
                "n_community_exchange_reactions": community_exchange_count,
                "n_species_exchange_reactions": species_exchange_count,
                "n_species_internal_reactions": species_internal_count,
                "n_unique_medium_exchange_ids": len(medium_ids),
            }
        ]
    )

    summary.to_csv(OUT_CSV, index=False)

    print(f"Saved: {OUT_CSV}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
