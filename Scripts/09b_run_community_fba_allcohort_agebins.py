from pathlib import Path
import importlib.util
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = PROJECT_ROOT / "Scripts"
INPUT_DIR = PROJECT_ROOT / "Results" / "inputs"
OUT_DIR = PROJECT_ROOT / "Results" / "fba"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SCENARIO_DIR = OUT_DIR / "scenarios"
SCENARIO_DIR.mkdir(parents=True, exist_ok=True)

ALLCOHORT_INPUT = INPUT_DIR / "allcohort_agebin_input_for_fba.csv"

OUT_SUMMARY = OUT_DIR / "allcohort_agebin_weighting_scenarios_summary.csv"
OUT_SPECIES = OUT_DIR / "allcohort_agebin_weighting_scenarios_species_biomass.csv"
OUT_EXCHANGE = OUT_DIR / "allcohort_agebin_weighting_scenarios_exchange_fluxes.csv"


def load_flux6_module():
    p = SCRIPTS_DIR / "06_flux_diagnostics.py"
    spec = importlib.util.spec_from_file_location("flux6", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def objective_for_weights(model, biomass_rxn_ids: set[str], species_ids: list[str], weights: dict[str, float], helper_mod):
    obj = {}
    for rid in biomass_rxn_ids:
        if rid not in model.reactions:
            continue
        species, _ = helper_mod.infer_species_and_local_id(rid, species_ids)
        w = float(weights.get(species, 0.0))
        if w > 0:
            obj[model.reactions.get_by_id(rid)] = w
    return obj


def run_one(mod, community, species_ids, biomass_rxn_ids, diet_name: str, medium: dict[str, float], scenario_name: str, weights: dict[str, float]):
    m = community.copy()
    m.objective = objective_for_weights(m, biomass_rxn_ids, species_ids, weights, mod)

    applied = mod.apply_medium_bounds(m, medium)
    for rid in ["EX_o2(e)", "EX_o2_e", "EX_o2"]:
        if rid in m.reactions:
            m.reactions.get_by_id(rid).lower_bound = 0.0

    sol = m.optimize()

    summary = {
        "scenario": scenario_name,
        "diet": diet_name,
        "status": sol.status,
        "community_objective": float(sol.objective_value) if sol.status == "optimal" else None,
        "applied_exchanges": len(applied),
    }

    species_rows = []
    if sol.status == "optimal":
        for rid in sorted(biomass_rxn_ids):
            if rid not in m.reactions:
                continue
            species, _ = mod.infer_species_and_local_id(rid, species_ids)
            species_rows.append(
                {
                    "scenario": scenario_name,
                    "diet": diet_name,
                    "species": species,
                    "biomass_reaction": rid,
                    "biomass_flux": float(sol.fluxes[rid]),
                    "weight": float(weights.get(species, 0.0)),
                }
            )

    exchange_rows = []
    if sol.status == "optimal":
        for ex_id in sorted(applied.keys()):
            f = float(sol.fluxes[ex_id])
            exchange_rows.append(
                {
                    "scenario": scenario_name,
                    "diet": diet_name,
                    "exchange_id": ex_id,
                    "flux": f,
                    "uptake_flux": max(0.0, -f),
                    "secretion_flux": max(0.0, f),
                    "abs_flux": abs(f),
                }
            )
    return summary, species_rows, exchange_rows


def save_if_nonempty(df: pd.DataFrame, path: Path):
    if df.empty:
        return
    df.to_csv(path, index=False)
    print(f"Saved: {path}")


def main():
    if not ALLCOHORT_INPUT.exists():
        raise FileNotFoundError(f"Missing required input: {ALLCOHORT_INPUT}. Run script 08b first.")

    mod = load_flux6_module()
    western = mod.load_medium(mod.WESTERN_CSV)
    fiber = mod.load_medium(mod.FIBER_CSV)
    model_files = sorted(mod.MODELS_DIR.glob("*.xml"))
    medium_ids = set(western.keys()) | set(fiber.keys())
    community, _, species_ids, biomass_rxn_ids = mod.build_community_model(model_files, medium_ids)

    age_df = pd.read_csv(ALLCOHORT_INPUT)
    need = {"age_group", "model_species_id", "normalized_weight"}
    if not need.issubset(age_df.columns):
        raise ValueError(f"Input missing required columns: {sorted(need)}")

    scenarios = []
    equal_weights = {sp: 1.0 / len(species_ids) for sp in species_ids}
    scenarios.append(("allcohort_equal_abundance", equal_weights))

    for grp, g in age_df.groupby("age_group"):
        w = dict(zip(g["model_species_id"].astype(str), g["normalized_weight"].astype(float)))
        scenarios.append((f"allcohort_agebin_{grp}", w))

    summary_all = []
    species_all = []
    exchange_all = []

    for scenario_name, weights in scenarios:
        for diet_name, medium in [("western", western), ("high_fiber", fiber)]:
            summary, species_rows, exchange_rows = run_one(
                mod, community, species_ids, biomass_rxn_ids, diet_name, medium, scenario_name, weights
            )
            summary_all.append(summary)
            species_all.extend(species_rows)
            exchange_all.extend(exchange_rows)
            print(f"{scenario_name} | {diet_name} | {summary['status']} | objective={summary['community_objective']}")

    df_summary = pd.DataFrame(summary_all)
    df_species = pd.DataFrame(species_all).sort_values(["scenario", "diet", "biomass_flux"], ascending=[True, True, False])
    df_exchange = pd.DataFrame(exchange_all).sort_values(["scenario", "diet", "abs_flux"], ascending=[True, True, False])

    for scenario_name, _ in scenarios:
        save_if_nonempty(
            df_summary[df_summary["scenario"] == scenario_name],
            SCENARIO_DIR / f"{scenario_name}_community_summary_by_diet.csv",
        )
        save_if_nonempty(
            df_species[df_species["scenario"] == scenario_name],
            SCENARIO_DIR / f"{scenario_name}_species_biomass_by_diet.csv",
        )
        save_if_nonempty(
            df_exchange[df_exchange["scenario"] == scenario_name],
            SCENARIO_DIR / f"{scenario_name}_community_exchange_fluxes_by_diet.csv",
        )

    df_summary.to_csv(OUT_SUMMARY, index=False)
    df_species.to_csv(OUT_SPECIES, index=False)
    df_exchange.to_csv(OUT_EXCHANGE, index=False)

    print(f"Saved: {OUT_SUMMARY}")
    print(f"Saved: {OUT_SPECIES}")
    print(f"Saved: {OUT_EXCHANGE}")
    print(f"Saved scenario files in: {SCENARIO_DIR}")


if __name__ == "__main__":
    main()
