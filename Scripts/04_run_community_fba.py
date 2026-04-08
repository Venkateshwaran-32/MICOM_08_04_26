from pathlib import Path
import re
import pandas as pd
import cobra
from cobra import Model, Reaction, Metabolite


# ---------- PATHS ----------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODELS_DIR = PROJECT_ROOT / "Models" / "vmh_agora_sbml"
MEDIA_DIR = PROJECT_ROOT / "Media"
OUT_DIR = PROJECT_ROOT / "Results" / "fba"
OUT_DIR.mkdir(parents=True, exist_ok=True)

WESTERN_CSV = MEDIA_DIR / "western.csv"
FIBER_CSV = MEDIA_DIR / "high_fiber.csv"

OUT_SUMMARY = OUT_DIR / "community_growth_by_diet.csv"
OUT_SPECIES = OUT_DIR / "community_species_objective_contributions_by_diet.csv"

# the objective for script 04 is to write a community-level summary plus a
# species-level table with biomass fluxes and contributions to the objective.

# ---------- HELPERS ----------
def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", text)


def load_medium(csv_path: Path) -> dict:
    df = pd.read_csv(csv_path)
    if not {"exchange_id", "max_uptake"}.issubset(df.columns):
        raise ValueError(f"{csv_path.name} must contain columns: exchange_id,max_uptake")
    return dict(zip(df["exchange_id"].astype(str), df["max_uptake"].astype(float)))


def get_objective_reactions(model: cobra.Model) -> dict:
    # In COBRApy, objective coefficients are defined on reactions.
    # Using reaction.objective_coefficient is robust across solver backends.
    obj = {
        rxn.id: float(rxn.objective_coefficient)
        for rxn in model.reactions
        if abs(float(rxn.objective_coefficient)) > 1e-12
    }
    if obj:
        return obj

    # fallback if objective parsing fails
    for rxn in model.reactions:
        if "biomass" in rxn.id.lower() or "biomass" in rxn.name.lower():
            return {rxn.id: 1.0}
    raise ValueError(f"Could not infer objective reaction for model {model.id}")


def prefixed_copy(model: cobra.Model, prefix: str) -> cobra.Model:
    m = model.copy()
    for met in list(m.metabolites):
        met.id = f"{prefix}__{met.id}"
    for rxn in list(m.reactions):
        rxn.id = f"{prefix}__{rxn.id}"
    for gene in list(m.genes):
        gene.id = f"{prefix}__{gene.id}"
    return m


def shared_met_id(original_ext_met_id: str) -> str:
    return f"u__{sanitize(original_ext_met_id)}"


def build_community_model(model_files, medium_exchange_ids):
    community = Model("community_model")

    # maps original extracellular metabolite ID -> shared metabolite object
    shared_mets = {}

    # maps exchange rxn ID from diet files -> original extracellular metabolite ID
    exchange_to_ext_met = {}

    # species -> objective rxn IDs (original IDs, pre-prefix)
    species_objectives = {}

    # species -> one representative biomass rxn (for reporting)
    species_main_biomass = {}

    for fp in model_files:
        base = cobra.io.read_sbml_model(str(fp))
        species = sanitize(fp.stem)

        # objective reactions in original model (before prefix)
        obj_rxn_map = get_objective_reactions(base)
        species_objectives[species] = obj_rxn_map
        species_main_biomass[species] = next(iter(obj_rxn_map.keys()))

        # map medium exchange IDs to extracellular metabolite IDs using this model
        for ex_id in medium_exchange_ids:
            if ex_id in base.reactions:
                rxn = base.reactions.get_by_id(ex_id)
                mets = list(rxn.metabolites.keys())
                if len(mets) == 1:
                    ext_met_id = mets[0].id
                    exchange_to_ext_met.setdefault(ex_id, ext_met_id)

        # collect extracellular metabolite IDs in this species
        ext_met_ids = [m.id for m in base.metabolites if m.compartment == "e"]

        # prefix and add species model to community
        pref = prefixed_copy(base, species)
        community.add_reactions(pref.reactions)
        # CROSS FEEDING CONNECTORS: add reversible reactions between species extracellular and shared pool
        # add reversible connector between species extracellular and shared pool
        for ext_id in ext_met_ids:
            local_id = f"{species}__{ext_id}"
            if local_id not in community.metabolites:
                continue

            if ext_id not in shared_mets:
                sm = Metabolite(
                    id=shared_met_id(ext_id),
                    name=f"shared_{ext_id}",
                    compartment="u",
                )
                community.add_metabolites([sm])
                shared_mets[ext_id] = sm

            conn_id = f"COMM__{species}__{sanitize(ext_id)}"
            if conn_id in community.reactions:
                continue
            # connector bounds are -1000 to 1000 because they allow reversible cross-feeding and are not meant to be restrictive
            # - diet bounds are small because they represent limited nutrient supply from the environment
        
            conn = Reaction(conn_id)
            conn.lower_bound = -1000.0
            conn.upper_bound = 1000.0
            conn.add_metabolites(
                {
                    community.metabolites.get_by_id(local_id): -1.0,
                    shared_mets[ext_id]: 1.0,
                }
            )
            community.add_reactions([conn])

    # add community-level exchange reactions with IDs matching diet CSV exchange IDs
    for ex_id, ext_met_id in exchange_to_ext_met.items():
        if ex_id in community.reactions:
            continue
        if ext_met_id not in shared_mets:
            continue

        ex = Reaction(ex_id)
        ex.lower_bound = 0.0
        ex.upper_bound = 1000.0
        ex.add_metabolites({shared_mets[ext_met_id]: -1.0})
        community.add_reactions([ex])

    # set objective = sum of all species objective reactions
    objective_dict = {}
    for species, rxn_map in species_objectives.items():
        for old_rxn_id, coef in rxn_map.items():
            new_id = f"{species}__{old_rxn_id}"
            if new_id in community.reactions:
                objective_dict[community.reactions.get_by_id(new_id)] = float(coef)

    community.objective = objective_dict

    return community, species_main_biomass


def run_diet(model: cobra.Model, diet_name: str, medium_dict: dict, species_main_biomass: dict):
    m = model.copy()

    applied = {rid: val for rid, val in medium_dict.items() if rid in m.reactions}
    m.medium = applied

    # enforce anaerobic gut assumption
    for rid in ["EX_o2(e)", "EX_o2_e", "EX_o2"]:
        if rid in m.reactions:
            m.reactions.get_by_id(rid).lower_bound = 0.0

    sol = m.optimize()

    summary = {
        "diet": diet_name,
        "status": sol.status,
        "community_objective": sol.objective_value if sol.status == "optimal" else None,
        "applied_exchanges": len(applied),
        "missing_exchanges": len(medium_dict) - len(applied),
    }

    species_rows = []
    for species, old_biomass_id in species_main_biomass.items():
        pref_biomass_id = f"{species}__{old_biomass_id}"
        flux = None
        objective_coefficient = 0.0
        objective_contribution = None
        objective_contribution_fraction = None
        if sol.status == "optimal" and pref_biomass_id in m.reactions:
            rxn = m.reactions.get_by_id(pref_biomass_id)
            flux = float(sol.fluxes[pref_biomass_id])
            objective_coefficient = float(rxn.objective_coefficient)
            objective_contribution = objective_coefficient * flux
            if summary["community_objective"] not in (None, 0):
                objective_contribution_fraction = objective_contribution / float(summary["community_objective"])
            else:
                objective_contribution_fraction = 0.0

        species_rows.append(
            {
                "diet": diet_name,
                "species": species,
                "biomass_reaction": pref_biomass_id,
                "objective_coefficient": objective_coefficient,
                "biomass_flux": flux,
                "objective_contribution": objective_contribution,
                "objective_contribution_fraction": objective_contribution_fraction,
            }
        )

    return summary, species_rows


def main():
    western = load_medium(WESTERN_CSV)
    fiber = load_medium(FIBER_CSV)
    medium_ids = set(western.keys()) | set(fiber.keys())

    model_files = sorted(MODELS_DIR.glob("*.xml"))
    print(f"Models found: {len(model_files)}")

    community, species_main_biomass = build_community_model(model_files, medium_ids)
    print(
        f"Community built: {len(community.reactions)} reactions, "
        f"{len(community.metabolites)} metabolites"
    )

    summary_rows = []
    species_rows_all = []

    for diet_name, medium in [("western", western), ("high_fiber", fiber)]:
        summary, species_rows = run_diet(community, diet_name, medium, species_main_biomass)
        summary_rows.append(summary)
        species_rows_all.extend(species_rows)
        print(
            f"{diet_name}: status={summary['status']}, "
            f"community_objective={summary['community_objective']}"
        )
        top = (
            pd.DataFrame(species_rows)
            .fillna({"objective_contribution": 0.0})
            .sort_values("objective_contribution", ascending=False)
            .head(5)
        )
        print("Top species contributions:")
        print(
            top[["species", "biomass_flux", "objective_coefficient", "objective_contribution"]]
            .to_string(index=False)
        )

    pd.DataFrame(summary_rows).to_csv(OUT_SUMMARY, index=False)
    pd.DataFrame(species_rows_all).to_csv(OUT_SPECIES, index=False)

    print(f"\nSaved: {OUT_SUMMARY}")
    print(f"Saved: {OUT_SPECIES}")


if __name__ == "__main__":
    main()
