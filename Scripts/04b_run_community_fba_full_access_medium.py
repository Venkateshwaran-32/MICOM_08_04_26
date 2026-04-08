import os
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
DEFAULT_UPTAKE_BOUND = 1000.0


# ------HELPERS------

# Clean a name so it is safe to reuse inside IDs.
def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", text) # basically replaces any character that is not a letter, number, or underscore with an underscore. This is useful for creating valid reaction IDs from arbitrary text.


# Read a diet CSV and turn it into {exchange_id: uptake_limit}.
def load_medium(csv_path: Path) -> dict:  
    df = pd.read_csv(csv_path)
    if not {"exchange_id", "max_uptake"}.issubset(df.columns):
        raise ValueError(f"{csv_path.name} must contain columns: exchange_id,max_uptake")
    return dict(zip(df["exchange_id"].astype(str), df["max_uptake"].astype(float)))
# basically reads the medium CSV file and returns a dictionary where the keys are exchange
# reaction IDs and the values are the maximum uptake rates.
# This is used to set the medium for the FBA simulations.


# Find which reaction(s) the model is trying to maximize, usually biomass.
def get_objective_reactions(model: cobra.Model) -> dict:
    obj = {
        rxn.id: float(rxn.objective_coefficient)
        for rxn in model.reactions
        if abs(float(rxn.objective_coefficient)) > 1e-12
    }
    if obj:
        return obj

    for rxn in model.reactions:
        if "biomass" in rxn.id.lower() or "biomass" in rxn.name.lower():
            return {rxn.id: 1.0}
    raise ValueError(f"Could not infer objective reaction for model {model.id}")


# Make a copy of one species model and prefix all IDs so species do not clash in the community model.
def prefixed_copy(model: cobra.Model, prefix: str) -> cobra.Model:
    m = model.copy()
    for met in list(m.metabolites):
        met.id = f"{prefix}__{met.id}"
    for rxn in list(m.reactions):
        rxn.id = f"{prefix}__{rxn.id}"
    for gene in list(m.genes):
        gene.id = f"{prefix}__{gene.id}"
    return m


# Create a shared-pool metabolite ID from an original extracellular metabolite ID.
def shared_met_id(original_ext_met_id: str) -> str:
    return f"u__{sanitize(original_ext_met_id)}"


def get_full_access_bound() -> float:
    raw = os.environ.get("FULL_ACCESS_UPTAKE_BOUND", str(DEFAULT_UPTAKE_BOUND))
    try:
        bound = float(raw)
    except ValueError as exc:
        raise ValueError(f"FULL_ACCESS_UPTAKE_BOUND must be numeric, got {raw!r}") from exc
    if bound <= 0:
        raise ValueError(f"FULL_ACCESS_UPTAKE_BOUND must be positive, got {bound}")
    return bound


def bound_label(bound: float) -> str:
    return str(int(bound)) if float(bound).is_integer() else str(bound).replace(".", "_")


def output_paths(bound: float) -> tuple[Path, Path]:
    if bound == DEFAULT_UPTAKE_BOUND:
        return (
            OUT_DIR / "community_growth_full_access.csv",
            OUT_DIR / "community_species_biomass_flux_full_access.csv",
        )
    suffix = f"_bound_{bound_label(bound)}"
    return (
        OUT_DIR / f"community_growth_full_access{suffix}.csv",
        OUT_DIR / f"community_species_biomass_flux_full_access{suffix}.csv",
    )


# Build one combined community model from all species models plus shared metabolites/connectors.
def build_community_model(model_files, medium_exchange_ids):
    community = Model("community_model_full_access")

    shared_mets = {}
    exchange_to_ext_met = {}
    species_objectives = {}
    species_main_biomass = {}

    for fp in model_files:
        base = cobra.io.read_sbml_model(str(fp))
        species = sanitize(fp.stem)

        obj_rxn_map = get_objective_reactions(base)
        species_objectives[species] = obj_rxn_map
        species_main_biomass[species] = next(iter(obj_rxn_map.keys()))

        for ex_id in medium_exchange_ids:
            if ex_id in base.reactions:
                rxn = base.reactions.get_by_id(ex_id)
                mets = list(rxn.metabolites.keys())
                if len(mets) == 1:
                    ext_met_id = mets[0].id
                    exchange_to_ext_met.setdefault(ex_id, ext_met_id)

        ext_met_ids = [m.id for m in base.metabolites if m.compartment == "e"]

        pref = prefixed_copy(base, species)
        community.add_reactions(pref.reactions)

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

    objective_dict = {}
    for species, rxn_map in species_objectives.items():
        for old_rxn_id, coef in rxn_map.items():
            new_id = f"{species}__{old_rxn_id}"
            if new_id in community.reactions:
                objective_dict[community.reactions.get_by_id(new_id)] = float(coef)

    community.objective = objective_dict

    return community, species_main_biomass


# Replace all diet uptake values with the same full-access bound while keeping the same exchange IDs.
def to_full_access_medium(medium_dict: dict, uptake_bound: float = 1000.0) -> dict:
    return {rid: float(uptake_bound) for rid in medium_dict}


# Run one diet scenario on the community model and collect community and species biomass outputs.
def run_diet(model: cobra.Model, diet_name: str, medium_dict: dict, species_main_biomass: dict):
    m = model.copy()

    applied = {rid: val for rid, val in medium_dict.items() if rid in m.reactions}
    m.medium = applied

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
        if sol.status == "optimal" and pref_biomass_id in m.reactions:
            flux = float(sol.fluxes[pref_biomass_id])

        species_rows.append(
            {
                "diet": diet_name,
                "species": species,
                "biomass_reaction": pref_biomass_id,
                "biomass_flux": flux,
            }
        )

    return summary, species_rows


def main():
    uptake_bound = get_full_access_bound()
    out_summary, out_species = output_paths(uptake_bound)

    western = to_full_access_medium(load_medium(WESTERN_CSV), uptake_bound)
    fiber = to_full_access_medium(load_medium(FIBER_CSV), uptake_bound)
    medium_ids = set(western.keys()) | set(fiber.keys())

    model_files = sorted(MODELS_DIR.glob("*.xml"))
    print(f"Models found: {len(model_files)}")
    print(f"Full-access uptake bound per exchange: {uptake_bound}")

    community, species_main_biomass = build_community_model(model_files, medium_ids)
    print(
        f"Community built: {len(community.reactions)} reactions, "
        f"{len(community.metabolites)} metabolites"
    )

    summary_rows = []
    species_rows_all = []

    for diet_name, medium in [("western_full_access", western), ("high_fiber_full_access", fiber)]:
        summary, species_rows = run_diet(community, diet_name, medium, species_main_biomass)
        summary["uptake_bound_per_exchange"] = uptake_bound
        for row in species_rows:
            row["uptake_bound_per_exchange"] = uptake_bound
        summary_rows.append(summary)
        species_rows_all.extend(species_rows)
        print(
            f"{diet_name}: status={summary['status']}, "
            f"community_objective={summary['community_objective']}"
        )

    pd.DataFrame(summary_rows).to_csv(out_summary, index=False)
    pd.DataFrame(species_rows_all).to_csv(out_species, index=False)

    print(f"\nSaved: {out_summary}")
    print(f"Saved: {out_species}")


if __name__ == "__main__":
    main()
