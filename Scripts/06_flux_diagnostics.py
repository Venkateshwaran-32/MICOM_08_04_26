from pathlib import Path
import importlib
import re
from typing import Dict, List, Optional, Tuple
import pandas as pd
import cobra
from cobra import Model, Reaction, Metabolite
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODELS_DIR = PROJECT_ROOT / "Models" / "vmh_agora_sbml"
MEDIA_DIR = PROJECT_ROOT / "Media"
OUT_DIR = PROJECT_ROOT / "Results" / "fba"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR = PROJECT_ROOT / "Results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

WESTERN_CSV = MEDIA_DIR / "western.csv"
FIBER_CSV = MEDIA_DIR / "high_fiber.csv"

OUT_COMM_EX = OUT_DIR / "community_exchange_fluxes_by_diet.csv"
OUT_COMM_CONN = OUT_DIR / "community_species_connector_fluxes_by_diet.csv"
OUT_INDIV = OUT_DIR / "individual_exchange_fluxes_by_diet.csv"
OUT_COMM_RXN = OUT_DIR / "community_reaction_fluxes_by_diet.csv"
OUT_PATHWAY_ACTIVITY = OUT_DIR / "community_pathway_activity_by_diet.csv"
OUT_BIOMASS_PATHWAYS = OUT_DIR / "community_biomass_associated_pathways_by_diet.csv"
OUT_PATHWAY_COMPARE = OUT_DIR / "community_pathway_diet_comparison.csv"

FLUX_EPS = 1e-9


def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", text)


def load_medium(csv_path: Path) -> Dict[str, float]:
    df = pd.read_csv(csv_path)
    if not {"exchange_id", "max_uptake"}.issubset(df.columns):
        raise ValueError(f"{csv_path.name} must contain columns: exchange_id,max_uptake")
    return dict(zip(df["exchange_id"].astype(str), df["max_uptake"].astype(float)))


def apply_medium_bounds(model: cobra.Model, medium: Dict[str, float]) -> Dict[str, float]:
    """
    Apply uptake bounds without using model.medium to avoid COBRA warnings
    for custom community exchange reactions.
    """
    # Match medium behavior: close uptake on all recognized exchange reactions first.
    for rxn in model.exchanges:
        rxn.lower_bound = 0.0

    applied = {}
    for rid, val in medium.items():
        if rid not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(rid)
        rxn.lower_bound = -float(val)
        applied[rid] = float(val)
    return applied


def parse_carbon_count(formula: Optional[str]) -> Optional[int]:
    if not formula or not isinstance(formula, str):
        return None
    total_c = 0
    for elem, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if elem != "C":
            continue
        total_c += int(n) if n else 1
    return total_c


def get_exchange_met(rxn: cobra.Reaction):
    mets = list(rxn.metabolites.keys())
    return mets[0] if len(mets) == 1 else None


def first_nonempty_text(value) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, (list, tuple, set)):
        items = [str(x).strip() for x in value if str(x).strip()]
        return " | ".join(items) if items else None
    text = str(value).strip()
    return text if text else None


def get_reaction_pathway(rxn: cobra.Reaction) -> str:
    candidate = first_nonempty_text(getattr(rxn, "subsystem", None))
    if candidate:
        return candidate

    ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
    for key in ["subsystem", "pathway", "kegg.pathway", "metacyc.reaction", "sbo"]:
        candidate = first_nonempty_text(ann.get(key))
        if candidate:
            return candidate
    return "Unassigned"


def infer_species_and_local_id(rxn_id: str, species_ids: List[str]) -> Tuple[Optional[str], Optional[str]]:
    # Longest-prefix match prevents ambiguity when one species ID is a prefix of another.
    for species in sorted(species_ids, key=len, reverse=True):
        prefix = f"{species}__"
        if rxn_id.startswith(prefix):
            return species, rxn_id[len(prefix) :]
    return None, None


def get_objective_reactions(model: cobra.Model) -> Dict[str, float]:
    obj = {
        rxn.id: float(rxn.objective_coefficient)
        for rxn in model.reactions
        if abs(float(rxn.objective_coefficient)) > 1e-12
    }
    if obj:
        return obj
    for rxn in model.reactions:
        nm = (rxn.name or "").lower()
        if "biomass" in rxn.id.lower() or "biomass" in nm:
            return {rxn.id: 1.0}
    raise ValueError(f"Could not infer objective reaction for model {model.id}")


def prefixed_copy(model: cobra.Model, prefix: str) -> cobra.Model:
    m = model.copy()
    # Fast bulk rename: set private IDs, then rebuild indices once.
    # Using property setters here is very slow on large models because each set
    # triggers index maintenance.
    for met in m.metabolites:
        met._id = f"{prefix}__{met.id}"
    m.metabolites._generate_index()

    for rxn in m.reactions:
        rxn._id = f"{prefix}__{rxn.id}"
    m.reactions._generate_index()

    for gene in m.genes:
        gene._id = f"{prefix}__{gene.id}"
    m.genes._generate_index()

    m.repair()
    return m


def shared_met_id(original_ext_met_id: str) -> str:
    return f"u__{sanitize(original_ext_met_id)}"


def build_community_model(model_files, medium_exchange_ids):
    community = Model("community_model")
    shared_mets = {}
    exchange_to_ext_met = {}
    species_objectives = {}
    connector_meta = {}
    species_ids = []
    biomass_rxn_ids = set()

    for fp in model_files:
        base = cobra.io.read_sbml_model(str(fp))
        species = sanitize(fp.stem)
        species_ids.append(species)

        obj_rxn_map = get_objective_reactions(base)
        species_objectives[species] = obj_rxn_map

        for ex_id in medium_exchange_ids:
            if ex_id in base.reactions:
                rxn = base.reactions.get_by_id(ex_id)
                mets = list(rxn.metabolites.keys())
                if len(mets) == 1:
                    exchange_to_ext_met.setdefault(ex_id, mets[0].id)

        ext_met_ids = [m.id for m in base.metabolites if m.compartment == "e"]
        pref = prefixed_copy(base, species)
        community.add_reactions(pref.reactions)

        for ext_id in ext_met_ids:
            local_id = f"{species}__{ext_id}"
            if local_id not in community.metabolites:
                continue

            if ext_id not in shared_mets:
                base_ext_met = base.metabolites.get_by_id(ext_id)
                sm = Metabolite(
                    id=shared_met_id(ext_id),
                    name=f"shared_{base_ext_met.name or ext_id}",
                    formula=base_ext_met.formula,
                    charge=base_ext_met.charge,
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
            base_ext_met = base.metabolites.get_by_id(ext_id)
            conn.add_metabolites(
                {
                    community.metabolites.get_by_id(local_id): -1.0,
                    shared_mets[ext_id]: 1.0,
                }
            )
            community.add_reactions([conn])
            connector_meta[conn_id] = {
                "species": species,
                "external_metabolite_id": ext_id,
                "external_metabolite_name": base_ext_met.name,
                "external_metabolite_formula": base_ext_met.formula,
                "carbon_count": parse_carbon_count(base_ext_met.formula),
            }

    for ex_id, ext_met_id in exchange_to_ext_met.items():
        if ex_id in community.reactions or ext_met_id not in shared_mets:
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
                biomass_rxn_ids.add(new_id)
    community.objective = objective_dict

    return community, connector_meta, species_ids, biomass_rxn_ids


def summarize_pathway_activity(df_rxn: pd.DataFrame) -> pd.DataFrame:
    if df_rxn.empty:
        return pd.DataFrame(
            columns=[
                "diet",
                "species",
                "pathway",
                "sum_abs_flux",
                "net_flux",
                "n_active_reactions",
                "n_biomass_reactions",
            ]
        )

    use = df_rxn[df_rxn["reaction_type"] == "species_internal"].copy()
    if use.empty:
        return pd.DataFrame(
            columns=[
                "diet",
                "species",
                "pathway",
                "sum_abs_flux",
                "net_flux",
                "n_active_reactions",
                "n_biomass_reactions",
            ]
        )

    out = (
        use.groupby(["diet", "species", "pathway"], as_index=False)
        .agg(
            sum_abs_flux=("abs_flux", "sum"),
            net_flux=("flux", "sum"),
            n_active_reactions=("reaction_id", "count"),
            n_biomass_reactions=("is_biomass_reaction", "sum"),
        )
        .sort_values(["diet", "species", "sum_abs_flux"], ascending=[True, True, False])
    )
    return out


def summarize_biomass_associated_pathways(df_rxn: pd.DataFrame) -> pd.DataFrame:
    if df_rxn.empty:
        return pd.DataFrame(
            columns=[
                "diet",
                "species",
                "pathway",
                "pathway_abs_flux",
                "species_internal_abs_flux_total",
                "pathway_flux_share",
                "species_biomass_flux_total",
                "biomass_association_score",
            ]
        )

    internal = df_rxn[df_rxn["reaction_type"] == "species_internal"].copy()
    if internal.empty:
        return pd.DataFrame(
            columns=[
                "diet",
                "species",
                "pathway",
                "pathway_abs_flux",
                "species_internal_abs_flux_total",
                "pathway_flux_share",
                "species_biomass_flux_total",
                "biomass_association_score",
            ]
        )

    total_internal = (
        internal.groupby(["diet", "species"], as_index=False)["abs_flux"]
        .sum()
        .rename(columns={"abs_flux": "species_internal_abs_flux_total"})
    )
    biomass_flux = (
        df_rxn[df_rxn["is_biomass_reaction"]]
        .groupby(["diet", "species"], as_index=False)["flux"]
        .sum()
        .rename(columns={"flux": "species_biomass_flux_total"})
    )
    by_path = (
        internal.groupby(["diet", "species", "pathway"], as_index=False)["abs_flux"]
        .sum()
        .rename(columns={"abs_flux": "pathway_abs_flux"})
    )

    out = by_path.merge(total_internal, on=["diet", "species"], how="left").merge(
        biomass_flux, on=["diet", "species"], how="left"
    )
    out["species_biomass_flux_total"] = out["species_biomass_flux_total"].fillna(0.0)
    out["pathway_flux_share"] = np.where(
        out["species_internal_abs_flux_total"] > 0,
        out["pathway_abs_flux"] / out["species_internal_abs_flux_total"],
        0.0,
    )
    out["biomass_association_score"] = out["pathway_flux_share"] * out["species_biomass_flux_total"]
    out = out.sort_values(["diet", "species", "biomass_association_score"], ascending=[True, True, False])
    return out


def compare_pathways_between_diets(df_pathway: pd.DataFrame) -> pd.DataFrame:
    if df_pathway.empty:
        return pd.DataFrame(
            columns=[
                "species",
                "pathway",
                "western_sum_abs_flux",
                "high_fiber_sum_abs_flux",
                "delta_high_fiber_minus_western",
                "log2_fc_high_fiber_vs_western",
            ]
        )

    pivot = (
        df_pathway.pivot_table(
            index=["species", "pathway"],
            columns="diet",
            values="sum_abs_flux",
            aggfunc="sum",
            fill_value=0.0,
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )

    if "western" not in pivot.columns:
        pivot["western"] = 0.0
    if "high_fiber" not in pivot.columns:
        pivot["high_fiber"] = 0.0

    eps = 1e-9
    pivot["delta_high_fiber_minus_western"] = pivot["high_fiber"] - pivot["western"]
    pivot["log2_fc_high_fiber_vs_western"] = np.log2((pivot["high_fiber"] + eps) / (pivot["western"] + eps))

    out = pivot.rename(
        columns={
            "western": "western_sum_abs_flux",
            "high_fiber": "high_fiber_sum_abs_flux",
        }
    ).sort_values("delta_high_fiber_minus_western", ascending=False)
    return out


def make_pathway_figures(df_pathway: pd.DataFrame, df_compare: pd.DataFrame):
    try:
        plt = importlib.import_module("matplotlib.pyplot")
    except Exception as e:
        print(f"Skipping pathway figures: matplotlib unavailable ({type(e).__name__})")
        return

    if not df_pathway.empty:
        top = (
            df_pathway.groupby("pathway", as_index=False)["sum_abs_flux"]
            .sum()
            .sort_values("sum_abs_flux", ascending=False)
            .head(20)["pathway"]
            .tolist()
        )
        heat = (
            df_pathway[df_pathway["pathway"].isin(top)]
            .groupby(["pathway", "diet"], as_index=False)["sum_abs_flux"]
            .sum()
            .pivot(index="pathway", columns="diet", values="sum_abs_flux")
            .fillna(0.0)
        )
        if not heat.empty:
            fig, ax = plt.subplots(figsize=(8, 8))
            im = ax.imshow(heat.values, aspect="auto")
            ax.set_yticks(range(len(heat.index)))
            ax.set_yticklabels(heat.index)
            ax.set_xticks(range(len(heat.columns)))
            ax.set_xticklabels(heat.columns, rotation=45, ha="right")
            ax.set_title("Top Pathway Activity by Diet")
            fig.colorbar(im, ax=ax, label="sum_abs_flux")
            fig.tight_layout()
            fig.savefig(FIG_DIR / "pathway_activity_heatmap_top20.png", dpi=180)
            plt.close(fig)

    if not df_compare.empty:
        top_delta = (
            df_compare.assign(abs_delta=lambda d: d["delta_high_fiber_minus_western"].abs())
            .sort_values("abs_delta", ascending=False)
            .head(20)
            .copy()
        )
        if not top_delta.empty:
            top_delta["label"] = top_delta["species"] + " | " + top_delta["pathway"]
            fig, ax = plt.subplots(figsize=(10, 7))
            ax.barh(top_delta["label"], top_delta["delta_high_fiber_minus_western"])
            ax.set_title("Top Pathway Flux Changes (high_fiber - western)")
            ax.set_xlabel("delta sum_abs_flux")
            ax.invert_yaxis()
            fig.tight_layout()
            fig.savefig(FIG_DIR / "pathway_diet_delta_top20.png", dpi=180)
            plt.close(fig)


def run_community_flux_diagnostics(
    community: cobra.Model,
    connector_meta: dict[str, dict],
    species_ids: list[str],
    biomass_rxn_ids: set[str],
    diet_name: str,
    medium: dict[str, float],
):
    m = community.copy()
    applied = apply_medium_bounds(m, medium)

    for rid in ["EX_o2(e)", "EX_o2_e", "EX_o2"]:
        if rid in m.reactions:
            m.reactions.get_by_id(rid).lower_bound = 0.0

    sol = m.optimize()
    if sol.status != "optimal":
        raise RuntimeError(f"Community optimization failed for {diet_name}: {sol.status}")

    ex_rows = []
    for ex_id in sorted(applied.keys()):
        rxn = m.reactions.get_by_id(ex_id)
        flux = float(sol.fluxes[ex_id])
        met = get_exchange_met(rxn)
        formula = met.formula if met is not None else None
        c_count = parse_carbon_count(formula)
        ex_rows.append(
            {
                "diet": diet_name,
                "exchange_id": ex_id,
                "uptake_bound": medium[ex_id],
                "flux": flux,
                "uptake_flux": max(0.0, -flux),
                "secretion_flux": max(0.0, flux),
                "abs_flux": abs(flux),
                "metabolite_id": (met.id if met is not None else None),
                "metabolite_formula": formula,
                "carbon_count": c_count,
                "carbon_uptake_flux": (max(0.0, -flux) * c_count) if c_count is not None else None,
                "at_uptake_bound": abs(max(0.0, -flux) - medium[ex_id]) <= 1e-6,
            }
        )

    conn_rows = []
    for rxn in m.reactions:
        if not rxn.id.startswith("COMM__"):
            continue
        flux = float(sol.fluxes[rxn.id])
        if abs(flux) <= FLUX_EPS:
            continue
        md = connector_meta.get(rxn.id, {})
        species = md.get("species")
        c_count = md.get("carbon_count")
        uptake = max(0.0, -flux)
        secretion = max(0.0, flux)
        conn_rows.append(
            {
                "diet": diet_name,
                "connector_reaction": rxn.id,
                "species": species,
                "external_metabolite_id": md.get("external_metabolite_id"),
                "external_metabolite_name": md.get("external_metabolite_name"),
                "external_metabolite_formula": md.get("external_metabolite_formula"),
                "flux": flux,
                "uptake_from_shared_flux": uptake,
                "secretion_to_shared_flux": secretion,
                "carbon_count": c_count,
                "carbon_uptake_from_shared_flux": (uptake * c_count) if c_count is not None else None,
                "carbon_secretion_to_shared_flux": (secretion * c_count) if c_count is not None else None,
                "abs_flux": abs(flux),
            }
        )

    rxn_rows = []
    for rxn in m.reactions:
        flux = float(sol.fluxes[rxn.id])
        if abs(flux) <= FLUX_EPS:
            continue

        species, local_id = infer_species_and_local_id(rxn.id, species_ids)
        if rxn.id.startswith("COMM__"):
            reaction_type = "connector_shared"
        elif species is None and rxn.id in applied:
            reaction_type = "community_exchange"
        elif species is not None and local_id is not None and local_id.startswith("EX_"):
            reaction_type = "species_exchange"
        elif species is not None:
            reaction_type = "species_internal"
        else:
            reaction_type = "other"

        is_biomass = (rxn.id in biomass_rxn_ids) or ("biomass" in rxn.id.lower()) or ("biomass" in (rxn.name or "").lower())
        rxn_rows.append(
            {
                "diet": diet_name,
                "reaction_id": rxn.id,
                "reaction_name": rxn.name,
                "species": species,
                "local_reaction_id": local_id,
                "reaction_type": reaction_type,
                "pathway": get_reaction_pathway(rxn),
                "is_biomass_reaction": bool(is_biomass),
                "flux": flux,
                "abs_flux": abs(flux),
            }
        )

    return ex_rows, conn_rows, rxn_rows


def run_individual_exchange_fluxes(model_files, diet_name: str, medium: dict[str, float]):
    rows = []
    medium_ids = set(medium.keys())

    for fp in sorted(model_files):
        model = cobra.io.read_sbml_model(str(fp))
        species = fp.stem

        for rid in ["EX_o2(e)", "EX_o2_e", "EX_o2"]:
            if rid in model.reactions:
                model.reactions.get_by_id(rid).lower_bound = 0.0

        applied = apply_medium_bounds(model, medium)

        sol = model.optimize()
        status = sol.status

        if status != "optimal":
            rows.append(
                {
                    "diet": diet_name,
                    "species": species,
                    "status": status,
                    "exchange_id": None,
                    "flux": None,
                    "uptake_flux": None,
                    "secretion_flux": None,
                    "abs_flux": None,
                    "metabolite_id": None,
                    "metabolite_formula": None,
                    "carbon_count": None,
                    "carbon_uptake_flux": None,
                    "in_diet_medium": None,
                    "applied_exchanges": len(applied),
                    "missing_exchanges": len(medium) - len(applied),
                }
            )
            continue

        for rxn in model.exchanges:
            flux = float(sol.fluxes[rxn.id])
            if abs(flux) <= FLUX_EPS:
                continue

            met = get_exchange_met(rxn)
            formula = met.formula if met is not None else None
            c_count = parse_carbon_count(formula)
            uptake = max(0.0, -flux)
            bound = medium.get(rxn.id) if rxn.id in medium else None

            rows.append(
                {
                    "diet": diet_name,
                    "species": species,
                    "status": status,
                    "exchange_id": rxn.id,
                    "uptake_bound": bound,
                    "flux": flux,
                    "uptake_flux": uptake,
                    "secretion_flux": max(0.0, flux),
                    "abs_flux": abs(flux),
                    "metabolite_id": (met.id if met is not None else None),
                    "metabolite_formula": formula,
                    "carbon_count": c_count,
                    "carbon_uptake_flux": (uptake * c_count) if c_count is not None else None,
                    "in_diet_medium": rxn.id in medium_ids,
                    "at_uptake_bound": (bound is not None and abs(uptake - bound) <= 1e-6),
                    "applied_exchanges": len(applied),
                    "missing_exchanges": len(medium) - len(applied),
                }
            )

    return rows


def main():
    western = load_medium(WESTERN_CSV)
    fiber = load_medium(FIBER_CSV)

    model_files = sorted(MODELS_DIR.glob("*.xml"))
    medium_ids = set(western.keys()) | set(fiber.keys())

    community, connector_meta, species_ids, biomass_rxn_ids = build_community_model(model_files, medium_ids)

    comm_ex_rows_all = []
    comm_conn_rows_all = []
    comm_rxn_rows_all = []
    indiv_rows_all = []

    for diet_name, medium in [("western", western), ("high_fiber", fiber)]:
        ex_rows, conn_rows, rxn_rows = run_community_flux_diagnostics(
            community, connector_meta, species_ids, biomass_rxn_ids, diet_name, medium
        )
        comm_ex_rows_all.extend(ex_rows)
        comm_conn_rows_all.extend(conn_rows)
        comm_rxn_rows_all.extend(rxn_rows)

        indiv_rows = run_individual_exchange_fluxes(model_files, diet_name, medium)
        indiv_rows_all.extend(indiv_rows)

    df_ex = pd.DataFrame(comm_ex_rows_all).sort_values(["diet", "uptake_flux"], ascending=[True, False])
    df_conn = pd.DataFrame(comm_conn_rows_all).sort_values(["diet", "abs_flux"], ascending=[True, False])
    df_rxn = pd.DataFrame(comm_rxn_rows_all).sort_values(["diet", "species", "abs_flux"], ascending=[True, True, False])
    df_ind = pd.DataFrame(indiv_rows_all).sort_values(["diet", "species", "abs_flux"], ascending=[True, True, False])
    df_pathway = summarize_pathway_activity(df_rxn)
    df_biomass_pathway = summarize_biomass_associated_pathways(df_rxn)
    df_compare = compare_pathways_between_diets(df_pathway)

    df_ex.to_csv(OUT_COMM_EX, index=False)
    df_conn.to_csv(OUT_COMM_CONN, index=False)
    df_rxn.to_csv(OUT_COMM_RXN, index=False)
    df_ind.to_csv(OUT_INDIV, index=False)
    df_pathway.to_csv(OUT_PATHWAY_ACTIVITY, index=False)
    df_biomass_pathway.to_csv(OUT_BIOMASS_PATHWAYS, index=False)
    df_compare.to_csv(OUT_PATHWAY_COMPARE, index=False)
    make_pathway_figures(df_pathway, df_compare)

    print(f"Saved: {OUT_COMM_EX} ({len(df_ex)} rows)")
    print(f"Saved: {OUT_COMM_CONN} ({len(df_conn)} rows)")
    print(f"Saved: {OUT_COMM_RXN} ({len(df_rxn)} rows)")
    print(f"Saved: {OUT_INDIV} ({len(df_ind)} rows)")
    print(f"Saved: {OUT_PATHWAY_ACTIVITY} ({len(df_pathway)} rows)")
    print(f"Saved: {OUT_BIOMASS_PATHWAYS} ({len(df_biomass_pathway)} rows)")
    print(f"Saved: {OUT_PATHWAY_COMPARE} ({len(df_compare)} rows)")
    print(f"Figures dir: {FIG_DIR}")


if __name__ == "__main__":
    main()
