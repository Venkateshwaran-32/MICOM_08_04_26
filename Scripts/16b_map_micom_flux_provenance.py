from pathlib import Path
import argparse
import re
from collections import Counter, defaultdict
from functools import lru_cache

import cobra
import pandas as pd


# -------------------------------------------------------------------
# Script 16b: build MICOM provenance summaries for selected species
# Inputs:
# - Results/micom/pathway_flux/proper_age_bins/reaction_fluxes_long_by_agegroup_diet.csv
# - Results/micom/growth/proper_age_bins/organism_growth_rates_by_agegroup_diet.csv
# - Results/micom/tradeoff_sensitivity/proper_age_bins/tradeoff_sensitivity_summary.csv
# - Media/high_fiber.csv
# - Media/western.csv
# Outputs:
# - Results/micom/provenance/proper_age_bins/*__pathway_overview.csv
# - Results/micom/provenance/proper_age_bins/*__reaction_steps.csv
# - Results/micom/provenance/proper_age_bins/*__exchange_summary.csv
# - Results/micom/provenance/proper_age_bins/*__crossfeeding_candidates.csv
# - Results/micom/provenance/proper_age_bins/*__provenance_edges.csv
# - Results/micom/provenance/proper_age_bins/*__summary.txt
# - Results/micom/provenance/proper_age_bins/micom_provenance_case_overview.csv
# Run:
# - HOME="$PWD" python Scripts/16b_map_micom_flux_provenance.py
# -------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODELS_DIR = PROJECT_ROOT / "Models" / "vmh_agora_sbml"
MICOM_DIR = PROJECT_ROOT / "Results" / "micom"
IN_FLUX = MICOM_DIR / "pathway_flux" / "proper_age_bins" / "reaction_fluxes_long_by_agegroup_diet.csv"
IN_GROWTH = MICOM_DIR / "growth" / "proper_age_bins" / "organism_growth_rates_by_agegroup_diet.csv"
IN_TRADEOFF = MICOM_DIR / "tradeoff_sensitivity" / "proper_age_bins" / "tradeoff_sensitivity_summary.csv"
MEDIA_DIR = PROJECT_ROOT / "Media"
WESTERN_CSV = MEDIA_DIR / "western.csv"
FIBER_CSV = MEDIA_DIR / "high_fiber.csv"
OUT_DIR = MICOM_DIR / "provenance" / "proper_age_bins"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CASE_OVERVIEW = OUT_DIR / "micom_provenance_case_overview.csv"

FLUX_EPS = 1e-9
DEFAULT_AGE_GROUPS = ["61_70", "71_80", "81_plus"]
DEFAULT_DIETS = ["high_fiber", "western"]
DEFAULT_SPECIES = [
    "Faecalibacterium_prausnitzii_M21_2__AGORA1_03",
    "Escherichia_coli_UTI89_UPEC_AGORA1_03",
]
DEFAULT_TRADEOFF = 0.5
TOP_N_PATHWAYS = 6
TOP_N_EXCHANGES = 15
COMMON_COFATOR_METS = {
    "h[c]", "h[e]", "h2o[c]", "h2o[e]", "pi[c]", "pi[e]", "ppi[c]", "ppi[e]",
    "atp[c]", "atp[e]", "adp[c]", "adp[e]", "amp[c]", "amp[e]",
    "nad[c]", "nad[e]", "nadh[c]", "nadh[e]", "nadp[c]", "nadp[e]", "nadph[c]", "nadph[e]",
    "coa[c]", "coa[e]", "co2[c]", "co2[e]", "o2[c]", "o2[e]",
}
EXCLUDED_EXTERNAL_METS = {
    "h(e)", "h2o(e)", "co2(e)", "pi(e)", "ppi(e)", "nh4(e)",
    "h[e]", "h2o[e]", "co2[e]", "pi[e]", "ppi[e]", "nh4[e]",
}
GENERIC_KEYWORDS = {
    "metabolism", "pathway", "biosynthesis", "degradation", "transport", "exchange",
    "reversible", "reaction", "putative", "protein", "enzyme", "and", "via",
    "for", "the", "with", "without", "cofactor", "associated",
}


def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", str(text))


def slug(text: str) -> str:
    return sanitize(text).strip("_")


def first_nonempty(value):
    if value is None:
        return None
    if isinstance(value, (list, tuple, set)):
        vals = [str(v).strip() for v in value if str(v).strip()]
        return " | ".join(vals) if vals else None
    s = str(value).strip()
    return s if s else None


def reaction_pathway(rxn: cobra.Reaction) -> str:
    cand = first_nonempty(getattr(rxn, "subsystem", None))
    if cand:
        return cand
    ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
    for key in ["subsystem", "pathway", "kegg.pathway", "metacyc.reaction", "sbo"]:
        cand = first_nonempty(ann.get(key))
        if cand:
            return cand
    return "Unassigned"


@lru_cache(maxsize=None)
def load_species_model(species_id: str) -> cobra.Model:
    model_files = sorted(MODELS_DIR.glob("*.xml"))
    for fp in model_files:
        if sanitize(fp.stem) == species_id:
            return cobra.io.read_sbml_model(str(fp))
    raise FileNotFoundError(f"Could not find SBML model for species id: {species_id}")


def directed_sides(rxn: cobra.Reaction, flux: float):
    if flux >= 0:
        left = [(met.id, -coef) for met, coef in rxn.metabolites.items() if coef < 0]
        right = [(met.id, coef) for met, coef in rxn.metabolites.items() if coef > 0]
    else:
        left = [(met.id, coef) for met, coef in rxn.metabolites.items() if coef > 0]
        right = [(met.id, -coef) for met, coef in rxn.metabolites.items() if coef < 0]
    return left, right


def side_to_text(items):
    if not items:
        return "0"
    pieces = []
    for met_id, coeff in sorted(items):
        if abs(coeff - 1.0) <= 1e-12:
            pieces.append(met_id)
        else:
            pieces.append(f"{coeff:g} {met_id}")
    return " + ".join(pieces)


def build_edges(records):
    edges = []
    for src in records:
        src_products = {met_id for met_id, _ in src["products"]}
        if not src_products:
            continue
        for dst in records:
            if src["reaction_id"] == dst["reaction_id"]:
                continue
            dst_substrates = {met_id for met_id, _ in dst["substrates"]}
            shared = sorted((src_products & dst_substrates) - COMMON_COFATOR_METS)
            if not shared:
                continue
            edges.append(
                {
                    "from_reaction_id": src["reaction_id"],
                    "to_reaction_id": dst["reaction_id"],
                    "shared_metabolites": " | ".join(shared),
                    "n_shared_metabolites": len(shared),
                }
            )
    return pd.DataFrame(edges)


def build_stepwise_pipeline(records):
    by_id = {rec["reaction_id"]: rec for rec in records}
    edges = build_edges(records)
    if edges.empty:
        return pd.DataFrame(
            [
                {
                    "step": i,
                    "reaction_id": rec["reaction_id"],
                    "reaction_name": rec["reaction_name"],
                    "flux": rec["flux"],
                    "equation_in_flux_direction": rec["equation_in_flux_direction"],
                    "incoming_from": "",
                    "shared_metabolites_from_previous": "",
                }
                for i, rec in enumerate(sorted(records, key=lambda r: (-r["abs_flux"], r["reaction_id"])), start=1)
            ]
        )

    edge_list = edges.to_dict("records")
    in_degree = Counter()
    adjacency = defaultdict(list)
    reverse_adjacency = defaultdict(list)
    for edge in edge_list:
        src = edge["from_reaction_id"]
        dst = edge["to_reaction_id"]
        adjacency[src].append(edge)
        reverse_adjacency[dst].append(edge)
        in_degree[dst] += 1

    candidate_starts = [
        rec["reaction_id"]
        for rec in sorted(records, key=lambda r: (-r["abs_flux"], r["reaction_id"]))
        if in_degree[rec["reaction_id"]] == 0
    ]
    current = candidate_starts[0] if candidate_starts else sorted(
        records, key=lambda r: (-r["abs_flux"], r["reaction_id"])
    )[0]["reaction_id"]

    visited = set()
    ordered_rows = []
    step = 1
    while current and current not in visited:
        visited.add(current)
        rec = by_id[current]
        incoming = sorted(
            reverse_adjacency.get(current, []),
            key=lambda e: (e["n_shared_metabolites"], abs(by_id[e["from_reaction_id"]]["flux"])),
            reverse=True,
        )
        prev_rxn = incoming[0]["from_reaction_id"] if incoming else ""
        prev_mets = incoming[0]["shared_metabolites"] if incoming else ""
        ordered_rows.append(
            {
                "step": step,
                "reaction_id": rec["reaction_id"],
                "reaction_name": rec["reaction_name"],
                "flux": rec["flux"],
                "equation_in_flux_direction": rec["equation_in_flux_direction"],
                "incoming_from": prev_rxn,
                "shared_metabolites_from_previous": prev_mets,
            }
        )
        step += 1

        outgoing = [e for e in adjacency.get(current, []) if e["to_reaction_id"] not in visited]
        if not outgoing:
            current = None
            continue
        outgoing.sort(
            key=lambda e: (e["n_shared_metabolites"], abs(by_id[e["to_reaction_id"]]["flux"])),
            reverse=True,
        )
        current = outgoing[0]["to_reaction_id"]

    remaining = [
        rec for rec in sorted(records, key=lambda r: (-r["abs_flux"], r["reaction_id"]))
        if rec["reaction_id"] not in visited
    ]
    for rec in remaining:
        ordered_rows.append(
            {
                "step": step,
                "reaction_id": rec["reaction_id"],
                "reaction_name": rec["reaction_name"],
                "flux": rec["flux"],
                "equation_in_flux_direction": rec["equation_in_flux_direction"],
                "incoming_from": "",
                "shared_metabolites_from_previous": "",
            }
        )
        step += 1
    return pd.DataFrame(ordered_rows)


def load_medium(csv_path: Path) -> set[str]:
    df = pd.read_csv(csv_path)
    if not {"exchange_id", "max_uptake"}.issubset(df.columns):
        raise ValueError(f"{csv_path.name} must contain columns exchange_id,max_uptake")
    return set(df.loc[df["max_uptake"].fillna(0.0) > 0, "exchange_id"].astype(str))


def normalize_external_met_id(text: str) -> str:
    if pd.isna(text):
        return ""
    s = str(text).strip()
    if not s:
        return ""
    s = s.replace("[", "(").replace("]", ")")
    return s


def tokenize(*parts) -> set[str]:
    text = " ".join([str(p).lower() for p in parts if pd.notna(p)])
    tokens = set(re.findall(r"[a-z0-9_]+", text))
    out = set()
    for tok in tokens:
        if len(tok) < 3 or tok in GENERIC_KEYWORDS:
            continue
        out.add(tok)
    return out


def external_metabolites(rxn: cobra.Reaction):
    rows = []
    for met, coef in rxn.metabolites.items():
        if getattr(met, "compartment", "") != "e":
            continue
        rows.append(
            {
                "metabolite_id": normalize_external_met_id(met.id),
                "metabolite_name": met.name or normalize_external_met_id(met.id),
                "coefficient": float(coef),
            }
        )
    return rows


def reaction_direction_label(rxn: cobra.Reaction, flux: float) -> str:
    if not rxn.id.startswith("EX_"):
        return "internal"
    ext = external_metabolites(rxn)
    if not ext:
        return "internal"
    coef = ext[0]["coefficient"]
    if flux * coef < 0:
        return "uptake"
    if flux * coef > 0:
        return "secretion"
    return "inactive"


def build_external_flux_table(df_flux: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in df_flux.itertuples(index=False):
        model = load_species_model(row.id)
        if row.reaction_id not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(row.reaction_id)
        ext = external_metabolites(rxn)
        if not ext:
            continue
        direction = reaction_direction_label(rxn, float(row.flux))
        for met in ext:
            rows.append(
                {
                    "age_group": row.age_group,
                    "diet": row.diet,
                    "species": row.id,
                    "reaction_id": row.reaction_id,
                    "reaction_name": row.reaction_name,
                    "pathway": row.pathway,
                    "flux": float(row.flux),
                    "abs_flux": float(row.abs_flux),
                    "direction": direction,
                    "external_metabolite_id": met["metabolite_id"],
                    "external_metabolite_name": met["metabolite_name"],
                }
            )
    return pd.DataFrame(rows)


def validate_tradeoff_context(df_tradeoff: pd.DataFrame, age_groups: list[str], diets: list[str], fraction: float):
    sub = df_tradeoff[
        (df_tradeoff["age_group"].isin(age_groups))
        & (df_tradeoff["diet"].isin(diets))
        & (df_tradeoff["tradeoff_fraction"].round(6) == round(fraction, 6))
    ].copy()
    if sub.empty:
        raise ValueError(f"Could not find tradeoff_fraction={fraction} in MICOM tradeoff sensitivity outputs.")
    bad = sub[sub["status"].astype(str) != "optimal"]
    if not bad.empty:
        raise ValueError("Tradeoff sensitivity summary contains non-optimal rows for the selected fraction.")


def build_pathway_overview(df_case_flux: pd.DataFrame, growth_rate: float) -> pd.DataFrame:
    internal = df_case_flux[
        (df_case_flux["abs_flux"] > FLUX_EPS)
        & (~df_case_flux["pathway"].isin(["Transport, extracellular", "Exchange/demand reaction"]))
    ].copy()
    if internal.empty:
        return pd.DataFrame(
            columns=[
                "age_group", "diet", "species", "pathway", "pathway_abs_flux",
                "species_internal_abs_flux_total", "species_growth_rate", "pathway_flux_share",
                "biomass_association_score", "top_reaction_id", "top_reaction_name", "top_reaction_abs_flux",
            ]
        )

    keys = ["age_group", "diet", "id"]
    total_internal = (
        internal.groupby(keys, as_index=False)["abs_flux"]
        .sum()
        .rename(columns={"id": "species", "abs_flux": "species_internal_abs_flux_total"})
    )
    by_path = (
        internal.groupby(keys + ["pathway"], as_index=False)["abs_flux"]
        .sum()
        .rename(columns={"id": "species", "abs_flux": "pathway_abs_flux"})
    )
    top_rxn = (
        internal.sort_values(["pathway", "abs_flux"], ascending=[True, False])
        .groupby("pathway", as_index=False)
        .first()[["pathway", "reaction_id", "reaction_name", "abs_flux"]]
        .rename(
            columns={
                "reaction_id": "top_reaction_id",
                "reaction_name": "top_reaction_name",
                "abs_flux": "top_reaction_abs_flux",
            }
        )
    )
    out = by_path.merge(total_internal, on=["age_group", "diet", "species"], how="left")
    out["species_growth_rate"] = float(growth_rate)
    out["pathway_flux_share"] = out["pathway_abs_flux"] / out["species_internal_abs_flux_total"]
    out["biomass_association_score"] = out["pathway_flux_share"] * out["species_growth_rate"]
    out = out.merge(top_rxn, on="pathway", how="left")
    out = out.sort_values("biomass_association_score", ascending=False).head(TOP_N_PATHWAYS).reset_index(drop=True)
    return out


def trace_internal_pathway(df_case_flux: pd.DataFrame, species: str, age_group: str, diet: str, pathway: str) -> pd.DataFrame:
    use = df_case_flux[
        (df_case_flux["pathway"] == pathway)
        & (~df_case_flux["pathway"].isin(["Transport, extracellular", "Exchange/demand reaction"]))
        & (df_case_flux["abs_flux"] > FLUX_EPS)
    ].copy()
    if use.empty:
        return pd.DataFrame()

    model = load_species_model(species)
    records = []
    for row in use.sort_values(["abs_flux", "reaction_id"], ascending=[False, True]).itertuples(index=False):
        if row.reaction_id not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(row.reaction_id)
        substrates, products = directed_sides(rxn, float(row.flux))
        records.append(
            {
                "reaction_id": row.reaction_id,
                "reaction_name": rxn.name or row.reaction_name,
                "flux": float(row.flux),
                "abs_flux": float(row.abs_flux),
                "equation_in_flux_direction": f"{side_to_text(substrates)} -> {side_to_text(products)}",
                "substrates": substrates,
                "products": products,
            }
        )
    if not records:
        return pd.DataFrame()

    pipeline = build_stepwise_pipeline(records)
    pipeline.insert(1, "age_group", age_group)
    pipeline.insert(2, "diet", diet)
    pipeline.insert(3, "species", species)
    pipeline.insert(4, "pathway", pathway)
    return pipeline


def build_exchange_summary(df_case_external: pd.DataFrame, medium_exchange_ids: set[str]) -> pd.DataFrame:
    if df_case_external.empty:
        return pd.DataFrame(
            columns=[
                "age_group", "diet", "species", "shared_metabolite_id", "shared_metabolite_name",
                "reaction_id", "reaction_name", "flux", "abs_flux", "direction",
                "medium_available_flag", "exchange_class", "rank_within_species",
            ]
        )

    sub = df_case_external[
        (df_case_external["abs_flux"] > FLUX_EPS)
        & (df_case_external["pathway"].isin(["Transport, extracellular", "Exchange/demand reaction"]))
    ].copy()
    if sub.empty:
        return pd.DataFrame()

    sub["medium_available_flag"] = sub["reaction_id"].astype(str).isin(medium_exchange_ids)
    sub["exchange_class"] = "transport_only"
    sub.loc[(sub["direction"] == "uptake") & (sub["medium_available_flag"]), "exchange_class"] = "diet_uptake"
    sub.loc[(sub["direction"] == "secretion") & (sub["medium_available_flag"]), "exchange_class"] = "diet_secretion"
    sub = sub.sort_values("abs_flux", ascending=False).copy()
    sub["rank_within_species"] = range(1, len(sub) + 1)
    keep = [
        "age_group", "diet", "species", "shared_metabolite_id", "shared_metabolite_name",
        "reaction_id", "reaction_name", "flux", "abs_flux", "direction",
        "medium_available_flag", "exchange_class", "rank_within_species",
    ]
    sub = sub.rename(
        columns={
            "external_metabolite_id": "shared_metabolite_id",
            "external_metabolite_name": "shared_metabolite_name",
        }
    )
    return sub.loc[:, keep].head(TOP_N_EXCHANGES).reset_index(drop=True)


def pathway_keyword_sets(df_case_flux: pd.DataFrame, pathway: str, species: str):
    model = load_species_model(species)
    sub = df_case_flux[
        (df_case_flux["pathway"] == pathway)
        & (~df_case_flux["pathway"].isin(["Transport, extracellular", "Exchange/demand reaction"]))
        & (df_case_flux["abs_flux"] > FLUX_EPS)
    ].copy()
    keyword_set = tokenize(pathway)
    for row in sub.itertuples(index=False):
        keyword_set |= tokenize(row.reaction_name, row.reaction_id)
        if row.reaction_id in model.reactions:
            rxn = model.reactions.get_by_id(row.reaction_id)
            for met in rxn.metabolites:
                keyword_set |= tokenize(met.id, met.name)
    return keyword_set


def rank_anchor_candidates(df_case_external: pd.DataFrame, df_case_flux: pd.DataFrame, pathway: str, medium_exchange_ids: set[str]) -> pd.DataFrame:
    uptakes = df_case_external[
        (df_case_external["pathway"] == "Exchange/demand reaction")
        & (df_case_external["direction"] == "uptake")
        & (df_case_external["abs_flux"] > FLUX_EPS)
        & (~df_case_external["external_metabolite_id"].isin(EXCLUDED_EXTERNAL_METS))
    ].copy()
    if uptakes.empty:
        return pd.DataFrame()

    keyword_set = pathway_keyword_sets(df_case_flux, pathway, str(uptakes["species"].iloc[0]))
    uptakes["text"] = (
        uptakes["external_metabolite_id"].astype(str).str.lower()
        + " "
        + uptakes["external_metabolite_name"].astype(str).str.lower()
        + " "
        + uptakes["reaction_name"].astype(str).str.lower()
    )
    uptakes["keyword_match_score"] = uptakes["text"].apply(
        lambda txt: sum(1 for token in keyword_set if token and token in txt)
    )
    uptakes["medium_available_flag"] = uptakes["reaction_id"].astype(str).isin(medium_exchange_ids)
    return uptakes.sort_values(
        ["keyword_match_score", "abs_flux", "medium_available_flag"],
        ascending=[False, False, False],
    ).reset_index(drop=True)


def build_crossfeeding_candidates(df_external: pd.DataFrame, age_group: str, diet: str, species: str, pathway: str, medium_exchange_ids: set[str]) -> pd.DataFrame:
    sub = df_external[(df_external["age_group"] == age_group) & (df_external["diet"] == diet)].copy()
    focus = sub[
        (sub["species"] == species)
        & (sub["pathway"] == "Exchange/demand reaction")
        & (sub["direction"] == "uptake")
        & (sub["abs_flux"] > FLUX_EPS)
        & (~sub["external_metabolite_id"].isin(EXCLUDED_EXTERNAL_METS))
    ].copy()
    if focus.empty:
        return pd.DataFrame()

    anchor_rank = rank_anchor_candidates(
        df_case_external=sub[sub["species"] == species].copy(),
        df_case_flux=df_flux_case_lookup[(age_group, diet, species)],
        pathway=pathway,
        medium_exchange_ids=medium_exchange_ids,
    )
    anchor_ranks = {}
    if not anchor_rank.empty:
        for i, met_id in enumerate(anchor_rank["external_metabolite_id"].astype(str).tolist(), start=1):
            anchor_ranks[met_id] = i

    rows = []
    for met_row in focus.itertuples(index=False):
        met_id = met_row.external_metabolite_id
        producers = sub[
            (sub["species"] != species)
            & (sub["external_metabolite_id"] == met_id)
            & (sub["pathway"] == "Exchange/demand reaction")
            & (sub["direction"] == "secretion")
            & (sub["abs_flux"] > FLUX_EPS)
        ].sort_values("abs_flux", ascending=False)
        if producers.empty:
            rows.append(
                {
                    "age_group": age_group,
                    "diet": diet,
                    "consumer_species": species,
                    "producer_species": "",
                    "shared_metabolite_id": met_id,
                    "shared_metabolite_name": met_row.external_metabolite_name,
                    "consumer_uptake_flux": float(met_row.abs_flux),
                    "producer_secretion_flux": 0.0,
                    "candidate_score": 0.0,
                    "medium_available_flag": met_row.reaction_id in medium_exchange_ids,
                    "linked_pathway": pathway,
                    "anchor_rank": anchor_ranks.get(met_id, 999),
                    "likely_origin": "direct_from_medium" if met_row.reaction_id in medium_exchange_ids else "internal_only",
                }
            )
            continue
        medium_flag = met_row.reaction_id in medium_exchange_ids
        for prod in producers.itertuples(index=False):
            candidate_score = min(float(met_row.abs_flux), float(prod.abs_flux))
            likely_origin = "crossfed_only"
            if medium_flag:
                likely_origin = "mixed_medium_plus_crossfeeding"
            rows.append(
                {
                    "age_group": age_group,
                    "diet": diet,
                    "consumer_species": species,
                    "producer_species": prod.species,
                    "shared_metabolite_id": met_id,
                    "shared_metabolite_name": met_row.external_metabolite_name,
                    "consumer_uptake_flux": float(met_row.abs_flux),
                    "producer_secretion_flux": float(prod.abs_flux),
                    "candidate_score": candidate_score,
                    "medium_available_flag": medium_flag,
                    "linked_pathway": pathway,
                    "anchor_rank": anchor_ranks.get(met_id, 999),
                    "likely_origin": likely_origin,
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(
        ["anchor_rank", "candidate_score", "consumer_uptake_flux"],
        ascending=[True, False, False],
    ).reset_index(drop=True)


def choose_anchor_metabolite(df_crossfeeding: pd.DataFrame, df_case_external: pd.DataFrame, pathway: str, medium_exchange_ids: set[str]) -> dict:
    if not df_crossfeeding.empty:
        ranked = df_crossfeeding.sort_values(
            ["anchor_rank", "candidate_score", "consumer_uptake_flux"],
            ascending=[True, False, False],
        ).reset_index(drop=True)
        top = ranked.iloc[0]
        return {
            "metabolite_id": top["shared_metabolite_id"],
            "metabolite_name": top["shared_metabolite_name"],
            "uptake_flux": float(top["consumer_uptake_flux"]),
            "crossfeeding_support_score": float(top["candidate_score"]),
            "likely_origin": top["likely_origin"],
            "top_producer_species": top["producer_species"],
            "top_producer_flux": float(top["producer_secretion_flux"]),
        }

    ranked_uptakes = rank_anchor_candidates(df_case_external, df_flux_case_lookup[(df_case_external["age_group"].iloc[0], df_case_external["diet"].iloc[0], df_case_external["species"].iloc[0])], pathway, medium_exchange_ids)
    if not ranked_uptakes.empty:
        top = ranked_uptakes.iloc[0]
        return {
            "metabolite_id": top["external_metabolite_id"],
            "metabolite_name": top["external_metabolite_name"],
            "uptake_flux": float(top["abs_flux"]),
            "crossfeeding_support_score": 0.0,
            "likely_origin": "direct_from_medium" if top["medium_available_flag"] else "internal_only_no_external_anchor",
            "top_producer_species": "",
            "top_producer_flux": 0.0,
        }

    return {
        "metabolite_id": "",
        "metabolite_name": "",
        "uptake_flux": 0.0,
        "crossfeeding_support_score": 0.0,
        "likely_origin": "internal_only_no_external_anchor",
        "top_producer_species": "",
        "top_producer_flux": 0.0,
    }


def build_provenance_edges(df_external: pd.DataFrame, age_group: str, diet: str, species: str, pathway: str, anchor: dict, medium_exchange_ids: set[str]) -> pd.DataFrame:
    rows = []
    met_id = anchor["metabolite_id"]
    if met_id:
        sub = df_external[
            (df_external["age_group"] == age_group)
            & (df_external["diet"] == diet)
            & (df_external["external_metabolite_id"] == met_id)
            & (df_external["pathway"] == "Exchange/demand reaction")
            & (df_external["abs_flux"] > FLUX_EPS)
        ].copy()
        env_rows = sub[
            sub["reaction_id"].astype(str).isin(medium_exchange_ids)
            & (sub["direction"] == "uptake")
        ]
        for row in env_rows.itertuples(index=False):
            rows.append(
                {
                    "age_group": age_group,
                    "diet": diet,
                    "species": species,
                    "pathway": pathway,
                    "metabolite_id": met_id,
                    "metabolite_name": anchor["metabolite_name"],
                    "layer_from": "environment",
                    "layer_to": "organism_exchange",
                    "species_from": "",
                    "species_to": row.species,
                    "reaction_id": row.reaction_id,
                    "reaction_name": row.reaction_name,
                    "flux": float(row.flux),
                    "abs_flux": float(row.abs_flux),
                    "role": "environment_exchange",
                }
            )
        producers = sub[(sub["species"] != species) & (sub["direction"] == "secretion")]
        consumers = sub[(sub["species"] == species) & (sub["direction"] == "uptake")]
        for row in producers.itertuples(index=False):
            rows.append(
                {
                    "age_group": age_group,
                    "diet": diet,
                    "species": species,
                    "pathway": pathway,
                    "metabolite_id": met_id,
                    "metabolite_name": anchor["metabolite_name"],
                    "layer_from": "organism_exchange",
                    "layer_to": "shared_metabolite",
                    "species_from": row.species,
                    "species_to": "",
                    "reaction_id": row.reaction_id,
                    "reaction_name": row.reaction_name,
                    "flux": float(row.flux),
                    "abs_flux": float(row.abs_flux),
                    "role": "organism_secretion",
                }
            )
        for row in consumers.itertuples(index=False):
            rows.append(
                {
                    "age_group": age_group,
                    "diet": diet,
                    "species": species,
                    "pathway": pathway,
                    "metabolite_id": met_id,
                    "metabolite_name": anchor["metabolite_name"],
                    "layer_from": "shared_metabolite",
                    "layer_to": "organism_exchange",
                    "species_from": "",
                    "species_to": row.species,
                    "reaction_id": row.reaction_id,
                    "reaction_name": row.reaction_name,
                    "flux": float(row.flux),
                    "abs_flux": float(row.abs_flux),
                    "role": "organism_uptake",
                }
            )
        for prod in producers.itertuples(index=False):
            for cons in consumers.itertuples(index=False):
                rows.append(
                    {
                        "age_group": age_group,
                        "diet": diet,
                        "species": species,
                        "pathway": pathway,
                        "metabolite_id": met_id,
                        "metabolite_name": anchor["metabolite_name"],
                        "layer_from": "organism_exchange",
                        "layer_to": "organism_exchange",
                        "species_from": prod.species,
                        "species_to": cons.species,
                        "reaction_id": "",
                        "reaction_name": f"Cross-feeding candidate via {anchor['metabolite_name'] or met_id}",
                        "flux": min(float(prod.abs_flux), float(cons.abs_flux)),
                        "abs_flux": min(float(prod.abs_flux), float(cons.abs_flux)),
                        "role": "crossfeeding_candidate",
                    }
                )

    pathway_rows = df_flux_case_lookup[(age_group, diet, species)]
    internal_hits = pathway_rows[
        (pathway_rows["pathway"] == pathway)
        & (~pathway_rows["pathway"].isin(["Transport, extracellular", "Exchange/demand reaction"]))
        & (pathway_rows["abs_flux"] > FLUX_EPS)
    ].sort_values("abs_flux", ascending=False)
    for row in internal_hits.itertuples(index=False):
        rows.append(
            {
                "age_group": age_group,
                "diet": diet,
                "species": species,
                "pathway": pathway,
                "metabolite_id": met_id,
                "metabolite_name": anchor["metabolite_name"],
                "layer_from": "organism_exchange",
                "layer_to": "species_internal",
                "species_from": species,
                "species_to": species,
                "reaction_id": row.reaction_id,
                "reaction_name": row.reaction_name,
                "flux": float(row.flux),
                "abs_flux": float(row.abs_flux),
                "role": "internal_processing",
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    role_order = {
        "environment_exchange": 1,
        "organism_secretion": 2,
        "crossfeeding_candidate": 3,
        "organism_uptake": 4,
        "internal_processing": 5,
    }
    out["role_order"] = out["role"].map(role_order).fillna(99)
    out = out.sort_values(["role_order", "abs_flux"], ascending=[True, False]).drop(columns=["role_order"])
    return out.reset_index(drop=True)


def write_summary_txt(out_path: Path, age_group: str, diet: str, species: str, pathway: str, growth_rate: float, anchor: dict, crossfeeding: pd.DataFrame):
    lines = [
        f"MICOM provenance summary",
        f"Tradeoff fraction: {DEFAULT_TRADEOFF}",
        f"Age group: {age_group}",
        f"Diet: {diet}",
        f"Species: {species}",
        f"Top pathway: {pathway}",
        f"Species growth rate: {growth_rate:.6f}",
    ]
    if anchor["metabolite_id"]:
        lines.extend(
            [
                f"Anchor metabolite: {anchor['metabolite_id']} ({anchor['metabolite_name']})",
                f"Anchor uptake flux: {anchor['uptake_flux']:.6f}",
                f"Likely origin: {anchor['likely_origin']}",
                f"Crossfeeding support score: {anchor['crossfeeding_support_score']:.6f}",
                f"Top producer species: {anchor['top_producer_species'] or 'none'}",
                f"Top producer flux: {anchor['top_producer_flux']:.6f}",
            ]
        )
    else:
        lines.append("Anchor metabolite: none")
        lines.append("Likely origin: internal_only_no_external_anchor")

    if not crossfeeding.empty:
        lines.append("")
        lines.append("Top crossfeeding candidates:")
        for row in crossfeeding.head(5).itertuples(index=False):
            lines.append(
                f"- {row.producer_species or 'none'} -> {row.shared_metabolite_id} -> {row.consumer_species}; "
                f"support={row.candidate_score:.6f}; origin={row.likely_origin}"
            )
    out_path.write_text("\n".join(lines) + "\n")


def case_prefix(age_group: str, diet: str, species: str) -> str:
    return f"{slug(age_group)}__{slug(diet)}__{slug(species)}"


def parse_args():
    parser = argparse.ArgumentParser(description="Build MICOM provenance summaries for selected age bins, diets, and species.")
    parser.add_argument("--age-groups", nargs="+", default=DEFAULT_AGE_GROUPS)
    parser.add_argument("--diets", nargs="+", default=DEFAULT_DIETS)
    parser.add_argument("--species", nargs="+", default=DEFAULT_SPECIES)
    parser.add_argument("--tradeoff-fraction", type=float, default=DEFAULT_TRADEOFF)
    return parser.parse_args()


df_flux_case_lookup = {}


def main():
    global df_flux_case_lookup
    args = parse_args()

    if not IN_FLUX.exists():
        raise FileNotFoundError(f"Missing MICOM reaction flux input: {IN_FLUX}")
    if not IN_GROWTH.exists():
        raise FileNotFoundError(f"Missing MICOM growth input: {IN_GROWTH}")
    if not IN_TRADEOFF.exists():
        raise FileNotFoundError(f"Missing MICOM tradeoff summary input: {IN_TRADEOFF}")

    df_flux = pd.read_csv(IN_FLUX)
    df_growth = pd.read_csv(IN_GROWTH)
    df_tradeoff = pd.read_csv(IN_TRADEOFF)
    validate_tradeoff_context(df_tradeoff, args.age_groups, args.diets, args.tradeoff_fraction)

    required_flux = {"id", "reaction_id", "flux", "age_group", "diet", "abs_flux", "reaction_name", "pathway"}
    required_growth = {"age_group", "diet", "id", "growth_rate", "grows_in_diet"}
    if not required_flux.issubset(df_flux.columns):
        raise ValueError(f"MICOM flux table must contain columns: {sorted(required_flux)}")
    if not required_growth.issubset(df_growth.columns):
        raise ValueError(f"MICOM growth table must contain columns: {sorted(required_growth)}")

    df_flux = df_flux[
        (df_flux["age_group"].isin(args.age_groups))
        & (df_flux["diet"].isin(args.diets))
    ].copy()
    df_growth = df_growth[
        (df_growth["age_group"].isin(args.age_groups))
        & (df_growth["diet"].isin(args.diets))
        & (df_growth["id"].isin(args.species))
        & (df_growth["grows_in_diet"])
    ].copy()

    media_lookup = {
        "western": load_medium(WESTERN_CSV),
        "high_fiber": load_medium(FIBER_CSV),
    }
    df_external = build_external_flux_table(df_flux)

    df_flux_focus = df_flux[df_flux["id"].isin(args.species)].copy()
    for key, sub in df_flux_focus.groupby(["age_group", "diet", "id"]):
        df_flux_case_lookup[key] = sub.copy()

    case_rows = []
    for row in df_growth.sort_values(["age_group", "diet", "id"]).itertuples(index=False):
        key = (row.age_group, row.diet, row.id)
        if key not in df_flux_case_lookup:
            continue
        df_case_flux = df_flux_case_lookup[key]
        df_case_external = df_external[
            (df_external["age_group"] == row.age_group)
            & (df_external["diet"] == row.diet)
            & (df_external["species"] == row.id)
        ].copy()
        medium_exchange_ids = media_lookup[row.diet]
        pathway_overview = build_pathway_overview(df_case_flux, float(row.growth_rate))
        top_pathway = pathway_overview["pathway"].iloc[0] if not pathway_overview.empty else ""
        reaction_steps = trace_internal_pathway(df_case_flux, row.id, row.age_group, row.diet, top_pathway) if top_pathway else pd.DataFrame()
        exchange_summary = build_exchange_summary(df_case_external, medium_exchange_ids)
        crossfeeding = build_crossfeeding_candidates(df_external, row.age_group, row.diet, row.id, top_pathway, medium_exchange_ids) if top_pathway else pd.DataFrame()
        anchor = choose_anchor_metabolite(crossfeeding, df_case_external, top_pathway, medium_exchange_ids) if top_pathway else {
            "metabolite_id": "",
            "metabolite_name": "",
            "uptake_flux": 0.0,
            "crossfeeding_support_score": 0.0,
            "likely_origin": "internal_only_no_external_anchor",
            "top_producer_species": "",
            "top_producer_flux": 0.0,
        }
        provenance_edges = build_provenance_edges(df_external, row.age_group, row.diet, row.id, top_pathway, anchor, medium_exchange_ids) if top_pathway else pd.DataFrame()

        prefix = case_prefix(row.age_group, row.diet, row.id)
        pathway_overview.to_csv(OUT_DIR / f"{prefix}__pathway_overview.csv", index=False)
        reaction_steps.to_csv(OUT_DIR / f"{prefix}__reaction_steps.csv", index=False)
        exchange_summary.to_csv(OUT_DIR / f"{prefix}__exchange_summary.csv", index=False)
        crossfeeding.to_csv(OUT_DIR / f"{prefix}__crossfeeding_candidates.csv", index=False)
        provenance_edges.to_csv(OUT_DIR / f"{prefix}__provenance_edges.csv", index=False)
        write_summary_txt(
            OUT_DIR / f"{prefix}__summary.txt",
            age_group=row.age_group,
            diet=row.diet,
            species=row.id,
            pathway=top_pathway,
            growth_rate=float(row.growth_rate),
            anchor=anchor,
            crossfeeding=crossfeeding,
        )

        top_score = float(pathway_overview["biomass_association_score"].iloc[0]) if not pathway_overview.empty else 0.0
        case_rows.append(
            {
                "tradeoff_fraction": args.tradeoff_fraction,
                "age_group": row.age_group,
                "diet": row.diet,
                "species": row.id,
                "growth_rate": float(row.growth_rate),
                "top_pathway": top_pathway,
                "top_pathway_score": top_score,
                "anchor_metabolite_id": anchor["metabolite_id"],
                "anchor_metabolite_name": anchor["metabolite_name"],
                "anchor_uptake_flux": anchor["uptake_flux"],
                "crossfeeding_support_score": anchor["crossfeeding_support_score"],
                "likely_origin": anchor["likely_origin"],
                "top_producer_species": anchor["top_producer_species"],
                "top_producer_flux": anchor["top_producer_flux"],
            }
        )
        print(f"Saved MICOM provenance case: {prefix}")

    if case_rows:
        pd.DataFrame(case_rows).sort_values(["species", "age_group", "diet"]).to_csv(OUT_CASE_OVERVIEW, index=False)
        print(f"Saved case overview: {OUT_CASE_OVERVIEW}")


if __name__ == "__main__":
    main()
