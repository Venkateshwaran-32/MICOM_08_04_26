from pathlib import Path
import argparse
import re
from collections import Counter, defaultdict
from typing import Optional

import cobra
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODELS_DIR = PROJECT_ROOT / "Models" / "vmh_agora_sbml"
FBA_DIR = PROJECT_ROOT / "Results" / "fba"
OUT_DIR = FBA_DIR / "provenance"
OUT_DIR.mkdir(parents=True, exist_ok=True)

IN_COMM_EX = FBA_DIR / "community_exchange_fluxes_by_diet.csv"
IN_CONN = FBA_DIR / "community_species_connector_fluxes_by_diet.csv"
IN_RXN = FBA_DIR / "community_reaction_fluxes_by_diet.csv"
IN_PATHWAY = FBA_DIR / "community_pathway_activity_by_diet.csv"
IN_BIOMASS = FBA_DIR / "community_biomass_associated_pathways_by_diet.csv"

FLUX_EPS = 1e-9
COMMON_COFATOR_METS = {
    "h[c]", "h[e]", "h2o[c]", "h2o[e]", "pi[c]", "pi[e]", "ppi[c]", "ppi[e]",
    "atp[c]", "atp[e]", "adp[c]", "adp[e]", "amp[c]", "amp[e]",
    "nad[c]", "nad[e]", "nadh[c]", "nadh[e]", "nadp[c]", "nadp[e]", "nadph[c]", "nadph[e]",
    "coa[c]", "coa[e]", "co2[c]", "co2[e]", "o2[c]", "o2[e]",
}


def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", str(text))


def slug(text: str) -> str:
    return sanitize(text).strip("_")


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


def make_connector_summary(df_conn: pd.DataFrame, diet: str, species: str, top_n: int) -> pd.DataFrame:
    sub = df_conn[(df_conn["diet"] == diet) & (df_conn["species"] == species)].copy()
    if sub.empty:
        return pd.DataFrame(
            columns=[
                "diet", "species", "shared_metabolite_id", "shared_metabolite_name",
                "connector_reaction_id", "flux", "abs_flux", "direction", "rank_within_species",
            ]
        )

    sub["direction"] = "mixed"
    sub.loc[sub["uptake_from_shared_flux"] > FLUX_EPS, "direction"] = "uptake"
    sub.loc[sub["secretion_to_shared_flux"] > FLUX_EPS, "direction"] = "secretion"
    sub = sub.sort_values("abs_flux", ascending=False).copy()
    sub["rank_within_species"] = range(1, len(sub) + 1)
    sub["shared_metabolite_id"] = sub["external_metabolite_id"]
    sub["shared_metabolite_name"] = sub["external_metabolite_name"]
    sub["connector_reaction_id"] = sub["connector_reaction"]
    keep = [
        "diet", "species", "shared_metabolite_id", "shared_metabolite_name",
        "connector_reaction_id", "flux", "abs_flux", "direction", "rank_within_species",
    ]
    return sub.loc[:, keep].head(top_n).reset_index(drop=True)


def make_crossfeeding_candidates(df_conn: pd.DataFrame, diet: str, focus_species: str) -> pd.DataFrame:
    sub = df_conn[df_conn["diet"] == diet].copy()
    rows = []
    for metabolite_id, met_df in sub.groupby("external_metabolite_id", dropna=False):
        producers = met_df[met_df["secretion_to_shared_flux"] > FLUX_EPS]
        consumers = met_df[met_df["uptake_from_shared_flux"] > FLUX_EPS]
        if producers.empty or consumers.empty:
            continue
        met_name = met_df["external_metabolite_name"].dropna().astype(str).iloc[0] if met_df["external_metabolite_name"].notna().any() else ""
        for prod in producers.itertuples(index=False):
            for cons in consumers.itertuples(index=False):
                if prod.species == cons.species:
                    continue
                rows.append(
                    {
                        "diet": diet,
                        "shared_metabolite_id": metabolite_id,
                        "shared_metabolite_name": met_name,
                        "producer_species": prod.species,
                        "consumer_species": cons.species,
                        "producer_flux": float(prod.secretion_to_shared_flux),
                        "consumer_flux": float(cons.uptake_from_shared_flux),
                        "candidate_score": min(float(prod.secretion_to_shared_flux), float(cons.uptake_from_shared_flux)),
                        "crossfeeding_candidate": True,
                        "involves_focus_species": prod.species == focus_species or cons.species == focus_species,
                    }
                )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values(
        ["involves_focus_species", "candidate_score", "shared_metabolite_id"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    return out


def build_species_biomass_support_summary(
    df_conn: pd.DataFrame,
    df_biomass: pd.DataFrame,
    df_rxn: pd.DataFrame,
    species: str,
    diet: str,
    top_n_metabolites: int = 12,
    top_n_pathways: int = 6,
) -> pd.DataFrame:
    conn_sub = df_conn[
        (df_conn["species"] == species)
        & (df_conn["diet"] == diet)
        & (df_conn["uptake_from_shared_flux"] > FLUX_EPS)
    ].copy()
    if conn_sub.empty:
        return pd.DataFrame(
            columns=[
                "diet", "species", "shared_metabolite_id", "shared_metabolite_name",
                "uptake_from_shared_flux", "top_producer_species", "top_producer_flux",
                "crossfeeding_support_score", "linked_pathway", "pathway_abs_flux",
                "pathway_flux_share", "biomass_association_score", "link_evidence",
            ]
        )

    top_metabolites = conn_sub.sort_values("uptake_from_shared_flux", ascending=False).head(top_n_metabolites).copy()
    biomass_sub = df_biomass[
        (df_biomass["species"] == species)
        & (df_biomass["diet"] == diet)
    ].copy()
    biomass_sub = biomass_sub.sort_values("biomass_association_score", ascending=False)
    preferred_pathways = biomass_sub[
        ~biomass_sub["pathway"].isin(["Transport, extracellular", "Exchange/demand reaction"])
    ].head(top_n_pathways).copy()
    if preferred_pathways.empty:
        preferred_pathways = biomass_sub.head(top_n_pathways).copy()

    species_rxn = df_rxn[
        (df_rxn["species"] == species)
        & (df_rxn["diet"] == diet)
        & (df_rxn["reaction_type"] == "species_internal")
        & (df_rxn["abs_flux"] > FLUX_EPS)
    ].copy()

    pathway_keywords = {
        "Nucleotide interconversion": ["adenosine", "guanosine", "cytidine", "uridine", "nucleoside", "purine", "pyrimidine", "adenyl"],
        "Pentose phosphate pathway": ["ribose", "r1p", "glyceraldehyde", "sedoheptulose", "pentose"],
        "Glycolysis/gluconeogenesis": ["glucose", "fructose", "glyceraldehyde", "pyruvate", "lactate"],
        "Pyruvate metabolism": ["pyruvate", "acetate", "lactate", "formate", "alanine"],
        "Alanine and aspartate metabolism": ["alanine", "aspartate", "glutamate", "asparagine"],
        "Valine, leucine, and isoleucine metabolism": ["valine", "leucine", "isoleucine", "2-oxobutanoate", "2obut"],
        "Glycerophospholipid metabolism": ["glycerol", "glycerol 3-phosphate", "glyc3p"],
    }

    def normalize_text(*parts) -> str:
        return " ".join(str(x).lower() for x in parts if pd.notna(x))

    summary_rows = []
    for met_row in top_metabolites.itertuples(index=False):
        met_text = normalize_text(met_row.external_metabolite_id, met_row.external_metabolite_name)
        met_all = df_conn[
            (df_conn["diet"] == diet)
            & (df_conn["external_metabolite_id"] == met_row.external_metabolite_id)
        ].copy()
        producers = met_all[met_all["secretion_to_shared_flux"] > FLUX_EPS].sort_values(
            "secretion_to_shared_flux", ascending=False
        )
        top_producer_species = producers.iloc[0]["species"] if not producers.empty else ""
        top_producer_flux = float(producers.iloc[0]["secretion_to_shared_flux"]) if not producers.empty else 0.0
        crossfeeding_support_score = min(float(met_row.uptake_from_shared_flux), top_producer_flux) if top_producer_flux > 0 else 0.0

        linked_any = False
        for path_row in preferred_pathways.itertuples(index=False):
            pathway = path_row.pathway
            keywords = pathway_keywords.get(pathway, [])
            keyword_match = any(k in met_text for k in keywords)
            pathway_rxn = species_rxn[species_rxn["pathway"] == pathway].copy()
            rxn_match = pathway_rxn[
                pathway_rxn["reaction_name"].astype(str).str.contains("|".join(map(re.escape, keywords)), case=False, na=False)
            ] if keywords else pathway_rxn.iloc[0:0]
            if not keyword_match and rxn_match.empty:
                continue

            linked_any = True
            evidence_bits = []
            if keyword_match:
                evidence_bits.append("metabolite_name_keyword")
            if not rxn_match.empty:
                evidence_bits.append("pathway_reaction_keyword")
            summary_rows.append(
                {
                    "diet": diet,
                    "species": species,
                    "shared_metabolite_id": met_row.external_metabolite_id,
                    "shared_metabolite_name": met_row.external_metabolite_name,
                    "uptake_from_shared_flux": float(met_row.uptake_from_shared_flux),
                    "top_producer_species": top_producer_species,
                    "top_producer_flux": top_producer_flux,
                    "crossfeeding_support_score": crossfeeding_support_score,
                    "linked_pathway": pathway,
                    "pathway_abs_flux": float(path_row.pathway_abs_flux),
                    "pathway_flux_share": float(path_row.pathway_flux_share),
                    "biomass_association_score": float(path_row.biomass_association_score),
                    "link_evidence": " + ".join(evidence_bits),
                }
            )

        if not linked_any:
            fallback = preferred_pathways.iloc[0]
            summary_rows.append(
                {
                    "diet": diet,
                    "species": species,
                    "shared_metabolite_id": met_row.external_metabolite_id,
                    "shared_metabolite_name": met_row.external_metabolite_name,
                    "uptake_from_shared_flux": float(met_row.uptake_from_shared_flux),
                    "top_producer_species": top_producer_species,
                    "top_producer_flux": top_producer_flux,
                    "crossfeeding_support_score": crossfeeding_support_score,
                    "linked_pathway": fallback["pathway"],
                    "pathway_abs_flux": float(fallback["pathway_abs_flux"]),
                    "pathway_flux_share": float(fallback["pathway_flux_share"]),
                    "biomass_association_score": float(fallback["biomass_association_score"]),
                    "link_evidence": "fallback_top_biomass_pathway",
                }
            )

    out = pd.DataFrame(summary_rows)
    if out.empty:
        return out
    out = out.sort_values(
        ["crossfeeding_support_score", "uptake_from_shared_flux", "biomass_association_score"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    return out


def build_clean_crossfeeding_biomass_summary(
    df_conn: pd.DataFrame,
    df_biomass: pd.DataFrame,
    species: str,
    diet: str,
    top_n_metabolites: int = 20,
) -> pd.DataFrame:
    conn_sub = df_conn[
        (df_conn["species"] == species)
        & (df_conn["diet"] == diet)
        & (df_conn["uptake_from_shared_flux"] > FLUX_EPS)
    ].copy()
    if conn_sub.empty:
        return pd.DataFrame()

    # Keep only metabolites with interpretable biochemical classes for a presentation-safe summary.
    explicit_map = {
        "adn[e]": ("nucleoside", "Nucleotide interconversion"),
        "gsn[e]": ("nucleoside", "Nucleotide interconversion"),
        "cytd[e]": ("nucleoside", "Nucleotide interconversion"),
        "uri[e]": ("nucleoside", "Nucleotide interconversion"),
        "dcyt[e]": ("nucleoside", "Nucleotide interconversion"),
        "asp_L[e]": ("amino_acid", "Alanine and aspartate metabolism"),
        "glu_L[e]": ("amino_acid", "Alanine and aspartate metabolism"),
        "gln_L[e]": ("amino_acid", "Alanine and aspartate metabolism"),
        "ser_L[e]": ("amino_acid", "Alanine and aspartate metabolism"),
        "ala_L[e]": ("amino_acid", "Alanine and aspartate metabolism"),
        "val_L[e]": ("branched_chain_amino_acid", "Valine, leucine, and isoleucine metabolism"),
        "leu_L[e]": ("branched_chain_amino_acid", "Valine, leucine, and isoleucine metabolism"),
        "ile_L[e]": ("branched_chain_amino_acid", "Valine, leucine, and isoleucine metabolism"),
        "2obut[e]": ("branched_chain_keto_acid", "Valine, leucine, and isoleucine metabolism"),
        "glyc3p[e]": ("glycerophospholipid_related", "Glycerophospholipid metabolism"),
    }
    excluded = {"h[e]", "h2o[e]", "co2[e]", "pi[e]", "nh4[e]"}

    biomass_lookup = (
        df_biomass[(df_biomass["species"] == species) & (df_biomass["diet"] == diet)]
        .set_index("pathway")
        .to_dict("index")
    )

    rows = []
    for met_row in conn_sub.sort_values("uptake_from_shared_flux", ascending=False).head(top_n_metabolites).itertuples(index=False):
        met_id = met_row.external_metabolite_id
        if met_id in excluded or met_id not in explicit_map:
            continue

        metabolite_class, linked_pathway = explicit_map[met_id]
        biomass_info = biomass_lookup.get(linked_pathway)
        if not biomass_info:
            continue

        met_all = df_conn[
            (df_conn["diet"] == diet)
            & (df_conn["external_metabolite_id"] == met_id)
        ].copy()
        producers = met_all[met_all["secretion_to_shared_flux"] > FLUX_EPS].sort_values(
            "secretion_to_shared_flux", ascending=False
        )
        if producers.empty:
            continue

        top_producer = producers.iloc[0]
        support_score = min(float(met_row.uptake_from_shared_flux), float(top_producer["secretion_to_shared_flux"]))
        confidence = "medium"
        if metabolite_class == "nucleoside" and support_score >= 5:
            confidence = "high"
        elif support_score < 1:
            confidence = "low"

        rows.append(
            {
                "diet": diet,
                "producer_species": top_producer["species"],
                "shared_metabolite_id": met_id,
                "shared_metabolite_name": met_row.external_metabolite_name,
                "consumer_species": species,
                "uptake_flux": float(met_row.uptake_from_shared_flux),
                "producer_flux": float(top_producer["secretion_to_shared_flux"]),
                "crossfeeding_support_score": support_score,
                "metabolite_class": metabolite_class,
                "consumer_top_pathway": linked_pathway,
                "pathway_flux_share": float(biomass_info["pathway_flux_share"]),
                "biomass_association_score": float(biomass_info["biomass_association_score"]),
                "evidence_type": "connector_flux + pathway_class_match",
                "confidence": confidence,
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values(
        ["confidence", "crossfeeding_support_score", "biomass_association_score"],
        ascending=[True, False, False],
    ).copy()
    confidence_order = {"high": 0, "medium": 1, "low": 2}
    out["confidence_order"] = out["confidence"].map(confidence_order)
    out = out.sort_values(
        ["confidence_order", "crossfeeding_support_score", "biomass_association_score"],
        ascending=[True, False, False],
    ).drop(columns=["confidence_order"])
    return out.reset_index(drop=True)


def build_pathway_overview(
    df_biomass: pd.DataFrame,
    df_rxn: pd.DataFrame,
    species: str,
    diet: str,
    top_n_pathways: int = 6,
) -> pd.DataFrame:
    biomass_sub = df_biomass[
        (df_biomass["species"] == species)
        & (df_biomass["diet"] == diet)
    ].copy()
    biomass_sub = biomass_sub.sort_values("biomass_association_score", ascending=False)
    biomass_sub = biomass_sub[
        ~biomass_sub["pathway"].isin(["Transport, extracellular", "Exchange/demand reaction"])
    ].head(top_n_pathways).copy()
    if biomass_sub.empty:
        return pd.DataFrame()

    rxn_sub = df_rxn[
        (df_rxn["species"] == species)
        & (df_rxn["diet"] == diet)
        & (df_rxn["reaction_type"] == "species_internal")
        & (df_rxn["abs_flux"] > FLUX_EPS)
    ].copy()
    top_rxn = (
        rxn_sub.sort_values(["pathway", "abs_flux"], ascending=[True, False])
        .groupby("pathway", as_index=False)
        .first()[["pathway", "local_reaction_id", "reaction_name", "abs_flux"]]
        .rename(columns={"local_reaction_id": "top_reaction_id", "reaction_name": "top_reaction_name", "abs_flux": "top_reaction_abs_flux"})
    )
    out = biomass_sub.merge(top_rxn, on="pathway", how="left")
    return out.reset_index(drop=True)


def pick_pathway(
    df_rxn: pd.DataFrame,
    df_biomass: pd.DataFrame,
    species: str,
    diet: str,
    metabolite_id: str,
    requested_pathway: Optional[str],
) -> str:
    if requested_pathway:
        return requested_pathway

    species_rxn = df_rxn[
        (df_rxn["species"] == species)
        & (df_rxn["diet"] == diet)
        & (df_rxn["reaction_type"] == "species_internal")
        & (df_rxn["abs_flux"] > FLUX_EPS)
    ].copy()
    if species_rxn.empty:
        raise ValueError("No active species_internal reactions found for the requested species and diet.")

    transport_like = {"Transport, extracellular", "Exchange/demand reaction"}
    hint_mask = species_rxn["reaction_name"].astype(str).str.contains(
        "adenosine|adenyl|purine|nucleoside", case=False, na=False
    )
    hinted = [
        pathway for pathway in species_rxn.loc[hint_mask, "pathway"].dropna().astype(str).tolist()
        if pathway not in transport_like
    ]
    if hinted:
        ranked = (
            pd.Series(hinted)
            .value_counts()
            .sort_values(ascending=False)
            .index
            .tolist()
        )
        hinted = ranked

    biomass_sub = df_biomass[(df_biomass["species"] == species) & (df_biomass["diet"] == diet)].copy()
    if biomass_sub.empty:
        if hinted:
            return hinted[0]
        return species_rxn.sort_values("abs_flux", ascending=False)["pathway"].iloc[0]

    biomass_sub = biomass_sub.sort_values("biomass_association_score", ascending=False)
    for pathway in biomass_sub["pathway"].astype(str):
        if pathway in hinted:
            return pathway

    if hinted:
        return hinted[0]
    non_transport = biomass_sub[~biomass_sub["pathway"].isin(transport_like)]
    if not non_transport.empty:
        return non_transport["pathway"].iloc[0]
    return biomass_sub["pathway"].iloc[0]


def trace_internal_pathway(
    df_rxn: pd.DataFrame,
    df_biomass: pd.DataFrame,
    species: str,
    diet: str,
    pathway: str,
) -> pd.DataFrame:
    use = df_rxn[
        (df_rxn["diet"] == diet)
        & (df_rxn["species"] == species)
        & (df_rxn["pathway"] == pathway)
        & (df_rxn["reaction_type"] == "species_internal")
        & (df_rxn["abs_flux"] > FLUX_EPS)
    ].copy()
    if use.empty:
        raise ValueError("No active species_internal reactions matched the chosen pathway.")

    model = load_species_model(species)
    records = []
    for row in use.sort_values(["abs_flux", "local_reaction_id"], ascending=[False, True]).itertuples(index=False):
        if row.local_reaction_id not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(row.local_reaction_id)
        substrates, products = directed_sides(rxn, float(row.flux))
        records.append(
            {
                "reaction_id": row.local_reaction_id,
                "reaction_name": rxn.name,
                "flux": float(row.flux),
                "abs_flux": float(row.abs_flux),
                "equation_in_flux_direction": f"{side_to_text(substrates)} -> {side_to_text(products)}",
                "substrates": substrates,
                "products": products,
            }
        )
    if not records:
        raise ValueError("Matched pathway reactions but could not map them back to the SBML model.")

    pipeline = build_stepwise_pipeline(records)
    biomass_adjacent = not df_biomass[
        (df_biomass["species"] == species) & (df_biomass["diet"] == diet) & (df_biomass["pathway"] == pathway)
    ].empty
    pipeline.insert(1, "diet", diet)
    pipeline.insert(2, "species", species)
    pipeline.insert(3, "pathway", pathway)
    pipeline["biomass_adjacent_flag"] = biomass_adjacent
    return pipeline


def build_provenance_edges(
    df_ex: pd.DataFrame,
    df_conn: pd.DataFrame,
    df_rxn: pd.DataFrame,
    df_biomass: pd.DataFrame,
    species: str,
    diet: str,
    metabolite_id: str,
    pathway: str,
) -> pd.DataFrame:
    rows = []
    metabolite_name = ""

    conn_met = df_conn[(df_conn["diet"] == diet) & (df_conn["external_metabolite_id"] == metabolite_id)].copy()
    if not conn_met.empty and conn_met["external_metabolite_name"].notna().any():
        metabolite_name = conn_met["external_metabolite_name"].dropna().astype(str).iloc[0]

    sanitized_ext = f"u__{sanitize(metabolite_id)}"
    ex_match = df_ex[(df_ex["diet"] == diet) & (df_ex["metabolite_id"] == sanitized_ext)].copy()
    for row in ex_match.itertuples(index=False):
        rows.append(
            {
                "diet": diet,
                "metabolite_id": metabolite_id,
                "metabolite_name": metabolite_name,
                "layer_from": "environment",
                "layer_to": "shared_pool",
                "species_from": "",
                "species_to": "",
                "reaction_id": row.exchange_id,
                "reaction_name": f"Community exchange for {metabolite_name or metabolite_id}",
                "pathway": "Exchange/demand reaction",
                "flux": float(row.flux),
                "abs_flux": float(row.abs_flux),
                "role": "environment_exchange",
            }
        )

    producers = conn_met[conn_met["secretion_to_shared_flux"] > FLUX_EPS]
    for row in producers.itertuples(index=False):
        rows.append(
            {
                "diet": diet,
                "metabolite_id": metabolite_id,
                "metabolite_name": metabolite_name,
                "layer_from": "species_connector",
                "layer_to": "shared_pool",
                "species_from": row.species,
                "species_to": "",
                "reaction_id": row.connector_reaction,
                "reaction_name": f"Connector secretion of {metabolite_name or metabolite_id}",
                "pathway": "Shared-pool connector",
                "flux": float(row.flux),
                "abs_flux": float(row.abs_flux),
                "role": "shared_pool_secretion",
            }
        )

    consumers = conn_met[conn_met["uptake_from_shared_flux"] > FLUX_EPS]
    for row in consumers.itertuples(index=False):
        rows.append(
            {
                "diet": diet,
                "metabolite_id": metabolite_id,
                "metabolite_name": metabolite_name,
                "layer_from": "shared_pool",
                "layer_to": "species_connector",
                "species_from": "",
                "species_to": row.species,
                "reaction_id": row.connector_reaction,
                "reaction_name": f"Connector uptake of {metabolite_name or metabolite_id}",
                "pathway": "Shared-pool connector",
                "flux": float(row.flux),
                "abs_flux": float(row.abs_flux),
                "role": "shared_pool_uptake",
            }
        )

    if not producers.empty and not consumers.empty:
        for prod in producers.itertuples(index=False):
            for cons in consumers.itertuples(index=False):
                if prod.species == cons.species:
                    continue
                rows.append(
                    {
                        "diet": diet,
                        "metabolite_id": metabolite_id,
                        "metabolite_name": metabolite_name,
                        "layer_from": "species_connector",
                        "layer_to": "species_connector",
                        "species_from": prod.species,
                        "species_to": cons.species,
                        "reaction_id": "",
                        "reaction_name": f"Cross-feeding candidate via {metabolite_name or metabolite_id}",
                        "pathway": "Shared-pool connector",
                        "flux": min(float(prod.secretion_to_shared_flux), float(cons.uptake_from_shared_flux)),
                        "abs_flux": min(float(prod.secretion_to_shared_flux), float(cons.uptake_from_shared_flux)),
                        "role": "crossfeeding_candidate",
                    }
                )

    species_rxn = df_rxn[
        (df_rxn["diet"] == diet)
        & (df_rxn["species"] == species)
        & (df_rxn["reaction_type"] == "species_internal")
        & (df_rxn["abs_flux"] > FLUX_EPS)
    ].copy()
    internal_mask = species_rxn["reaction_name"].astype(str).str.contains("adenosine|adenyl|purine|nucleoside", case=False, na=False)
    internal_hits = species_rxn[internal_mask].sort_values("abs_flux", ascending=False)
    for row in internal_hits.itertuples(index=False):
        rows.append(
            {
                "diet": diet,
                "metabolite_id": metabolite_id,
                "metabolite_name": metabolite_name,
                "layer_from": "species_connector",
                "layer_to": "species_internal",
                "species_from": species,
                "species_to": species,
                "reaction_id": row.local_reaction_id,
                "reaction_name": row.reaction_name,
                "pathway": row.pathway,
                "flux": float(row.flux),
                "abs_flux": float(row.abs_flux),
                "role": "internal_processing",
            }
        )

    biomass_match = df_biomass[
        (df_biomass["species"] == species)
        & (df_biomass["diet"] == diet)
        & (df_biomass["pathway"] == pathway)
    ].copy()
    for row in biomass_match.itertuples(index=False):
        rows.append(
            {
                "diet": diet,
                "metabolite_id": metabolite_id,
                "metabolite_name": metabolite_name,
                "layer_from": "species_internal",
                "layer_to": "biomass_context",
                "species_from": species,
                "species_to": species,
                "reaction_id": "",
                "reaction_name": f"Biomass-adjacent pathway summary for {pathway}",
                "pathway": pathway,
                "flux": float(row.biomass_association_score),
                "abs_flux": float(row.biomass_association_score),
                "role": "biomass_adjacent",
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    role_order = {
        "environment_exchange": 1,
        "shared_pool_secretion": 2,
        "crossfeeding_candidate": 3,
        "shared_pool_uptake": 4,
        "internal_processing": 5,
        "biomass_adjacent": 6,
    }
    out["role_order"] = out["role"].map(role_order).fillna(99)
    out = out.sort_values(["role_order", "abs_flux"], ascending=[True, False]).drop(columns=["role_order"])
    return out.reset_index(drop=True)


def write_summary(
    out_path: Path,
    species: str,
    diet: str,
    metabolite_id: str,
    pathway: str,
    connector_summary: pd.DataFrame,
    crossfeeding: pd.DataFrame,
    provenance_edges: pd.DataFrame,
):
    lines = [
        f"Species: {species}",
        f"Diet: {diet}",
        f"Anchor metabolite: {metabolite_id}",
        f"Selected pathway: {pathway}",
        "",
    ]

    if not connector_summary.empty:
        lines.append("Top connector fluxes for the species:")
        for row in connector_summary.head(8).itertuples(index=False):
            lines.append(
                f"- {row.shared_metabolite_id} ({row.shared_metabolite_name}): {row.direction}, flux={row.flux:.6f}"
            )
        lines.append("")

    cf_focus = crossfeeding[
        (crossfeeding["shared_metabolite_id"] == metabolite_id)
        & (crossfeeding["consumer_species"] == species)
    ].copy()
    if not cf_focus.empty:
        lines.append("Cross-feeding candidates for the anchor metabolite:")
        for row in cf_focus.sort_values("candidate_score", ascending=False).head(5).itertuples(index=False):
            lines.append(
                f"- producer={row.producer_species}, consumer={row.consumer_species}, score={row.candidate_score:.6f}"
            )
        lines.append("")

    if not provenance_edges.empty:
        lines.append("Provenance edge roles included:")
        for role in provenance_edges["role"].drop_duplicates().tolist():
            lines.append(f"- {role}")

    out_path.write_text("\n".join(lines))


def write_species_summary(
    out_path: Path,
    species: str,
    diet: str,
    support_summary: pd.DataFrame,
    pathway_overview: pd.DataFrame,
):
    lines = [
        f"Species-level biomass support summary",
        f"Species: {species}",
        f"Diet: {diet}",
        "",
    ]

    if not pathway_overview.empty:
        lines.append("Top biomass-adjacent pathways:")
        for row in pathway_overview.itertuples(index=False):
            lines.append(
                f"- {row.pathway}: score={row.biomass_association_score:.6f}, share={row.pathway_flux_share:.3f}, top reaction={row.top_reaction_id}"
            )
        lines.append("")

    if not support_summary.empty:
        lines.append("Top shared metabolite supports:")
        seen = set()
        for row in support_summary.itertuples(index=False):
            key = (row.shared_metabolite_id, row.linked_pathway)
            if key in seen:
                continue
            seen.add(key)
            lines.append(
                f"- {row.shared_metabolite_id} -> {row.linked_pathway}: uptake={row.uptake_from_shared_flux:.6f}, producer={row.top_producer_species or 'none'}, support={row.crossfeeding_support_score:.6f}"
            )
            if len(seen) >= 12:
                break

    out_path.write_text("\n".join(lines))


def write_clean_species_summary(
    out_path: Path,
    species: str,
    diet: str,
    clean_summary: pd.DataFrame,
    pathway_overview: pd.DataFrame,
):
    lines = [
        "Presentation-safe crossfeeding to biomass summary",
        f"Species: {species}",
        f"Diet: {diet}",
        "",
    ]
    if not clean_summary.empty:
        lines.append("Crossfeeding candidates retained:")
        for row in clean_summary.itertuples(index=False):
            lines.append(
                f"- {row.producer_species} -> {row.shared_metabolite_id} -> {row.consumer_species}; pathway={row.consumer_top_pathway}; support={row.crossfeeding_support_score:.3f}; confidence={row.confidence}"
            )
        lines.append("")

    if not pathway_overview.empty:
        lines.append("Dominant biomass-adjacent pathways:")
        for row in pathway_overview.head(4).itertuples(index=False):
            lines.append(
                f"- {row.pathway}: score={row.biomass_association_score:.3f}, top reaction={row.top_reaction_id}"
            )
    out_path.write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(
        description="Map a species-aware flux provenance chain from community exchange to connector fluxes to internal pathway activity."
    )
    parser.add_argument("--species", default="Alistipes_shahii_WAL_8301_AGORA1_03")
    parser.add_argument("--diet", default="high_fiber")
    parser.add_argument("--metabolite", default="adn[e]")
    parser.add_argument("--pathway", default=None)
    parser.add_argument("--top-n-connectors", type=int, default=15)
    args = parser.parse_args()

    df_ex = pd.read_csv(IN_COMM_EX)
    df_conn = pd.read_csv(IN_CONN)
    df_rxn = pd.read_csv(IN_RXN)
    df_pathway = pd.read_csv(IN_PATHWAY)
    df_biomass = pd.read_csv(IN_BIOMASS)

    connector_summary = make_connector_summary(df_conn, args.diet, args.species, args.top_n_connectors)
    crossfeeding = make_crossfeeding_candidates(df_conn, args.diet, args.species)
    pathway = pick_pathway(df_rxn, df_biomass, args.species, args.diet, args.metabolite, args.pathway)
    reaction_steps = trace_internal_pathway(df_rxn, df_biomass, args.species, args.diet, pathway)
    support_summary = build_species_biomass_support_summary(df_conn, df_biomass, df_rxn, args.species, args.diet)
    clean_summary = build_clean_crossfeeding_biomass_summary(df_conn, df_biomass, args.species, args.diet)
    pathway_overview = build_pathway_overview(df_biomass, df_rxn, args.species, args.diet)
    provenance_edges = build_provenance_edges(
        df_ex=df_ex,
        df_conn=df_conn,
        df_rxn=df_rxn,
        df_biomass=df_biomass,
        species=args.species,
        diet=args.diet,
        metabolite_id=args.metabolite,
        pathway=pathway,
    )

    base_species = slug(args.species)
    base_met = slug(args.metabolite)
    base_pathway = slug(pathway)
    connector_out = OUT_DIR / f"{slug(args.diet)}__{base_species}__connector_summary.csv"
    crossfeeding_out = OUT_DIR / f"{slug(args.diet)}__{base_species}__crossfeeding_candidates.csv"
    provenance_out = OUT_DIR / f"{slug(args.diet)}__{base_species}__{base_met}__provenance_edges.csv"
    reaction_out = OUT_DIR / f"{slug(args.diet)}__{base_species}__{base_pathway}__reaction_steps.csv"
    summary_out = OUT_DIR / f"{slug(args.diet)}__{base_species}__{base_met}__summary.txt"
    support_out = OUT_DIR / f"{slug(args.diet)}__{base_species}__biomass_support_summary.csv"
    clean_support_out = OUT_DIR / f"{slug(args.diet)}__{base_species}__clean_crossfeeding_biomass_summary.csv"
    pathway_overview_out = OUT_DIR / f"{slug(args.diet)}__{base_species}__pathway_overview.csv"
    species_summary_out = OUT_DIR / f"{slug(args.diet)}__{base_species}__biomass_support_summary.txt"
    clean_species_summary_out = OUT_DIR / f"{slug(args.diet)}__{base_species}__clean_crossfeeding_biomass_summary.txt"

    connector_summary.to_csv(connector_out, index=False)
    crossfeeding.to_csv(crossfeeding_out, index=False)
    provenance_edges.to_csv(provenance_out, index=False)
    reaction_steps.to_csv(reaction_out, index=False)
    support_summary.to_csv(support_out, index=False)
    clean_summary.to_csv(clean_support_out, index=False)
    pathway_overview.to_csv(pathway_overview_out, index=False)
    write_summary(summary_out, args.species, args.diet, args.metabolite, pathway, connector_summary, crossfeeding, provenance_edges)
    write_species_summary(species_summary_out, args.species, args.diet, support_summary, pathway_overview)
    write_clean_species_summary(clean_species_summary_out, args.species, args.diet, clean_summary, pathway_overview)

    print(f"Saved connector summary: {connector_out}")
    print(f"Saved crossfeeding candidates: {crossfeeding_out}")
    print(f"Saved provenance edges: {provenance_out}")
    print(f"Saved reaction steps: {reaction_out}")
    print(f"Saved biomass support summary: {support_out}")
    print(f"Saved clean biomass support summary: {clean_support_out}")
    print(f"Saved pathway overview: {pathway_overview_out}")
    print(f"Saved summary: {summary_out}")
    print(f"Saved species summary: {species_summary_out}")
    print(f"Saved clean species summary: {clean_species_summary_out}")
    print("")
    print("Selected pathway:")
    print(pathway)
    print("")
    print("Clean biomass support rows:")
    print(clean_summary.to_string(index=False))
    print("")
    print("Top biomass support rows:")
    print(support_summary.head(12).to_string(index=False))
    print("")
    print("Top provenance edges:")
    print(provenance_edges.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
