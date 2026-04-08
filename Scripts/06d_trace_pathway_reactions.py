from pathlib import Path
import argparse
import re
from collections import Counter, defaultdict

import cobra
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODELS_DIR = PROJECT_ROOT / "Models" / "vmh_agora_sbml"
FBA_DIR = PROJECT_ROOT / "Results" / "fba"
IN_RXN_CSV = FBA_DIR / "community_reaction_fluxes_by_diet.csv"
OUT_DIR = FBA_DIR / "pathway_traces"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FLUX_EPS = 1e-9
COMMON_COFATOR_METS = {
    "h[c]", "h[e]", "h2o[c]", "h2o[e]", "pi[c]", "pi[e]", "ppi[c]", "ppi[e]",
    "atp[c]", "atp[e]", "adp[c]", "adp[e]", "amp[c]", "amp[e]",
    "nad[c]", "nad[e]", "nadh[c]", "nadh[e]", "nadp[c]", "nadp[e]", "nadph[c]", "nadph[e]",
    "coa[c]", "coa[e]", "co2[c]", "co2[e]", "o2[c]", "o2[e]",
}


def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", text)


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


def find_intermediates(records):
    produced_by = defaultdict(set)
    consumed_by = defaultdict(set)
    occurrence_counter = Counter()

    for rec in records:
        rid = rec["reaction_id"]
        for met_id, _ in rec["substrates"]:
            consumed_by[met_id].add(rid)
            occurrence_counter[met_id] += 1
        for met_id, _ in rec["products"]:
            produced_by[met_id].add(rid)
            occurrence_counter[met_id] += 1

    rows = []
    for met_id in sorted(set(produced_by) | set(consumed_by)):
        producers = sorted(produced_by.get(met_id, set()))
        consumers = sorted(consumed_by.get(met_id, set()))
        if met_id in COMMON_COFATOR_METS:
            continue
        if not producers or not consumers:
            continue
        rows.append(
            {
                "metabolite_id": met_id,
                "n_occurrences": occurrence_counter[met_id],
                "produced_by": " | ".join(producers),
                "consumed_by": " | ".join(consumers),
            }
        )
    return pd.DataFrame(rows)


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
    out_degree = Counter()
    adjacency = defaultdict(list)
    reverse_adjacency = defaultdict(list)

    for edge in edge_list:
        src = edge["from_reaction_id"]
        dst = edge["to_reaction_id"]
        adjacency[src].append(edge)
        reverse_adjacency[dst].append(edge)
        out_degree[src] += 1
        in_degree[dst] += 1

    # Prefer a starting point that has no incoming pathway edges; otherwise use the
    # reaction with the largest abs flux as the anchor.
    candidate_starts = [
        rec["reaction_id"]
        for rec in sorted(records, key=lambda r: (-r["abs_flux"], r["reaction_id"]))
        if in_degree[rec["reaction_id"]] == 0
    ]
    if candidate_starts:
        current = candidate_starts[0]
    else:
        current = sorted(records, key=lambda r: (-r["abs_flux"], r["reaction_id"]))[0]["reaction_id"]

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

        outgoing = [
            e for e in adjacency.get(current, [])
            if e["to_reaction_id"] not in visited
        ]
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


def write_pipeline_text(path: Path, diet: str, species: str, pathway: str, pipeline_table: pd.DataFrame):
    lines = [
        f"Diet: {diet}",
        f"Species: {species}",
        f"Pathway: {pathway}",
        "",
        "Stepwise reaction pipeline",
        "This ordering is based on metabolite connectivity within the active pathway reactions.",
        "It is a practical trace, not a claim of unique biochemical order.",
        "",
    ]
    for row in pipeline_table.itertuples(index=False):
        lines.append(f"Step {row.step}: {row.reaction_id} ({row.reaction_name})")
        lines.append(f"  Flux: {row.flux:.6f}")
        if row.incoming_from:
            lines.append(f"  Connected from: {row.incoming_from}")
            lines.append(f"  Shared intermediate(s): {row.shared_metabolites_from_previous}")
        lines.append(f"  Direction used in solution: {row.equation_in_flux_direction}")
        lines.append("")
    path.write_text("\n".join(lines))


def trace_pathway(diet: str, species: str, pathway: str, flux_eps: float):
    df = pd.read_csv(IN_RXN_CSV)
    use = df[
        (df["diet"].astype(str) == diet)
        & (df["species"].astype(str) == species)
        & (df["pathway"].astype(str) == pathway)
        & (df["reaction_type"].astype(str) == "species_internal")
        & (df["abs_flux"].astype(float) > flux_eps)
    ].copy()

    if use.empty:
        raise ValueError("No active species_internal reactions matched the requested diet/species/pathway.")

    model = load_species_model(species)
    records = []
    for row in use.sort_values(["abs_flux", "local_reaction_id"], ascending=[False, True]).itertuples(index=False):
        local_id = row.local_reaction_id
        if local_id not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(local_id)
        substrates, products = directed_sides(rxn, float(row.flux))
        records.append(
            {
                "diet": diet,
                "species": species,
                "pathway": pathway,
                "reaction_id": local_id,
                "reaction_name": rxn.name,
                "flux": float(row.flux),
                "abs_flux": float(row.abs_flux),
                "equation_in_flux_direction": f"{side_to_text(substrates)} -> {side_to_text(products)}",
                "substrates": substrates,
                "products": products,
            }
        )

    if not records:
        raise ValueError("Matched rows in diagnostics CSV but could not map local reaction IDs back to the SBML model.")

    reaction_table = pd.DataFrame(
        [
            {
                "diet": rec["diet"],
                "species": rec["species"],
                "pathway": rec["pathway"],
                "reaction_id": rec["reaction_id"],
                "reaction_name": rec["reaction_name"],
                "flux": rec["flux"],
                "abs_flux": rec["abs_flux"],
                "equation_in_flux_direction": rec["equation_in_flux_direction"],
            }
            for rec in records
        ]
    )
    intermediate_table = find_intermediates(records)
    edge_table = build_edges(records)
    pipeline_table = build_stepwise_pipeline(records)
    return reaction_table, intermediate_table, edge_table, pipeline_table


def slug(text: str) -> str:
    return sanitize(text).strip("_")


def main():
    parser = argparse.ArgumentParser(
        description="Trace active reactions and intermediate metabolites for one diet/species/pathway from community FBA diagnostics."
    )
    parser.add_argument("--diet", required=True, help="Diet name exactly as it appears in community_reaction_fluxes_by_diet.csv")
    parser.add_argument("--species", required=True, help="Sanitized species id from the diagnostics CSV")
    parser.add_argument("--pathway", required=True, help="Pathway name from the diagnostics CSV")
    parser.add_argument("--flux-eps", type=float, default=FLUX_EPS, help="Minimum abs flux to keep a reaction active")
    args = parser.parse_args()

    reaction_table, intermediate_table, edge_table, pipeline_table = trace_pathway(
        diet=args.diet,
        species=args.species,
        pathway=args.pathway,
        flux_eps=args.flux_eps,
    )

    base = f"{slug(args.diet)}__{slug(args.species)}__{slug(args.pathway)}"
    rxn_out = OUT_DIR / f"{base}_reactions.csv"
    met_out = OUT_DIR / f"{base}_intermediates.csv"
    edge_out = OUT_DIR / f"{base}_edges.csv"
    pipeline_out = OUT_DIR / f"{base}_pipeline.csv"
    pipeline_txt_out = OUT_DIR / f"{base}_pipeline.txt"

    reaction_table.to_csv(rxn_out, index=False)
    intermediate_table.to_csv(met_out, index=False)
    edge_table.to_csv(edge_out, index=False)
    pipeline_table.to_csv(pipeline_out, index=False)
    write_pipeline_text(pipeline_txt_out, args.diet, args.species, args.pathway, pipeline_table)

    print(f"Saved reactions: {rxn_out}")
    print(f"Saved intermediates: {met_out}")
    print(f"Saved edges: {edge_out}")
    print(f"Saved pipeline: {pipeline_out}")
    print(f"Saved pipeline text: {pipeline_txt_out}")
    print("")
    print("Stepwise pipeline:")
    print(pipeline_table.head(15).to_string(index=False))
    print("")
    print("Intermediate metabolites:")
    if intermediate_table.empty:
        print("None found.")
    else:
        print(intermediate_table.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
