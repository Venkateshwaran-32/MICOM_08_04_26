from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_MICOM = PROJECT_ROOT / "Results" / "micom"
GROWTH_DIR = RESULTS_MICOM / "growth"
PATHWAY_DIR = RESULTS_MICOM / "pathway_flux"
LYS_DIR = RESULTS_MICOM / "lysine_butyrate"
SUMMARY_DIR = RESULTS_MICOM / "summary_stats"
SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

IN_ANNOT = PATHWAY_DIR / "reaction_annotations_by_organism.csv"
IN_FLUX = PATHWAY_DIR / "reaction_fluxes_long_by_agegroup_diet.csv"
IN_PATHWAY = PATHWAY_DIR / "pathway_flux_by_agegroup_diet.csv"
IN_LYS = LYS_DIR / "lysine_related_fluxes_long.csv"

OUT_CSV = SUMMARY_DIR / "micom_pathway_space_summary.csv"


def build_rows() -> list[dict]:
    rows = []

    annot = pd.read_csv(IN_ANNOT)
    rows.extend(
        [
            {
                "section": "overall_micom_reaction_pathway_space",
                "metric": "n_organisms",
                "value": int(annot["id"].nunique()),
                "source_file": IN_ANNOT.name,
                "notes": "Unique organisms with reaction annotations.",
            },
            {
                "section": "overall_micom_reaction_pathway_space",
                "metric": "n_organism_specific_reaction_rows",
                "value": int(len(annot)),
                "source_file": IN_ANNOT.name,
                "notes": "One row per organism-reaction annotation pair.",
            },
            {
                "section": "overall_micom_reaction_pathway_space",
                "metric": "n_unique_reaction_ids_ignore_organism",
                "value": int(annot["reaction_id"].nunique()),
                "source_file": IN_ANNOT.name,
                "notes": "Unique reaction IDs after collapsing across organisms.",
            },
            {
                "section": "overall_micom_reaction_pathway_space",
                "metric": "n_annotated_pathways",
                "value": int(annot["pathway"].fillna("NA").nunique()),
                "source_file": IN_ANNOT.name,
                "notes": "Unique pathway labels in reaction annotations.",
            },
        ]
    )

    flux = pd.read_csv(IN_FLUX)
    rows.extend(
        [
            {
                "section": "overall_micom_reaction_pathway_space",
                "metric": "n_flux_rows",
                "value": int(len(flux)),
                "source_file": IN_FLUX.name,
                "notes": "Long-format organism-reaction flux rows across all age groups and diets.",
            },
            {
                "section": "overall_micom_reaction_pathway_space",
                "metric": "n_agegroup_diet_blocks",
                "value": int(flux[["age_group", "diet"]].drop_duplicates().shape[0]),
                "source_file": IN_FLUX.name,
                "notes": "Unique age-group by diet runs.",
            },
        ]
    )

    pathway = pd.read_csv(IN_PATHWAY)
    rows.append(
        {
            "section": "overall_micom_reaction_pathway_space",
            "metric": "n_pathways_with_nonzero_summarized_flux",
            "value": int(pathway["pathway"].nunique()),
            "source_file": IN_PATHWAY.name,
            "notes": "Unique pathways with summarized nonzero flux signal across runs.",
        }
    )

    lys = pd.read_csv(IN_LYS)
    rows.extend(
        [
            {
                "section": "lysine_keyword_filtered_subset",
                "metric": "n_rows",
                "value": int(len(lys)),
                "source_file": IN_LYS.name,
                "notes": "Rows retained by lysine-related keyword filter.",
            },
            {
                "section": "lysine_keyword_filtered_subset",
                "metric": "n_unique_reaction_ids",
                "value": int(lys["reaction_id"].nunique()),
                "source_file": IN_LYS.name,
                "notes": "Unique reaction IDs in the lysine-focused subset.",
            },
            {
                "section": "lysine_keyword_filtered_subset",
                "metric": "n_unique_organism_reaction_pairs",
                "value": int(lys[["id", "reaction_id"]].drop_duplicates().shape[0]),
                "source_file": IN_LYS.name,
                "notes": "Unique organism-reaction combinations in the lysine-focused subset.",
            },
            {
                "section": "lysine_keyword_filtered_subset",
                "metric": "n_pathways_touched",
                "value": int(lys["pathway"].fillna("NA").nunique()),
                "source_file": IN_LYS.name,
                "notes": "Unique pathway labels represented in the lysine-focused subset.",
            },
        ]
    )

    return rows


def main():
    for p in [IN_ANNOT, IN_FLUX, IN_PATHWAY, IN_LYS]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required input: {p}")

    out = pd.DataFrame(build_rows())
    out.to_csv(OUT_CSV, index=False)
    print(f"Saved: {OUT_CSV}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
