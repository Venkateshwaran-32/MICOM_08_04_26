# Codex Handover README

![Project cover](venkatproject.png)

## TL;DR
- Main orientation guide for understanding how the branches fit together:
  - `PIPELINE_STORY.md`
- Primary analysis: SG90-weighted community FBA (`sg90_*` files).
- Extension analysis: all-cohort age-bin FBA (`allcohort_agebin_*` files, bins `21_40` to `81_90`).
- MICOM branch now uses supervisor-approved `proper_age_bins`: `21_40`, `41_60`, `61_70`, `71_80`, `81_plus`.
- Canonical MICOM outputs are written under `Results/micom/*/proper_age_bins/` and `Results/figures/micom/growth/proper_age_bins/`.
- Older coarse-bin MICOM growth figures (`20_40`, `40_60`, `60_plus`) are retained only as legacy reference under `Results/figures/micom/growth/old_age_bin_till_60plus/`.
- Additional exploratory branches now exist for:
  - full-access community FBA (`04b`, `04c`)
  - community model structure / active reaction summaries (`06b`, `06c`)
  - species-level provenance / cross-feeding-to-biomass summaries (`06d`, `06f`, `06g`)
  - beta-age weighted community FBA (`08d`, `09c`)
- Most recent complete run branch: `08b -> 11d` plus `08c` (all-cohort age-bin inputs, scenarios, comparisons, figures).
- Next recommended step: make `Scripts/99_validate_pipeline_outputs.py` Python-3.9-compatible or run it under Python 3.10+, then use it as a final QA gate before reporting/cleanup.

## 1) Project Goal
Model gut microbial community metabolism under diet and abundance changes using:
- COBRApy weighted community FBA
- MICOM cooperative-tradeoff community modeling

Primary study focus: SG90.
Extension branch: all-cohort age-bin analysis (21-90).

If you need the project logic in plain language rather than script order, read `PIPELINE_STORY.md` first. The README remains the handover/index document; the pipeline story guide is the narrative map.

## 2) Canonical Inputs
- Models: `Models/vmh_agora_sbml/*.xml`
- Diet media: `Media/western.csv`, `Media/high_fiber.csv`
- Supplementary metadata: `Data/Supplementary/supplementary_data_file_1.xlsx`
- Supplementary taxonomy: `Data/Supplementary/taxonomic_profiles_filtered.csv`
- Model mapping manifest: `Metadata/models_manifest.csv/models_manifest.csv`

## 3) Pipeline Branches

### A. Original COBRApy Baseline (`01`-`07`)
- `01_load_qc_models.py`: SBML QC -> `Results/qc/model_qc_summary.csv`
- `02_single_species_fba_default.py`: single-species default growth -> `Results/fba/single_species_growth_default.csv`
- `03_single_species_fba_diets.py`: single-species western/high_fiber -> `Results/fba/single_species_growth_by_diet.csv`
- `04_run_community_fba.py`: unweighted community by diet -> `Results/fba/community_growth_by_diet.csv`
- `04b_run_community_fba_full_access_medium.py`: community FBA with all allowed diet exchanges set to `1000`
  - outputs: `Results/fba/community_growth_full_access.csv`, `Results/fba/community_species_biomass_flux_full_access.csv`
- `04c_plot_full_access_biomass_flux.py`: ranked biomass-flux plot for `04b`
- `05_run_fba_full_access_medium.py`: full-access species capability -> `Results/fba/individual_growth_full_access_medium.csv`
- `06_flux_diagnostics.py`: exchange/connector/reaction/pathway diagnostics
- `06b_community_model_summary.py`: structural summary of the built community model
  - output: `Results/fba/community_model_summary.csv`
- `06c_plot_active_reaction_counts.py`: plots active reaction counts by species and focused heatmap
- `06d_trace_pathway_reactions.py`: orders active internal pathway reactions by metabolite connectivity
  - outputs written under `Results/fba/pathway_traces/`
- `06f_map_species_flux_provenance.py`: builds species-aware connector, cross-feeding candidate, provenance-edge, pathway-overview, and cleaned crossfeeding-to-biomass summary tables
  - outputs written under `Results/fba/provenance/`
- `06g_plot_species_crossfeeding_biomass_summary.py`: renders provenance summary figures from the cleaned `06f` outputs
  - outputs written under `Results/figures/provenance/`
- `07_make_figures.py`: pathway + SG90 summary figures

### B. SG90 Weighted FBA (`08`-`11`)
- `08_prepare_abundance_inputs.py`
  - SG90 median + SG90 age-group inputs
- `08d_prepare_beta_age_weights.py`
  - maps paper beta coefficients to modeled species and creates `Results/inputs/beta_age_input_for_fba.csv`
- `09_run_community_fba_weighted.py`
  - scenarios: `sg90_median`, `equal_abundance`, optional `sg90_age_*`
- `09c_run_community_fba_beta_age_weighted.py`
  - scenarios: `beta_age_positive_only`, `beta_age_shifted`, `beta_age_raw_signed`
- `10_compare_weighted_vs_equal_fluxes.py`
- `11_plot_weighted_vs_equal.py`

### C. MICOM Branch (`13`-`16`)
- `13_prepare_micom_inputs.py`
- `14_run_micom_by_agegroup.py`
- `14b_rank_micom_growth_contributions.py`
- `14c_plot_micom_top_growers.py`
- `14d_tradeoff_sensitivity_micom.py`
  - sweeps cooperative-tradeoff fractions `0.0` to `1.0` and records community growth plus number of growing taxa
- `15_plot_lysine_pathways.py`
  - curated lysine biosynthesis and lysine-to-butyrate summaries/figures
  - also writes a reaction-step table for slide/report use
- `15d_rank_species_lysine_butyrate_contributions.py`
  - ranks species-level MICOM lysine and lysine-to-butyrate candidate contributions by age group and diet
  - writes canonical all-species and growers-only tables plus a ranked panel figure
- `15e_rank_species_lysine_exchange_contributions.py`
  - ranks species-level MICOM net lysine exchange with the shared environment by age group and diet
  - separates exporters from importers using `EX_lys_L(e)` rather than internal pathway sums
- `15f_rank_species_lysine_butyrate_exchange_contributions.py`
  - shows net lysine and net butyrate exchange with the shared environment in the same MICOM figure
  - uses direct exchange reactions `EX_lys_L(e)` and `EX_but(e)` with bidirectional bars around zero
  - default presentation is growers-only so the panel layout stays aligned with the `TF = 0.5` MICOM growth figure
- `16_find_differential_pathways.py`
- `16b_map_micom_flux_provenance.py`
  - builds MICOM provenance tables for selected species by age group and diet using the canonical `proper_age_bins` flux outputs
  - outputs written under `Results/micom/provenance/proper_age_bins/`
- `16c_plot_micom_provenance_summary.py`
  - renders compact MICOM provenance summary figures under `Results/figures/micom/provenance/`
- current canonical MICOM age bins: `21_40`, `41_60`, `61_70`, `71_80`, `81_plus`

### D. All-Cohort Age-Bin FBA Branch (`08b`-`11d`, `08c`)
- `08b_prepare_allcohort_agebin_inputs.py`
  - age bins: `21_40`, `41_60`, `61_70`, `71_80`, `81_90`
- `09b_run_community_fba_allcohort_agebins.py`
  - scenarios: `allcohort_agebin_*` + `allcohort_equal_abundance`
- `10b_compare_allcohort_agebins.py`
- `11b_plot_allcohort_agebins.py`
- `11c_plot_allcohort_top3_species.py`
- `11d_plot_nonzero_and_log_biomass.py`
- `08c_plot_agebin_cohort_coverage.py`

## 4) Current Status (as of latest run)

### All-cohort age-bin branch has been executed end-to-end.
Key outputs exist:
- Inputs:
  - `Data/Processed/allcohort_agebin_median_abundance_by_taxon.csv`
  - `Results/inputs/allcohort_agebin_input_for_fba.csv`
  - `Data/Processed/allcohort_agebin_input_summary.csv`
- Scenario outputs:
  - `Results/fba/scenarios/allcohort_agebin_*_community_summary_by_diet.csv`
  - `Results/fba/scenarios/allcohort_agebin_*_species_biomass_by_diet.csv`
  - `Results/fba/scenarios/allcohort_agebin_*_community_exchange_fluxes_by_diet.csv`
- Comparisons:
  - `Results/fba/allcohort_agebin_vs_equal_exchange_comparisons.csv`
  - `Results/fba/scenarios/allcohort_agebin_*_vs_allcohort_equal_exchange_flux_comparison_by_diet.csv`
- Figures:
  - `Results/figures/allcohort_agebin_vs_equal_flux_delta_top30.png`
  - `Results/figures/allcohort_agebin_vs_equal_flux_delta_heatmap_top40.png`
  - `Results/figures/allcohort_agebin_top3_species_by_diet.png`
  - `Results/figures/allcohort_agebin_nonzero_species_count_by_diet.png`
  - `Results/figures/allcohort_agebin_species_biomass_log10_western.png`
  - `Results/figures/allcohort_agebin_species_biomass_log10_high_fiber.png`
  - `Results/figures/allcohort_agebin_cohort_sample_counts_used.png`

### Additional exploratory outputs now exist / are supported
- Full-access community branch:
  - `Results/fba/community_growth_full_access.csv`
  - `Results/fba/community_species_biomass_flux_full_access.csv`
  - `Results/figures/community_species_biomass_flux_full_access_ranked.png`
- Community diagnostics summary branch:
  - `Results/fba/community_model_summary.csv`
  - `Results/figures/active_reaction_counts_by_species_diet.png`
  - `Results/figures/active_reaction_counts_main_species_heatmap.png`
- COBRApy provenance / cross-feeding summary branch:
  - `Results/fba/provenance/high_fiber__Alistipes_shahii_WAL_8301_AGORA1_03__clean_crossfeeding_biomass_summary.csv`
  - `Results/fba/provenance/high_fiber__Alistipes_shahii_WAL_8301_AGORA1_03__pathway_overview.csv`
  - `Results/fba/provenance/high_fiber__Alistipes_shahii_WAL_8301_AGORA1_03__adn_e__provenance_edges.csv`
  - `Results/fba/provenance/western__Faecalibacterium_prausnitzii_M21_2__AGORA1_03__clean_crossfeeding_biomass_summary.csv`
  - `Results/fba/provenance/western__Faecalibacterium_prausnitzii_M21_2__AGORA1_03__pathway_overview.csv`
  - `Results/fba/provenance/western__Faecalibacterium_prausnitzii_M21_2__AGORA1_03__dcyt_e__provenance_edges.csv`
  - `Results/figures/provenance/alistipes_shahii_crossfeeding_biomass_summary.png`
  - `Results/figures/provenance/faecalibacterium_prausnitzii_crossfeeding_biomass_summary.png`
- Beta-age weighting branch:
  - `Data/Processed/beta_age_taxon_mapping.csv`
  - `Results/inputs/beta_age_input_for_fba.csv`
  - scenario outputs written to `Results/fba/scenarios/beta_age_*`
- MICOM branch:
  - `Results/micom/growth/proper_age_bins/community_growth_summary_by_agegroup_diet.csv`
  - `Results/micom/growth/proper_age_bins/organism_growth_rates_by_agegroup_diet.csv`
  - `Results/micom/growth/proper_age_bins/micom_growth_contribution_ranking_by_agegroup_diet.csv`
  - `Results/micom/pathway_flux/proper_age_bins/reaction_fluxes_long_by_agegroup_diet.csv`
  - `Results/micom/pathway_flux/proper_age_bins/reaction_annotations_by_organism.csv`
  - `Results/micom/tradeoff_sensitivity/proper_age_bins/tradeoff_sensitivity_summary.csv`
  - `Results/micom/tradeoff_sensitivity/proper_age_bins/tradeoff_sensitivity_species_growth.csv`
  - `Results/micom/provenance/proper_age_bins/micom_provenance_case_overview.csv`
  - case-level provenance tables for `Faecalibacterium_prausnitzii_M21_2__AGORA1_03` and `Escherichia_coli_UTI89_UPEC_AGORA1_03` in `61_70`, `71_80`, and `81_plus` at tradeoff fraction `0.5`
  - `Results/micom/growth/proper_age_bins/scaled_uptake_community_growth_by_agegroup_diet.csv`
  - `Results/micom/growth/proper_age_bins/scaled_uptake_species_growth_by_agegroup_diet.csv`
  - `Results/micom/pathway_flux/pathway_flux_by_agegroup_diet.csv`
  - `Results/micom/pathway_flux/differential_pathways_by_agegroup_diet.csv`
  - `Results/micom/lysine_butyrate/lysine_related_fluxes_long.csv`
  - `Results/micom/lysine_butyrate/lysine_biosynthesis_flux_summary.csv`
  - `Results/micom/lysine_butyrate/lysine_to_butyrate_flux_summary.csv`
  - `Results/micom/lysine_butyrate/lysine_to_butyrate_reaction_step_summary.csv`
  - `Results/micom/lysine_butyrate/species_lysine_butyrate_contributions_by_agegroup_diet.csv`
  - `Results/micom/lysine_butyrate/species_lysine_butyrate_contributions_ranked_growers_only.csv`
  - `Results/micom/lysine_butyrate/species_lysine_exchange_contributions_by_agegroup_diet.csv`
  - `Results/micom/lysine_butyrate/species_lysine_exchange_contributions_ranked_growers_only.csv`
  - `Results/micom/lysine_butyrate/species_lysine_butyrate_exchange_contributions_by_agegroup_diet.csv`
  - `Results/micom/lysine_butyrate/species_lysine_butyrate_exchange_contributions_ranked_growers_only.csv`
  - `Results/figures/micom/pathway_flux/differential_pathways_heatmap_top25.png`
  - `Results/figures/micom/growth/proper_age_bins/micom_all_growing_contributors_by_agegroup_diet_proper_age_bins.png`
  - `Results/figures/micom/growth/proper_age_bins/micom_all_growing_species_growth_rates_by_agegroup_diet_proper_age_bins.png`
  - `Results/figures/micom/growth/proper_age_bins/micom_input_abundances_by_agegroup_proper_age_bins.png`
  - `Results/figures/micom/growth/proper_age_bins/micom_tradeoff_species_growth_curves_proper_age_bins.png`
  - `Results/figures/micom/growth/proper_age_bins/micom_tradeoff_community_growth_curves_proper_age_bins.png`
  - `Results/figures/micom/lysine_butyrate/lysine_biosynthesis_flux_by_agegroup_diet.png`
  - `Results/figures/micom/lysine_butyrate/lysine_to_butyrate_flux_by_agegroup_diet.png`
  - `Results/figures/micom/lysine_butyrate/lysine_to_butyrate_reaction_steps_by_agegroup_diet.png`
  - `Results/figures/micom/lysine_butyrate/species_lysine_butyrate_contributions_by_agegroup_diet.png`
  - `Results/figures/micom/lysine_butyrate/species_lysine_exchange_contributions_by_agegroup_diet.png`
  - `Results/figures/micom/lysine_butyrate/species_lysine_butyrate_exchange_contributions_by_agegroup_diet.png`
  - `Results/figures/micom/provenance/Escherichia_coli_UTI89_UPEC_AGORA1_03__micom_provenance_summary.png`
  - `Results/figures/micom/provenance/Faecalibacterium_prausnitzii_M21_2__AGORA1_03__micom_provenance_summary.png`
  - legacy coarse-bin growth figures are archived under `Results/figures/micom/growth/old_age_bin_till_60plus/`

## 5) Important Interpretation Notes
- Weighted LP community FBA can produce sparse optima (few winners, many zero-growth taxa).
- Zero growth in weighted community output != model broken.
  - Always cross-check with full-access capability outputs.
- COBRApy and MICOM are parallel modeling branches, not a single linear replacement.
  - Both have age-bin analyses, but they answer them under different optimization frameworks.
- Lysine/butanoate and provenance are interpretation layers on top of primary modeling outputs.
  - They should not be treated as separate primary pipelines.
- MICOM branch is used to reduce winner-take-all behavior via cooperative tradeoff.
- The MICOM tradeoff fraction directly affects both community growth and how many taxa receive nonzero growth.
  - `14d_tradeoff_sensitivity_micom.py` was added to test this explicitly across multiple fractions.
- For current reporting, treat the `proper_age_bins` MICOM branch as canonical.
  - Older coarse-bin MICOM growth figures use legacy bins `20-39`, `40-59`, and `60+` and should not be mixed with `proper_age_bins`.
- The MICOM provenance pilot is currently scoped to `61_70`, `71_80`, and `81_plus` at tradeoff fraction `0.5`.
  - Interpret these outputs as inferred exchange-to-pathway provenance, not the same connector-explicit provenance used in the COBRApy branch.
- The regenerated lysine/butyrate MICOM summaries now also use `proper_age_bins` as their source of truth.
  - Do not interpret older lysine/butyrate tables that still show `20_40`, `40_60`, or `60_plus`.
- Species-level `sum_abs_flux` lysine/butyrate tables reflect pathway activity magnitude, not direct export to the shared environment.
  - For shared-environment interpretation, use the lysine exchange branch built from `EX_lys_L(e)`.
- For direct extracellular interpretation of both metabolites together, use the combined exchange branch built from `EX_lys_L(e)` and `EX_but(e)`.
- `Scripts/15_plot_lysine_pathways.py` now uses curated reaction IDs/names for lysine-to-butyrate candidate steps instead of only broad text fragments.
- In the current MICOM runs, curated lysine-to-butyrate step summaries retain all age groups, but only downstream butanoate steps show nonzero flux.
  - Upstream lysine-specific steps currently remain at zero in these solutions.
- In this project, western and high_fiber media use the same exchange IDs and differ mainly in uptake magnitudes.
  - Therefore, `04b` full-access runs largely remove the original diet contrast.
- Community model structure counts (reaction/metabolite totals) are diet-independent.
  - Diet affects the optimized flux solution, not the underlying model size.
- SG90 and all-cohort outputs are intentionally separated by filename prefix:
  - SG90: `sg90_*`
  - all-cohort: `allcohort_agebin_*`

## 6) Known Pitfalls / Gotchas
- `community_species_biomass_flux_by_diet.csv` was previously overwritten by a medium table and removed.
  - Do not rely on old copies of that file.
- Some plots previously had readability issues; latest scripts prioritize cleaner labels/layout.
- In all-cohort bins, some bins have only 1-2 species with positive growth under weighted LP objective.
- Beta-age coefficients include negative values.
  - These should not be used as raw LP weights without a transformation step.
- Broad text searches for `but*` reactions can overcall butyrate-related chemistry.
  - Prefer curated reaction IDs/names for lysine-to-butyrate inspection; pathway-level plots alone are too coarse for mechanistic claims.
- `Scripts/99_validate_pipeline_outputs.py` uses Python 3.10+ type-union syntax.
  - In the current `.venv` (`Python 3.9.6`), it fails before running checks unless the environment is upgraded or the annotations are backported.

## 7) Suggested Next Codex Actions
1. Make `Scripts/99_validate_pipeline_outputs.py` runnable in the project environment, then use it as the final QA gate.
2. Keep presentation-first plotting defaults (short labels, no overlapping text).
3. If reporting age trends, use all-cohort branch for bins and SG90 branch as focused validation.
4. Archive non-essential scenario files only after validation (move, do not delete).
5. Start MICOM branch walkthrough from `Scripts/13_prepare_micom_inputs.py`, using `proper_age_bins` outputs as the source of truth.

## 8) Quick Run Commands (macOS / zsh)
From project root:
- Activate the environment:
  - `source .venv/bin/activate`
- If `cobra` or `matplotlib` try to write cache/config files to blocked locations, run scripts with:
  - `HOME="$PWD" python Scripts/<script_name>.py`

- All-cohort branch:
  - `HOME="$PWD" python Scripts/08b_prepare_allcohort_agebin_inputs.py`
  - `HOME="$PWD" python Scripts/09b_run_community_fba_allcohort_agebins.py`
  - `HOME="$PWD" python Scripts/10b_compare_allcohort_agebins.py`
  - `HOME="$PWD" python Scripts/11b_plot_allcohort_agebins.py`
  - `HOME="$PWD" python Scripts/11c_plot_allcohort_top3_species.py`
  - `HOME="$PWD" python Scripts/11d_plot_nonzero_and_log_biomass.py`
  - `HOME="$PWD" python Scripts/08c_plot_agebin_cohort_coverage.py`

- MICOM branch:
  - `HOME="$PWD" python Scripts/13_prepare_micom_inputs.py`
  - `HOME="$PWD" python Scripts/14_run_micom_by_agegroup.py`
  - `HOME="$PWD" python Scripts/14b_rank_micom_growth_contributions.py`
  - `HOME="$PWD" python Scripts/14c_plot_micom_top_growers.py`
  - `HOME="$PWD" python Scripts/14d_tradeoff_sensitivity_micom.py`
  - `HOME="$PWD" python Scripts/15_plot_lysine_pathways.py`
  - `HOME="$PWD" python Scripts/15d_rank_species_lysine_butyrate_contributions.py`
  - `HOME="$PWD" python Scripts/15e_rank_species_lysine_exchange_contributions.py`
  - `HOME="$PWD" python Scripts/15f_rank_species_lysine_butyrate_exchange_contributions.py`
  - `HOME="$PWD" python Scripts/16_find_differential_pathways.py`
  - `HOME="$PWD" python Scripts/16b_map_micom_flux_provenance.py`
  - `HOME="$PWD" python Scripts/16c_plot_micom_provenance_summary.py`

- Full-access community branch:
  - `HOME="$PWD" python Scripts/04b_run_community_fba_full_access_medium.py`
  - `HOME="$PWD" python Scripts/04c_plot_full_access_biomass_flux.py`

- Community diagnostics summary branch:
  - `HOME="$PWD" python Scripts/06_flux_diagnostics.py`
  - `HOME="$PWD" python Scripts/06b_community_model_summary.py`
  - `HOME="$PWD" python Scripts/06c_plot_active_reaction_counts.py`

- Beta-age weighting branch:
  - `HOME="$PWD" python Scripts/08d_prepare_beta_age_weights.py`
  - `HOME="$PWD" python Scripts/09c_run_community_fba_beta_age_weighted.py`

## 9) Environment
- Python env: `.venv`
- `.venv` currently reports `Python 3.9.6`
- Core deps: `cobra`, `pandas`, `matplotlib`
- MICOM deps: `micom`
