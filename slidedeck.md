# Slide Deck Outline: Host-Microbiome Interaction Modelling Project

## Slide 1: Project Title, Aim, and Main Question
**Purpose**
Open with the big-picture aim of the project.

**Key points**
- Project aim: model how gut microbial community metabolism changes under diet and abundance differences.
- Main frameworks used:
  - COBRApy-based FBA
  - MICOM cooperative-tradeoff community modelling
- Core question:
  - how diet and abundance structure alter community growth, exchange behaviour, and pathway activity

**Suggested visual**
- [INSERT TITLE VISUAL OR SIMPLE PROJECT ARC]
- Source: no file required; create manually in slides
- If unavailable: use a plain text subtitle `baseline FBA -> weighted community FBA -> MICOM -> pathway interpretation`

---

## Slide 2: Why This Project Matters
**Purpose**
Frame the biological and modelling motivation.

**Key points**
- Gut microbes interact through shared metabolites, competition, and cross-feeding.
- Diet changes nutrient availability and therefore the feasible metabolic solution space.
- Relative abundance matters because a community objective should not always treat all taxa equally.
- Pathway-level interpretation is needed because growth summaries alone do not explain mechanism.

**Suggested visual**
- [INSERT CONCEPT DIAGRAM: diet -> abundance -> community metabolism -> pathway outputs]
- Source: no direct file; create manually
- If unavailable: use a 4-box process diagram made in PowerPoint/NotebookLM

**Speaker notes**
- Emphasize that each new branch was added to answer a limitation in the previous one.

---

## Slide 3: Input Data and Model Resources
**Purpose**
Show the upstream resources feeding the whole pipeline.

**Key points**
- Genome-scale metabolic models:
  - `Models/vmh_agora_sbml/*.xml`
- Diet media:
  - `Media/western.csv`
  - `Media/high_fiber.csv`
- Supplementary metadata and taxonomy:
  - age information
  - taxonomic abundance profiles
- Mapping manifest:
  - links taxonomic labels to local model files

**Suggested table**
- [INSERT TABLE: input resources]
- Source files:
  - `Models/vmh_agora_sbml/*.xml`
  - `Media/western.csv`
  - `Media/high_fiber.csv`
  - `Data/Supplementary/supplementary_data_file_1.xlsx`
  - `Data/Supplementary/taxonomic_profiles_filtered.csv`
  - `Metadata/models_manifest.csv/models_manifest.csv`
- If unavailable: manually create a 3-column table with `resource`, `role`, `used in branch`

**Definitions to include**
- `medium`: allowed environmental uptake bounds
- `exchange reaction`: uptake/secretion interface between model and environment

---

## Slide 4: End-to-End Pipeline Map
**Purpose**
Orient the audience before going branch by branch.

**Key points**
- Baseline COBRApy branch:
  - `01` to `07`
- Weighted community FBA branches:
  - SG90
  - all-cohort age-bin
  - beta-age
- MICOM branch:
  - `13` to `16`
- Interpretation and diagnostic branches:
  - full-access checks
  - flux diagnostics
  - pathway tracing
  - provenance and cross-feeding summaries

**Suggested visual**
- [INSERT PIPELINE MAP]
- Source references:
  - `README.md`
  - `Scripts/01_load_qc_models.py`
  - `Scripts/09_run_community_fba_weighted.py`
  - `Scripts/14_run_micom_by_agegroup.py`
- If unavailable: manually create a horizontal flow chart with branch labels

**Speaker notes**
- Treat SG90 as the primary weighted branch and `proper_age_bins` as the canonical MICOM branch.

---

## Slide 5: Baseline COBRApy, Part 1 — QC and Single-Species FBA
**Purpose**
Explain the first layer of the COBRApy workflow.

**Key points**
- `01_load_qc_models.py` checks that SBML models load correctly and records model size and objective information.
- `02_single_species_fba_default.py` tests default single-species growth.
- `03_single_species_fba_diets.py` compares growth under `western` and `high_fiber` media.
- These steps establish which models are usable and whether diet constraints change single-species feasibility.

**Suggested table**
- [INSERT TABLE: QC and single-species outputs]
- Source files:
  - `Results/qc/model_qc_summary.csv`
  - `Results/fba/single_species_growth_default.csv`
  - `Results/fba/single_species_growth_by_diet.csv`
- Recommended columns:
  - QC: `file`, `model_id`, `reactions`, `metabolites`, `genes`
  - Growth: species/model id, diet, growth value
- If unavailable: summarize in 3 bullets instead of a table

**Likely supervisor question**
- Why start with single-species FBA if the project is about communities?

**Answer to defend**
- It is the baseline capability check. Community results are easier to interpret if single-species feasibility and diet response are already known.

---

## Slide 6: Baseline COBRApy, Part 2 — Community Model Construction
**Purpose**
Show how multiple species models are assembled into one community.

**Key points**
- `04_run_community_fba.py` creates a community model by:
  - prefixing species-specific reactions and metabolites
  - linking extracellular metabolites to a shared pool
  - adding community-level exchange reactions
  - optimizing the community biomass objective
- This lets species remain internally distinct while interacting through shared extracellular metabolites.

**Suggested visual**
- [INSERT COMMUNITY ASSEMBLY DIAGRAM]
- Source references:
  - `Scripts/04_run_community_fba.py`
  - `Results/fba/community_model_summary.csv`
- If unavailable: draw manually with 2 species boxes, one shared pool, and exchange arrows

**Definitions to include**
- `biomass reaction`: growth proxy for one species
- `community objective`: optimization target over the assembled multi-species model
- `connector reaction`: reversible link between species extracellular space and shared pool

**Likely supervisor question**
- Why are connector reactions necessary?

**Answer to defend**
- They enable cross-feeding through a shared metabolite pool without collapsing all species into one undifferentiated network.

---

## Slide 7: Baseline COBRApy, Part 3 — What This Baseline Establishes
**Purpose**
Summarize what the baseline community branch is useful for.

**Key points**
- Baseline outputs show:
  - community growth by diet
  - initial species biomass distribution
  - whether the community model solves under the given media
- Full-access diagnostics (`04b`, `04c`, `05`) test capability under relaxed uptake bounds.
- These runs help distinguish:
  - true incapability
  - diet limitation
  - competition-driven suppression

**Suggested visual**
- [INSERT FIGURE: baseline or full-access branch]
- Source files:
  - `Results/fba/community_growth_by_diet.csv`
  - `Results/fba/community_growth_full_access.csv`
  - `Results/figures/community_species_biomass_flux_full_access_ranked.png`
- If unavailable: create a simple 2-row comparison table: `diet-constrained` vs `full-access`

**Speaker notes**
- Make clear that full-access runs are diagnostic controls, not the main biological answer.

---

## Slide 8: Why Weighted Community FBA Was Added
**Purpose**
Explain why the project moved beyond unweighted community FBA.

**Key points**
- Equal-abundance community FBA ignores observed abundance structure.
- A weighted biomass objective lets abundant taxa influence the solution more strongly.
- This branch was added to ask whether abundance structure changes:
  - which species dominate growth
  - which community exchange fluxes appear
  - which pathway summaries become prominent

**Suggested visual**
- [INSERT OBJECTIVE COMPARISON DIAGRAM]
- Source references:
  - `Scripts/08_prepare_abundance_inputs.py`
  - `Scripts/09_run_community_fba_weighted.py`
- If unavailable: manually show two formulas or two bullet blocks: `equal-abundance` vs `weighted objective`

**Definitions to include**
- `weighted objective`: species biomass reactions have different coefficients
- `sparse optimum`: solution where only a few taxa carry most growth and many remain at zero

**Likely supervisor question**
- Does zero growth imply a broken species model?

**Answer to defend**
- No. Weighted LP community FBA can produce sparse optima. Full-access checks are needed before calling it model failure.

---

## Slide 9: SG90 Weighted Branch — Inputs and Scenario Logic
**Purpose**
Introduce the main weighted-FBA branch.

**Key points**
- `08_prepare_abundance_inputs.py` generates SG90-weighted inputs.
- `09_run_community_fba_weighted.py` applies those weights to biomass reactions.
- Main scenarios:
  - `sg90_median`
  - `equal_abundance`
  - optional `sg90_age_*`
- `10` and `11` compare weighted vs equal-abundance community flux behaviour.

**Suggested table**
- [INSERT TABLE: SG90 scenarios]
- Source files:
  - `Results/inputs/sg90_median_input_for_fba.csv`
  - `Results/inputs/agegroup_median_input_for_fba.csv`
  - `Results/fba/scenarios/sg90_median_community_summary_by_diet.csv`
  - `Results/fba/scenarios/equal_abundance_community_summary_by_diet.csv`
- Recommended columns:
  - `scenario`, `diet`, `community_objective`, `applied_exchanges`
- If unavailable: manually create a compact scenario logic table

**Speaker notes**
- Present SG90 as the primary weighted branch for the project story.

---

## Slide 10: SG90 Weighted Branch — Main Results
**Purpose**
Show what the SG90 weighted branch actually contributes.

**Key points**
- Main outputs:
  - community summary by diet
  - species biomass by diet
  - exchange fluxes by diet
  - weighted vs equal-abundance comparisons
- These outputs help identify dominant growers under SG90-informed weighting.
- They also show how weighting changes exchange behaviour relative to equal-abundance assumptions.

**Suggested visual**
- [INSERT FIGURE: SG90 weighted vs equal]
- Source files:
  - `Results/figures/sg90_vs_equal_flux_delta_top25.png`
  - `Results/fba/scenarios/sg90_median_species_biomass_by_diet.csv`
  - `Results/fba/scenarios/sg90_median_vs_equal_exchange_flux_comparison_by_diet.csv`
- If unavailable: use a manually created summary table of top species and top delta exchanges

**Interpretation caution**
- These outputs support pattern recognition and comparative interpretation, not direct pathway mechanism claims.

---

## Slide 11: All-Cohort Age-Bin Branch — Why It Was Added
**Purpose**
Explain why the project expanded beyond SG90.

**Key points**
- `08b` to `11d` extend the weighted-FBA logic across multiple age bins.
- This branch asks whether age-group abundance structure changes community metabolic behaviour.
- It broadens the project from a focused SG90 branch to a wider age-bin comparison.
- Current age bins in this branch:
  - `21_40`, `41_60`, `61_70`, `71_80`, `81_90`

**Suggested visual**
- [INSERT AGE-BIN WORKFLOW PANEL]
- Source references:
  - `Scripts/08b_prepare_allcohort_agebin_inputs.py`
  - `Scripts/09b_run_community_fba_allcohort_agebins.py`
  - `Scripts/10b_compare_allcohort_agebins.py`
- If unavailable: use a simple age-bin branch flow built manually

**Speaker notes**
- Position this as the most recent end-to-end COBRApy branch.

---

## Slide 12: All-Cohort Age-Bin Branch — Main Results
**Purpose**
Show the key outputs from the all-cohort extension.

**Key points**
- The branch has been executed end-to-end with scenario, comparison, and figure outputs.
- Example community objective values differ across age bins and diets.
- In the current outputs, some age bins show only a small number of taxa with positive growth under the weighted LP objective.
- This branch is useful for age-trend comparison, but the sparse-growth behaviour needs careful interpretation.

**Suggested visual**
- [INSERT FIGURE OR TABLE: all-cohort age-bin results]
- Source files:
  - `Results/fba/allcohort_agebin_weighting_scenarios_summary.csv`
  - `Results/figures/allcohort_agebin_vs_equal_flux_delta_top30.png`
  - `Results/figures/allcohort_agebin_top3_species_by_diet.png`
  - `Results/figures/allcohort_agebin_nonzero_species_count_by_diet.png`
- Recommended table columns:
  - `scenario`, `diet`, `community_objective`
- If unavailable: manually summarize age-bin objective values and nonzero-growth trends

**Concrete anchor**
- `Results/fba/allcohort_agebin_weighting_scenarios_summary.csv` shows distinct objective values across bins and diets.

**Interpretation caution**
- Age-bin differences in objective or nonzero growth do not automatically imply stronger biological robustness without follow-up tracing.

---

## Slide 13: Beta-Age Branch — Why It Exists and What to Say Carefully
**Purpose**
Cover the beta-age extension without over-weighting it.

**Key points**
- `08d_prepare_beta_age_weights.py` maps published beta coefficients to modeled species.
- `09c_run_community_fba_beta_age_weighted.py` runs transformed weighting scenarios.
- This branch exists because raw beta coefficients include negative values and cannot be used directly as LP weights.
- It is an exploratory branch testing how transformed coefficient-based weighting changes the solution.

**Suggested table**
- [INSERT TABLE: beta-age branch summary]
- Source files:
  - `Results/inputs/beta_age_input_for_fba.csv`
  - `Results/fba/scenarios/beta_age_positive_only_community_summary_by_diet.csv`
  - `Results/fba/scenarios/beta_age_shifted_community_summary_by_diet.csv`
- Recommended columns:
  - `scenario`, `diet`, `community_objective`
- If unavailable: use 3 bullets summarizing why the branch exists and the two scenario names

**Likely supervisor question**
- Why not use the raw coefficients directly?

**Answer to defend**
- Negative values are incompatible with the intended weighting role in the LP objective, so a transformation step is required before using them as model weights.

---

## Slide 14: COBRApy Pathway Diagnostics — From Reaction Fluxes to Pathway Summaries
**Purpose**
Show how pathway-level interpretation starts in the COBRApy branch.

**Key points**
- `06_flux_diagnostics.py` exports:
  - exchange fluxes
  - connector fluxes
  - reaction-level fluxes
  - pathway activity summaries
- Pathway summaries group active species-internal reactions by pathway annotation.
- These outputs help identify which pathway classes are most active by species and diet.

**Suggested table**
- [INSERT TABLE: pathway activity summary]
- Source files:
  - `Results/fba/community_pathway_activity_by_diet.csv`
  - `Results/fba/community_pathway_diet_comparison.csv`
- Recommended columns:
  - `diet`, `species`, `pathway`, `sum_abs_flux`, `net_flux`, `n_active_reactions`
- If unavailable: create a small manual table using 3 to 5 example rows

**Critical metric definitions**
- `sum_abs_flux`: total absolute flux through pathway-assigned active reactions
- `net_flux`: signed sum of pathway reaction fluxes
- `n_active_reactions`: number of pathway reactions carrying non-negligible flux
- `n_biomass_reactions`: number of biomass-associated reactions in the summarized set

**Likely supervisor question**
- Why use `sum_abs_flux` rather than `net_flux` alone?

**Answer to defend**
- Net flux can cancel internally. `sum_abs_flux` better captures total pathway usage magnitude.

---

## Slide 15: Reaction Tracing — How `06d` Reconstructs Active Pathway Steps
**Purpose**
Explain how pathway summaries are turned into reaction-level traces.

**Key points**
- `06d_trace_pathway_reactions.py` filters active reactions for a chosen species, diet, and pathway.
- It uses:
  - solution flux direction
  - substrate/product orientation
  - shared intermediates between reactions
- Output tables provide a practical reaction ordering, not a claim of a unique biochemical order.

**Suggested visual**
- [INSERT TABLE OR FIGURE: pathway trace example]
- Source files:
  - `Results/fba/pathway_traces/western__Alistipes_shahii_WAL_8301_AGORA1_03__Glycolysis_gluconeogenesis_pipeline.csv`
  - `Results/fba/pathway_traces/western__Alistipes_shahii_WAL_8301_AGORA1_03__Glycolysis_gluconeogenesis_curated_pipeline.png`
- Recommended columns:
  - `step`, `reaction_id`, `reaction_name`, `flux`, `incoming_from`, `shared_metabolites_from_previous`
- If unavailable: copy 3 to 5 reaction steps manually into the slide

**Interpretation caution**
- The trace is a solution-informed route reconstruction, not definitive proof of the only pathway order.

---

## Slide 16: Provenance and Cross-Feeding — How `06f` Connects Metabolites to Biomass Support
**Purpose**
Explain the provenance branch because it is likely to trigger detailed questions.

**Key points**
- `06f_map_species_flux_provenance.py` links:
  - producer species
  - consumer species
  - shared metabolites
  - pathway context
  - biomass-supporting evidence
- `06g_plot_species_crossfeeding_biomass_summary.py` visualizes cleaned provenance summaries.
- This branch helps move from “a metabolite is exchanged” to “this metabolite may support biomass-related metabolism in a specific consumer.”

**Suggested visual**
- [INSERT FIGURE OR TABLE: provenance summary]
- Source files:
  - `Results/fba/provenance/high_fiber__Alistipes_shahii_WAL_8301_AGORA1_03__clean_crossfeeding_biomass_summary.csv`
  - `Results/figures/provenance/alistipes_shahii_crossfeeding_biomass_summary.png`
  - `Results/figures/provenance/faecalibacterium_prausnitzii_crossfeeding_biomass_summary.png`
- Recommended columns:
  - `producer_species`, `shared_metabolite_name`, `consumer_species`, `crossfeeding_support_score`, `biomass_association_score`, `confidence`
- If unavailable: build a manual top-3 candidate table

**Critical metric definitions**
- `crossfeeding_support_score`: heuristic ranking based on overlap between secretion and uptake support for a shared metabolite
- `biomass_association_score`: heuristic ranking linking a shared metabolite to biomass-associated pathway support in the consumer
- `provenance edges`: reaction-to-reaction links supported by shared intermediates in the traced solution

**Likely supervisor question**
- What exactly does `biomass_association_score` mean?

**Answer to defend**
- It is a workflow-specific ranking metric, not a universal biological constant. It summarizes how strongly the current solution links a shared metabolite to biomass-associated pathway support in the consumer species.

---

## Slide 17: Why MICOM Was Introduced After Weighted FBA
**Purpose**
Justify the transition to MICOM.

**Key points**
- Weighted LP community FBA can produce winner-take-all behaviour.
- MICOM was added to reduce that behaviour and distribute growth more realistically across community members.
- `13_prepare_micom_inputs.py` converts abundance tables into MICOM-ready community inputs.
- `14_run_micom_by_agegroup.py` builds one MICOM community per age group and diet.

**Suggested visual**
- [INSERT COMPARISON DIAGRAM: weighted LP vs MICOM]
- Source references:
  - `Scripts/13_prepare_micom_inputs.py`
  - `Scripts/14_run_micom_by_agegroup.py`
  - `README.md`
- If unavailable: manually create a two-column comparison of optimization logic

**Definitions to include**
- `cooperative tradeoff`: maximize community growth while allowing controlled redistribution across taxa
- `tradeoff fraction`: fraction of maximal community growth retained during cooperative tradeoff

---

## Slide 18: MICOM Main Outputs — Growth and Tradeoff Sensitivity
**Purpose**
Show what MICOM adds at the community-growth level.

**Key points**
- Main MICOM outputs include:
  - community growth summary by age group and diet
  - organism growth rates
  - growth contribution rankings
  - tradeoff sensitivity summaries
- `14d_tradeoff_sensitivity_micom.py` tests fractions from `0.0` to `1.0`.
- In the sampled outputs, tradeoff fraction changes community growth substantially, while the fraction of growing taxa can remain stable in some age-group/diet blocks.

**Suggested visual**
- [INSERT FIGURE OR TABLE: MICOM tradeoff sensitivity]
- Source files:
  - `Results/micom/growth/proper_age_bins/community_growth_summary_by_agegroup_diet.csv`
  - `Results/micom/growth/proper_age_bins/organism_growth_rates_by_agegroup_diet.csv`
  - `Results/micom/tradeoff_sensitivity/proper_age_bins/tradeoff_sensitivity_summary.csv`
  - `Results/figures/micom/growth/proper_age_bins/micom_tradeoff_community_growth_curves_proper_age_bins.png`
  - `Results/figures/micom/growth/proper_age_bins/micom_tradeoff_species_growth_curves_proper_age_bins.png`
- Recommended table columns:
  - `tradeoff_fraction`, `age_group`, `diet`, `community_growth_rate`, `fraction_growing_in_diet`
- If unavailable: manually plot or tabulate 3 to 4 representative rows

**Concrete anchor**
- `Results/micom/tradeoff_sensitivity/proper_age_bins/tradeoff_sensitivity_summary.csv` shows strong growth-rate changes across fractions.

**Likely supervisor question**
- Why use `0.5` as a main setting?

**Answer to defend**
- It is a compromise value examined in the sensitivity sweep, balancing retained community growth with broader taxon participation.

---

## Slide 19: MICOM Pathway Analysis and Lysine/Butyrate Case Study
**Purpose**
Show the strongest pathway-centered part of the MICOM branch.

**Key points**
- `16_find_differential_pathways.py` collapses reaction-level MICOM flux into pathway summaries and ranks pathways that differ across age groups.
- `15_plot_lysine_pathways.py` focuses on:
  - lysine biosynthesis
  - candidate lysine-to-butyrate routes
- Curated reaction IDs are used instead of broad keyword matching to reduce false positives.
- Current interpretation note:
  - downstream butanoate-related steps can show nonzero flux even when upstream lysine-specific steps remain zero

**Suggested visual**
- [INSERT FIGURE OR TABLE: MICOM pathway and lysine case study]
- Source files:
  - `Results/micom/pathway_flux/differential_pathways_by_agegroup_diet.csv`
  - `Results/figures/micom/pathway_flux/differential_pathways_heatmap_top25.png`
  - `Results/micom/lysine_butyrate/lysine_to_butyrate_reaction_step_summary.csv`
  - `Results/figures/micom/lysine_butyrate/lysine_to_butyrate_reaction_steps_by_agegroup_diet.png`
  - `Results/figures/micom/lysine_butyrate/lysine_to_butyrate_flux_by_agegroup_diet.png`
- Recommended table columns:
  - pathway table: `diet`, `pathway`, `range_sum_abs_flux`, `std_sum_abs_flux`
  - lysine table: `age_group`, `diet`, `reaction_step`, `sum_abs_flux`
- If unavailable: manually summarize the top 3 differential pathways and 3 lysine-route steps

**Critical metric definitions**
- pathway `sum_abs_flux`: total absolute pathway flux after grouping by age group and diet
- `range_sum_abs_flux`: max minus min pathway activity across age groups within a diet
- `std_sum_abs_flux`: variability of pathway activity across age groups within a diet

**Likely supervisor question**
- Why is the lysine-to-butyrate interpretation cautious?

**Answer to defend**
- Because downstream butanoate activity alone does not prove the full lysine-derived route is active. Reaction-level support for upstream lysine-specific steps is still required.

---

## Slide 20: Current Status, Main Limitations, and Next Steps
**Purpose**
End with a realistic and defensible project position.

**Key points**
- Completed branches:
  - baseline COBRApy
  - SG90 weighted FBA
  - all-cohort age-bin FBA
  - beta-age exploratory branch
  - MICOM proper-age-bin branch
  - pathway tracing and provenance summaries
- Main limitations:
  - sparse optima in weighted LP FBA
  - pathway annotation granularity
  - provenance and biomass-link metrics are heuristic, not causal proof
  - pathway-level summaries need reaction-level confirmation for strong mechanistic claims
- Next steps:
  - finalize validation and QA
  - tighten reporting around canonical outputs
  - deepen targeted mechanistic checks for specific metabolites and pathways

**Suggested visual**
- [INSERT SUMMARY SLIDE VISUAL]
- Source references:
  - `README.md`
  - `Scripts/99_validate_pipeline_outputs.py`
- If unavailable: use a two-column manually created slide: `what is established` vs `what still needs caution`

**Speaker notes**
- End on a strong but conservative message: the project now has both community-level and pathway-level views, but the strongest mechanistic claims still require reaction-level scrutiny.

---

## Appendix: Metric Cheat Sheet for Backup or Speaker Notes
- `community_objective`:
  optimized objective value of the community model under a given scenario and diet
- `biomass_flux`:
  flux through a species biomass reaction, used as a growth proxy
- `sum_abs_flux`:
  total absolute flux through a grouped set of reactions
- `net_flux`:
  signed sum of grouped reaction fluxes
- `n_active_reactions`:
  number of grouped reactions with non-negligible flux
- `crossfeeding_support_score`:
  heuristic score ranking producer-consumer support for a shared metabolite
- `biomass_association_score`:
  heuristic score ranking how strongly a shared metabolite links to biomass-associated pathway support
- `tradeoff_fraction`:
  fraction of maximal community growth retained in MICOM cooperative tradeoff
