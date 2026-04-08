# Speaker Notes: Supervisor Q&A Defense Guide

Use this with `slidedeck.md`. The goal is not to memorize every number, but to know:
- which file to cite
- which columns matter
- what a vague column name actually means
- what not to overclaim

## General rule for answers
- Start with what the metric is operationally in this project.
- Then point to the source CSV.
- Then say what to look for in the key columns.
- End with one caution if the metric is heuristic or summary-level.

---

## Slide 5: QC and Single-Species FBA
**Likely question**
- How do you know the models were usable before building communities?

**Strong answer**
- I checked model loading and basic model structure first in `Results/qc/model_qc_summary.csv`.
- The useful columns are:
  - `file`: which SBML file was loaded
  - `model_id`: model identifier inside the SBML
  - `reactions`, `metabolites`, `genes`: basic model size checks
  - `objective`: what the model is optimizing, usually a biomass-like objective
- If a row has an error instead of those fields, that model failed QC and should not be trusted downstream.
- I then checked single-species growth in `Results/fba/single_species_growth_default.csv` and `Results/fba/single_species_growth_by_diet.csv` to confirm whether species could grow before any community interaction was introduced.

**What to look out for**
- A model can load successfully but still have low or zero growth under a given diet.
- That is different from a model failing to load at QC stage.

---

## Slide 6: Community Construction
**Likely question**
- Why did you build a shared pool and connector reactions instead of just merging all models together?

**Strong answer**
- The shared pool preserves species identity while allowing metabolite exchange between species.
- The logic is implemented in `Scripts/04_run_community_fba.py`, and the structural summary can be supported with `Results/fba/community_model_summary.csv`.
- In this project, connector reactions are the reversible links between a species extracellular metabolite and the shared community metabolite pool.
- That matters because it lets me distinguish:
  - species-specific internal metabolism
  - shared extracellular cross-feeding
- If I merged all species reactions without that structure, I would lose clean attribution of who secreted or consumed what.

**What to look out for**
- Structural counts are diet-independent.
- Diet changes the optimized flux solution, not the community model size itself.

---

## Slide 7: Full-Access Diagnostic Logic
**Likely question**
- How do you separate model incapability from competition or diet limitation?

**Strong answer**
- I use the full-access branch as a diagnostic control, not as the main biological result.
- The relevant files are:
  - `Results/fba/community_growth_full_access.csv`
  - `Results/fba/community_species_biomass_flux_full_access.csv`
  - `Results/fba/individual_growth_full_access_medium.csv`
- If a species can grow in full-access but not in the diet-constrained community run, that suggests suppression by diet limitation or competition rather than outright incapability.

**What to look out for**
- Because full-access opens many exchanges strongly, it removes much of the original diet contrast.
- So it is a capability check, not a realistic ecological condition.

---

## Slide 8: Weighted FBA Rationale
**Likely question**
- Why did you move from equal-abundance community FBA to weighted community FBA?

**Strong answer**
- Equal-abundance FBA assumes all taxa influence the community objective equally.
- In this project I wanted the objective to reflect abundance structure, so weighted community FBA assigns different biomass weights to different species.
- The main supporting files are:
  - `Results/inputs/sg90_median_input_for_fba.csv`
  - `Results/fba/scenarios/sg90_median_community_summary_by_diet.csv`
  - `Results/fba/scenarios/equal_abundance_community_summary_by_diet.csv`
- The comparison is not just conceptual. I can inspect differences in:
  - `community_objective`
  - species biomass rankings
  - exchange flux rankings

**What to look out for**
- Zero growth in a weighted LP solution does not automatically mean the model is broken.
- It can reflect sparse optima.

---

## Slide 9: SG90 Inputs and Scenarios
**Likely question**
- What exactly do the scenario files mean?

**Strong answer**
- The scenario files separate different objective assumptions.
- The key files are under `Results/fba/scenarios/`.
- For example:
  - `sg90_median_community_summary_by_diet.csv`
  - `equal_abundance_community_summary_by_diet.csv`
  - `sg90_median_species_biomass_by_diet.csv`
- Important columns:
  - `scenario`: which weighting assumption was used
  - `diet`: `western` or `high_fiber`
  - `status`: whether optimization solved
  - `community_objective`: optimized objective value under that scenario
  - `applied_exchanges`: number of medium exchanges successfully applied

**Plain-language explanation**
- `community_objective` is not “total biology” or “true health.”
- It is the optimized value of the mathematical objective under that specific model setup.

---

## Slide 10: SG90 Weighted Results
**Likely question**
- How do you justify saying weighting changes the result?

**Strong answer**
- I compare weighted and equal-abundance outputs directly using:
  - `Results/fba/scenarios/sg90_median_vs_equal_exchange_flux_comparison_by_diet.csv`
  - `Results/figures/sg90_vs_equal_flux_delta_top25.png`
- For biomass outputs I check:
  - `Results/fba/scenarios/sg90_median_species_biomass_by_diet.csv`
- The point is not just that absolute values differ.
- The stronger claim is that changing the objective weights changes which species dominate and which community exchange fluxes become large.

**What to look out for**
- Do not claim pathway mechanism from these comparison files alone.
- These are still objective-level and exchange-level summaries.

---

## Slide 12: All-Cohort Age-Bin Results
**Likely question**
- What exactly is changing across age bins?

**Strong answer**
- The clean starting point is `Results/fba/allcohort_agebin_weighting_scenarios_summary.csv`.
- The important columns are:
  - `scenario`: which age-bin weighting scenario
  - `diet`
  - `community_objective`
  - `applied_exchanges`
- Example from the current file:
  - `allcohort_agebin_71_80` has a higher `community_objective` than `allcohort_agebin_81_90` in both diets
- Then I check:
  - `Results/figures/allcohort_agebin_top3_species_by_diet.png`
  - `Results/figures/allcohort_agebin_nonzero_species_count_by_diet.png`
- That lets me ask both:
  - how strong the optimized community objective is
  - how many taxa actually receive nonzero growth

**Plain-language explanation**
- `community_objective` here means the optimized weighted community growth target, not a direct biological measurement from patients.

**What to look out for**
- Some age bins show sparse solutions with only a few growing taxa.
- That is an optimization outcome, not proof that only those taxa matter biologically.

---

## Slide 13: Beta-Age Branch
**Likely question**
- Why is the beta-age branch exploratory?

**Strong answer**
- Because the starting coefficients are not directly ready for the LP objective.
- The relevant files are:
  - `Results/inputs/beta_age_input_for_fba.csv`
  - `Results/fba/scenarios/beta_age_positive_only_community_summary_by_diet.csv`
  - `Results/fba/scenarios/beta_age_shifted_community_summary_by_diet.csv`
- The main issue is that raw beta coefficients can be negative, and negative raw objective weights are not appropriate for the role they are supposed to play here.
- So the branch is about transformed weighting scenarios, not a final biological conclusion.

**What to look out for**
- Present it as a sensitivity or exploratory branch.
- Do not present it as equally canonical to SG90 or all-cohort age-bin.

---

## Slide 14: Pathway Diagnostics Columns
**Likely question**
- What do `sum_abs_flux`, `net_flux`, and `n_active_reactions` actually mean?

**Strong answer**
- The source file is `Results/fba/community_pathway_activity_by_diet.csv`.
- The important columns are:
  - `diet`
  - `species`
  - `pathway`
  - `sum_abs_flux`
  - `net_flux`
  - `n_active_reactions`
  - `n_biomass_reactions`

**Plain-language column explanations**
- `sum_abs_flux`:
  - add up the absolute values of all active reaction fluxes in that pathway
  - this is a pathway activity magnitude score
  - it tells me how much pathway traffic exists overall
- `net_flux`:
  - signed sum of fluxes in that pathway
  - this gives a rough directional tendency
  - it can be small even when the pathway is busy because positive and negative contributions can cancel
- `n_active_reactions`:
  - number of pathway reactions carrying non-negligible flux in the solution
- `n_biomass_reactions`:
  - number of biomass-associated reactions counted in that summarized set

**How to answer “why not just use net flux?”**
- Because `net_flux` can hide activity through cancellation.
- `sum_abs_flux` is better for pathway usage magnitude.
- `net_flux` is complementary, not sufficient by itself.

**What to look out for**
- Pathway labels depend on model annotations.
- High pathway activity is not the same thing as proof of a full biochemical route.

---

## Slide 15: Reaction Tracing
**Likely question**
- How do you move from pathway labels to actual reactions?

**Strong answer**
- I use `Results/fba/pathway_traces/...pipeline.csv` files generated by `06d_trace_pathway_reactions.py`.
- For example:
  - `Results/fba/pathway_traces/western__Alistipes_shahii_WAL_8301_AGORA1_03__Glycolysis_gluconeogenesis_pipeline.csv`
- Important columns are:
  - `step`
  - `reaction_id`
  - `reaction_name`
  - `flux`
  - `equation_in_flux_direction`
  - `incoming_from`
  - `shared_metabolites_from_previous`

**Plain-language column explanations**
- `equation_in_flux_direction`:
  - the reaction written in the direction the optimized solution actually used
- `incoming_from`:
  - which earlier reaction connects into this one through shared intermediates
- `shared_metabolites_from_previous`:
  - the metabolite links used to stitch the trace together

**What to look out for**
- This is a practical route reconstruction from the solution.
- It is not a claim that biology must follow exactly that one order outside the model.

---

## Slide 16: Provenance and Cross-Feeding Metrics
**Likely question**
- What do `crossfeeding_support_score` and `biomass_association_score` mean?

**Strong answer**
- Use `Results/fba/provenance/high_fiber__Alistipes_shahii_WAL_8301_AGORA1_03__clean_crossfeeding_biomass_summary.csv` as the worked example.
- Important columns are:
  - `producer_species`
  - `shared_metabolite_name`
  - `consumer_species`
  - `uptake_flux`
  - `producer_flux`
  - `crossfeeding_support_score`
  - `consumer_top_pathway`
  - `pathway_flux_share`
  - `biomass_association_score`
  - `evidence_type`
  - `confidence`

**Plain-language column explanations**
- `crossfeeding_support_score`:
  - a project-specific support score for a producer-consumer-metabolite link
  - operationally, it reflects how much producer secretion and consumer uptake overlap in the solved community
  - bigger values mean stronger candidate support in this solution
- `consumer_top_pathway`:
  - the main consumer pathway linked to the metabolite class in the workflow summary
- `pathway_flux_share`:
  - how much of the consumer’s pathway-level activity is represented by that linked pathway signal
  - it is a share-like contextual number, not a standalone biological truth
- `biomass_association_score`:
  - a project-specific heuristic linking a shared metabolite to biomass-associated pathway support in the consumer
  - bigger values mean stronger workflow evidence that the metabolite may support biomass-related metabolism
- `evidence_type`:
  - short label describing why the row was retained, such as flux support plus pathway-class matching
- `confidence`:
  - qualitative confidence tier assigned by the workflow

**Best short answer if pressed**
- These are ranking heuristics built from the solved flux pattern and pathway linkage logic.
- They are useful for prioritizing mechanistic hypotheses, not for making direct causal claims.

**What to look out for**
- If `crossfeeding_support_score` is high but `biomass_association_score` is weak, that suggests exchange evidence without strong biomass-link evidence.
- If both are stronger, that row is a better follow-up candidate.

---

## Slide 18: MICOM Tradeoff Sensitivity
**Likely question**
- What does `tradeoff_fraction` actually do, and how do you know `0.5` is reasonable?

**Strong answer**
- The source file is `Results/micom/tradeoff_sensitivity/proper_age_bins/tradeoff_sensitivity_summary.csv`.
- Important columns are:
  - `tradeoff_fraction`
  - `age_group`
  - `diet`
  - `community_growth_rate`
  - `n_organisms`
  - `n_growing_in_diet`
  - `fraction_growing_in_diet`
  - `applied_exchanges`
  - `missing_exchanges`
  - `note`

**Plain-language column explanations**
- `tradeoff_fraction`:
  - the fraction of the maximum community growth retained while redistributing growth across taxa
- `community_growth_rate`:
  - the optimized MICOM community growth under that tradeoff setting
- `n_growing_in_diet`:
  - count of organisms with nonzero growth above the project threshold
- `fraction_growing_in_diet`:
  - `n_growing_in_diet / n_organisms`
- `missing_exchanges`:
  - diet exchanges from the input medium that did not map onto MICOM community exchanges
- `note`:
  - records if a fallback run was needed, for example rerunning with relaxed minimum growth

**How to defend `0.5`**
- I did not pick it blindly.
- I tested the whole sweep from `0.0` to `1.0`.
- In the current file, for example in `21_40` high-fiber, `community_growth_rate` rises strongly as the fraction increases, while `fraction_growing_in_diet` stays at `0.7`.
- So `0.5` is a middle-ground setting that retains substantial growth without claiming it is uniquely optimal.

**What to look out for**
- A stable `fraction_growing_in_diet` across fractions does not mean the tradeoff is irrelevant.
- Growth magnitude can still change a lot.

---

## Slide 19: MICOM Differential Pathways and Lysine Case
**Likely question**
- What do `range_sum_abs_flux` and `std_sum_abs_flux` mean, and why is the lysine-to-butyrate story cautious?

**Strong answer**
- For pathway differences across age groups, use:
  - `Results/micom/pathway_flux/differential_pathways_by_agegroup_diet.csv`
- Important columns:
  - `diet`
  - `pathway`
  - `max_sum_abs_flux`
  - `min_sum_abs_flux`
  - `range_sum_abs_flux`
  - `std_sum_abs_flux`

**Plain-language column explanations**
- `max_sum_abs_flux`:
  - largest pathway activity seen across age groups within that diet
- `min_sum_abs_flux`:
  - smallest pathway activity seen across age groups within that diet
- `range_sum_abs_flux`:
  - `max - min`
  - used here as the main difference score across age groups
- `std_sum_abs_flux`:
  - standard deviation of pathway activity across age groups
  - another way to describe variability

**How to defend the lysine caution**
- Use:
  - `Results/micom/lysine_butyrate/lysine_to_butyrate_reaction_step_summary.csv`
  - `Results/figures/micom/lysine_butyrate/lysine_to_butyrate_reaction_steps_by_agegroup_diet.png`
- The caution is that downstream butanoate-related steps can show nonzero flux even if the upstream lysine-specific steps remain zero.
- So I should not claim a full lysine-derived butyrate route is active unless the reaction-step table supports upstream lysine-specific activity too.

**What to look out for**
- Broad pathway labels like `Lysine metabolism` are weaker evidence than curated reaction-step summaries.
- Reaction-level evidence is stronger than pathway-name-only evidence.

---

## Short Answers for Common Pushback
**“Is this metric a real biological quantity?”**
- Usually not directly.
- Many of these columns are workflow summaries or heuristic ranking scores built from solved model fluxes.

**“What is the safest file to show for pathway summaries?”**
- Start with `Results/fba/community_pathway_activity_by_diet.csv` for COBRApy.
- Start with `Results/micom/pathway_flux/differential_pathways_by_agegroup_diet.csv` for MICOM.
- Then drop to reaction-level files if the question becomes mechanistic.

**“What is the safest file to show for cross-feeding?”**
- Use the cleaned provenance summary under `Results/fba/provenance/*clean_crossfeeding_biomass_summary.csv`.
- It is already filtered into a more interpretable summary than raw connector tables.

**“What is the biggest overclaim risk in this project?”**
- Treating pathway-level activity or heuristic scores as direct causal biological proof.
- The safe language is:
  - supports
  - suggests
  - prioritizes
  - is consistent with
- Avoid:
  - proves
  - confirms
  - demonstrates definitively

---

## Quick File-to-Question Map
- Model usability:
  - `Results/qc/model_qc_summary.csv`
- Single-species growth:
  - `Results/fba/single_species_growth_default.csv`
  - `Results/fba/single_species_growth_by_diet.csv`
- Community and weighted FBA summaries:
  - `Results/fba/scenarios/*_community_summary_by_diet.csv`
- Age-bin comparison:
  - `Results/fba/allcohort_agebin_weighting_scenarios_summary.csv`
- Pathway summaries:
  - `Results/fba/community_pathway_activity_by_diet.csv`
- Reaction traces:
  - `Results/fba/pathway_traces/*_pipeline.csv`
- Cross-feeding and biomass support:
  - `Results/fba/provenance/*clean_crossfeeding_biomass_summary.csv`
- MICOM tradeoff defense:
  - `Results/micom/tradeoff_sensitivity/proper_age_bins/tradeoff_sensitivity_summary.csv`
- MICOM pathway difference defense:
  - `Results/micom/pathway_flux/differential_pathways_by_agegroup_diet.csv`
- Lysine/butyrate caution:
  - `Results/micom/lysine_butyrate/lysine_to_butyrate_reaction_step_summary.csv`
