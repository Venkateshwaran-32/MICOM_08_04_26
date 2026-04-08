# Pipeline Story: What The Results Actually Show

## What This Project Is Trying To Answer
The project started with a broad question:

How does gut microbial community metabolism change with diet and aging, and which species, pathways, and metabolites matter most?

The important part is that the repo did not answer this in one jump. Each branch was added because the previous branch gave a partial answer but also exposed a limitation.

The cleanest mental model is:

`QC/models -> single-species COBRApy -> baseline community COBRApy -> weighted community COBRApy -> all-cohort age-bin COBRApy -> MICOM proper-age-bins -> pathway interpretation -> provenance -> future subject-level branch`

What follows is the same story as the old pipeline guide, but grounded in the actual current outputs.

## Stage 1: Were The Models Even Usable?
This is the purpose of the earliest COBRApy checks.

Main scripts:
- `Scripts/01_load_qc_models.py`
- `Scripts/02_single_species_fba_default.py`
- `Scripts/03_single_species_fba_diets.py`

Canonical files:
- `Results/qc/model_qc_summary.csv`
- `Results/fba/single_species_growth_default.csv`
- `Results/fba/single_species_growth_by_diet.csv`

What the results show:
- All 10 AGORA models in this project have positive single-species growth under both `western` and `high_fiber`.
- In `Results/fba/single_species_growth_by_diet.csv`, all 10 species have `growth_western > 0` and all 10 have `growth_high_fiber > 0`.
- Some species grow better than others as single models. For example:
  - `Parabacteroides merdae` is the strongest single-species grower in both diets at about `0.862`.
  - `E. coli` and `Klebsiella` are also strong as single models.
  - `Faecalibacterium`, `Alistipes shahii`, and `Ruminococcus torques` do grow, but at lower levels.

What you are allowed to conclude:
- the SBML models load
- the model set is technically usable
- later zero-growth in communities is not because the species models are intrinsically dead

Why the project moved on:
- single-species growth says nothing about competition, sharing, or community structure
- the next step had to ask what happens when all taxa are put into one community objective

## Stage 2: What Happened In The First Community COBRApy Run?
This is the baseline unweighted community branch.

Main script:
- `Scripts/04_run_community_fba.py`

Canonical file:
- `Results/fba/community_growth_by_diet.csv`

What the results show:
- The community model solves in both diets.
- The total community objective is:
  - `4.296417` in `western`
  - `4.168428` in `high_fiber`
- `130` exchanges were applied and `34` were missing in both diets.

What this stage proved:
- the community model construction works
- the project can optimize one shared community objective under both diets

What this stage did not answer well enough:
- this table is mainly a top-level objective summary
- it does not yet tell a satisfying story about which species are winning and whether the solution is biologically plausible
- plain community FBA is known to push biomass into a few winners, so this stage was not enough for interpretation

Why the project moved on:
- the next question became whether changing the objective weights would change who wins

## Stage 3: Did Abundance Weighting Change The Community Answer?
This is the equal-abundance versus SG90-weighted COBRApy branch.

Main scripts:
- `Scripts/08_prepare_abundance_inputs.py`
- `Scripts/09_run_community_fba_weighted.py`
- `Scripts/10_compare_weighted_vs_equal_fluxes.py`
- `Scripts/11_plot_weighted_vs_equal.py`

Canonical files:
- `Results/fba/scenarios/equal_abundance_species_biomass_by_diet.csv`
- `Results/fba/scenarios/sg90_median_species_biomass_by_diet.csv`

What the equal-abundance results show:
- The equal-abundance COBRApy solution is very sparse.
- In both diets, `Alistipes shahii` dominates with biomass flux about `3.694`.
- `Faecalibacterium prausnitzii` is the only other clearly positive species:
  - about `0.474` in `high_fiber`
  - about `0.602` in `western`
- The remaining taxa are zero or numerical near-zero.

What the SG90-weighted results show:
- Weighting changes the winners.
- In both diets, `E. coli` becomes the top grower at `2.562616`.
- `Faecalibacterium prausnitzii` becomes second:
  - `0.732950` in `high_fiber`
  - `0.805065` in `western`
- Most other taxa remain zero or effectively zero.

What you are allowed to conclude:
- abundance weighting materially changes which species dominate the COBRApy solution
- the identity of the winners is sensitive to the weighting scheme
- SG90 weighting is scientifically meaningful because it changes the answer, not just the labels

What you are not allowed to conclude:
- that weighting solved the realism problem
- that a sparse optimum with two winners is automatically a faithful biological community

Why the project moved on:
- the repo now had a better weighted branch, but it was still too sparse
- the next question became whether age structure changes the answer in a meaningful way

## Stage 4: What Changed When Age Bins Were Added In COBRApy?
This is the all-cohort age-bin weighted COBRApy branch.

Main scripts:
- `Scripts/08b_prepare_allcohort_agebin_inputs.py`
- `Scripts/09b_run_community_fba_allcohort_agebins.py`
- `Scripts/10b_compare_allcohort_agebins.py`
- `Scripts/11b_plot_allcohort_agebins.py`
- `Scripts/11c_plot_allcohort_top3_species.py`
- `Scripts/11d_plot_nonzero_and_log_biomass.py`

Canonical files:
- `Results/fba/allcohort_agebin_weighting_scenarios_summary.csv`
- `Results/fba/allcohort_agebin_top3_species_by_diet.csv`

What the results show:
- The age-bin branch changes the winners by age group.
- In `61_70`, `Faecalibacterium prausnitzii` is the clear top grower:
  - `1.407890` in `high_fiber`
  - `1.546083` in `western`
- In `71_80`, `E. coli` becomes the top grower:
  - `2.562616` in both diets
  - `Faecalibacterium` becomes second
- In `81_90`, `E. coli` remains first and `Faecalibacterium` remains second.
- Community objectives also differ across bins:
  - `61_70`: about `1.13` to `1.24`
  - `71_80`: about `1.32` to `1.34`
  - `81_90`: about `0.80` to `0.82`

What this stage added to the story:
- older age bins are not just a repeat of the SG90 median result
- the >60 bins are where the species pattern becomes most interesting
- the switch from `Faecalibacterium` dominance in `61_70` to `E. coli` dominance in `71_80` and `81_90` becomes a real signal worth following

What still remained unsatisfying:
- the age-bin COBRApy solutions are still sparse
- the older-bin pattern might be real, but the LP objective still makes it hard to know whether the pattern is dominated by winner-take-all optimization

Why the project moved on:
- MICOM was added to keep the age-bin question but relax the extreme sparsity of weighted COBRApy

## Stage 5: What Did MICOM Change?
This is the MICOM proper-age-bins branch, which is now the canonical MICOM branch.

Main scripts:
- `Scripts/13_prepare_micom_inputs.py`
- `Scripts/14_run_micom_by_agegroup.py`
- `Scripts/14b_rank_micom_growth_contributions.py`
- `Scripts/14d_tradeoff_sensitivity_micom.py`

Canonical files:
- `Results/micom/growth/proper_age_bins/organism_growth_rates_by_agegroup_diet.csv`
- `Results/micom/tradeoff_sensitivity/proper_age_bins/tradeoff_sensitivity_summary.csv`

Why MICOM mattered:
- The key issue with COBRApy weighted FBA was not that it found nothing.
- The issue was that it found age-bin patterns inside solutions that were still too sparse.
- MICOM keeps the same biological ingredients but uses cooperative tradeoff, which distributes growth across more taxa.

What the tradeoff results show at the canonical `fraction = 0.5`:
- `21_40`: `7/10` taxa grow
- `41_60`: `5/10` taxa grow
- `61_70`: `5/10` taxa grow
- `71_80`: `7/10` taxa grow
- `81_plus`: `9/10` taxa grow

This is the crucial result:
- in the oldest bin, almost the whole modeled community receives positive growth
- that is a very different community structure from the sparse COBRApy optima

What the MICOM growth results show above 60:
- `61_70`
  - `Faecalibacterium prausnitzii` is still the top grower
  - `Ruminococcus torques` is second
  - `E. coli` is present but much lower
- `71_80`
  - `E. coli` becomes the top grower
  - `Faecalibacterium` is second
  - `Alistipes onderdonkii` and `Bacteroides dorei` also grow strongly
- `81_plus`
  - `E. coli` remains first
  - `Faecalibacterium` remains second
  - `Parabacteroides merdae` and `Bacteroides dorei` become substantial additional growers

What you are allowed to conclude:
- the older-age signal seen in COBRApy is not just erased by switching frameworks
- MICOM keeps the same broad older-bin story, especially the importance of `E. coli` and `Faecalibacterium`
- MICOM also makes the community more distributed, which makes later pathway interpretation more credible

What you are not allowed to conclude:
- that MICOM is automatically “correct” and COBRApy is “wrong”
- that MICOM numbers should be directly compared to COBRApy biomass flux values as if they were the same optimization problem

Why the project moved on:
- once the older-age MICOM pattern looked more stable and interpretable, the next question became which pathways and metabolites might explain it

## Stage 6: What Did The Lysine / Butanoate Work Actually Find?
This is an interpretation layer on top of MICOM reaction flux outputs.

Main scripts:
- `Scripts/15_plot_lysine_pathways.py`
- `Scripts/15b_plot_lysine_agegroup_metabolism.py`
- `Scripts/15c_compare_lysine_vs_butanoate.py`

Canonical files:
- `Results/micom/lysine_butyrate/lysine_biosynthesis_flux_summary.csv`
- `Results/micom/lysine_butyrate/lysine_to_butyrate_flux_summary.csv`
- `Results/micom/lysine_butyrate/lysine_to_butyrate_reaction_step_summary.csv`

What the results show:
- `Lysine metabolism` pathway-level flux is present across age bins.
- For example, in `lysine_biosynthesis_flux_summary.csv`, `Lysine metabolism` is strongly nonzero across bins and diets.
- But the curated lysine-to-butyrate summary is much narrower.
- In `lysine_to_butyrate_flux_summary.csv`, both the `Butanoate metabolism` and `Lysine metabolism` entries are `0` through `61_70`.
- In `lysine_to_butyrate_reaction_step_summary.csv`, downstream butanoate steps become nonzero only in older groups:
  - `71_80 high_fiber`: `46.607253`
  - `71_80 western`: `37.584627`
  - older bins continue to show nonzero downstream butanoate activity

What this means scientifically:
- the current evidence does not support a strong claim that the full lysine-to-butyrate route is active from start to finish
- what it does support is a narrower claim:
  - older MICOM bins show downstream butanoate-related activity
  - but upstream lysine-specific proof is still missing in the curated step summaries

Why this mattered:
- this is exactly the point where the project shifted from “which taxa grow?” to “what mechanism is worth tracing more carefully?”
- that question naturally led to provenance

## Stage 7: What Is Provenance Actually Explaining?
Provenance is the attempt to connect a growth-relevant pathway signal to a plausible supporting metabolite and then ask where that metabolite came from.

There are two versions in this repo:
- COBRApy provenance
- MICOM provenance pilot

### COBRApy provenance
Main scripts:
- `Scripts/06d_trace_pathway_reactions.py`
- `Scripts/06f_map_species_flux_provenance.py`
- `Scripts/06g_plot_species_crossfeeding_biomass_summary.py`

Canonical outputs:
- `Results/fba/provenance/*`
- `Results/figures/provenance/*`

Why this branch is stronger mechanistically:
- the COBRApy community has an explicit shared-pool connector structure
- that lets the provenance outputs trace producer -> shared metabolite -> consumer more directly

What the current COBRApy provenance examples show:
- `Alistipes shahii` under `high_fiber`
  - top biomass-adjacent pathways include glycolysis and nucleotide interconversion
  - the strongest cross-feeding candidate retained is `Adenosine` from `Parabacteroides merdae`
- `Faecalibacterium prausnitzii` under `western`
  - analogous summary outputs exist for its connector and pathway story

What this branch is trying to prove:
- not “this single metabolite caused biomass”
- but “this is a plausible metabolite-level support story for the pathway state associated with biomass”

### MICOM provenance pilot
Main scripts:
- `Scripts/16b_map_micom_flux_provenance.py`
- `Scripts/16c_plot_micom_provenance_summary.py`

Canonical file:
- `Results/micom/provenance/proper_age_bins/micom_provenance_case_overview.csv`

What the current MICOM pilot shows above 60:
- The selected focal species are `E. coli` and `Faecalibacterium`, because they are the most important older-bin MICOM growers.
- For all 12 current MICOM pilot cases, the top pathway is `Nucleotide interconversion`.

For `E. coli`:
- `61_70`: anchor metabolite is `Adenosine`
- `71_80` and `81_plus`: anchor metabolite is `2-deoxyadenosine`
- the likely origin is usually `mixed_medium_plus_crossfeeding`
- the top producer changes with context:
  - `Parabacteroides merdae` in `61_70`
  - `Bacteroides xylanisolvens` in `71_80`
  - `Alistipes shahii` in `81_plus`

For `Faecalibacterium`:
- the anchor metabolite is consistently `Hypoxanthine`
- `61_70`: mixed medium plus crossfeeding, with `Parabacteroides merdae` as top producer
- `71_80`: direct from medium
- `81_plus`: mixed again, with producer support appearing in the oldest bin

What you are allowed to conclude:
- older-bin MICOM growth for these species is associated with a consistent nucleotide-interconversion signal
- there are plausible external metabolite candidates that may support that pathway state
- the source of support can shift between medium-only and mixed medium-plus-crossfeeding

What you are not allowed to conclude:
- that MICOM provenance proves direct transfer between named species
- that the anchor metabolite alone explains the whole growth phenotype

The practical difference is:
- COBRApy provenance is closer to a mechanistic connector trace
- MICOM provenance is an inferred support map

## What The Repo Has Actually Found So Far
If you ignore script numbers and look only at the results, the current story is:

1. All 10 species models are individually viable under both diets.
2. Plain and weighted COBRApy community FBA are both sparse, but the winning taxa depend strongly on the objective weights.
3. In age-bin COBRApy, the most interesting shift happens above 60:
   - `61_70` is `Faecalibacterium`-led
   - `71_80` and `81_90` become `E. coli`-led with `Faecalibacterium` second
4. MICOM preserves that >60 pattern while distributing growth across many more taxa, especially in `81_plus`.
5. The most important older-bin MICOM species right now are `E. coli` and `Faecalibacterium`, with other taxa becoming relevant in the oldest bin.
6. The lysine/butanoate work is suggestive but not a full mechanistic proof:
   - downstream butanoate activity appears in older bins
   - upstream lysine-specific confirmation is still missing
7. Provenance is now being used to explain which external metabolite may support a top biomass-adjacent pathway and whether that support looks medium-derived, cross-fed, or mixed.

## Why Subject-Level FBA Is The Next Logical Step
Age-bin medians were useful because they exposed a clear pattern above 60. Without them, the repo would not have found the `Faecalibacterium -> E. coli` shift and would not have known where to focus pathway interpretation.

But age-bin medians are still grouped communities.

The next big question is:

Do these same growth, pathway, and provenance signals still appear when the modeling is done at individual-subject level instead of age-bin median level?

That is why subject-level FBA is the logical next step.

What subject-level work should test:
- whether the >60 signal is consistent across individuals or driven by a subset of subjects
- whether `E. coli` and `Faecalibacterium` remain the main older-subject taxa
- whether the nucleotide-interconversion and metabolite-support stories survive outside the age-bin medians

That next step is not a new biological question. It is the validation step for the story the current repo has already built.
