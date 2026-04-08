# AGENTS.md

## Purpose
This repository contains host-microbiome interaction modelling work centered on flux balance analysis (FBA), weighted community FBA, and a MICOM branch for cooperative-tradeoff community modelling.

Your job is to help maintain the repository documentation, especially `README.md`, without creating noisy or unnecessary edits.

## Primary Rule
Only update `README.md` when there are significant, scientifically meaningful changes to the results or FBA workflow.

Do not update `README.md` for small, cosmetic, or unrelated repository changes.

## What Counts As Significant
Treat a change as significant if one or more of the following is true:

- New FBA scenarios, branches, or analyses were added.
- Existing FBA outputs materially changed in a way that affects interpretation.
- New result files changed the story told in the README.
- Core pipeline scripts changed expected inputs, outputs, or execution order.
- New figures or summary tables replaced the current canonical outputs.
- MICOM outputs became important enough that the README should describe them differently.
- The recommended next steps, interpretation notes, or project status are no longer accurate.

## What Does Not Count As Significant
Do not update `README.md` for changes like these unless they clearly alter interpretation:

- Minor refactors or code cleanup in scripts.
- Renaming variables or internal helper functions.
- Formatting-only changes.
- Regenerated outputs that are numerically identical or effectively equivalent.
- Small file timestamp changes.
- Cache files, `.pyc` files, `.DS_Store`, or local environment changes.
- Experimental scratch outputs that are not part of the main reported analysis.

## Source Of Truth
When deciding whether the README needs updating, inspect these areas first:

- `README.md`
- `Scripts/`
- `Results/fba/`
- `Results/fba/scenarios/`
- `Results/figures/`
- `Results/micom/`
- `Results/inputs/`
- `Results/qc/`
- `Data/Processed/`

Prioritize these scripts and branches:

- Baseline COBRApy branch: `Scripts/01_load_qc_models.py` to `Scripts/07_make_figures.py`
- SG90 weighted FBA branch: `Scripts/08_prepare_abundance_inputs.py`, `Scripts/09_run_community_fba_weighted.py`, `Scripts/10_compare_weighted_vs_equal_fluxes.py`, `Scripts/11_plot_weighted_vs_equal.py`
- Beta-age branch: `Scripts/08d_prepare_beta_age_weights.py`, `Scripts/09c_run_community_fba_beta_age_weighted.py`
- All-cohort age-bin branch: `Scripts/08b_prepare_allcohort_agebin_inputs.py`, `Scripts/09b_run_community_fba_allcohort_agebins.py`, `Scripts/10b_compare_allcohort_agebins.py`, `Scripts/11b_plot_allcohort_agebins.py`, `Scripts/11c_plot_allcohort_top3_species.py`, `Scripts/11d_plot_nonzero_and_log_biomass.py`
- MICOM branch: `Scripts/13_prepare_micom_inputs.py` to `Scripts/16_find_differential_pathways.py`
- Validation step: `Scripts/99_validate_pipeline_outputs.py`

## README Update Decision Process
Before editing `README.md`, do this:

1. Compare current repository state against what `README.md` currently claims.
2. Check whether new or changed outputs in `Results/` materially alter:
   - project status
   - main findings
   - active branches
   - canonical output files
   - interpretation notes
   - recommended next steps
3. Check whether script changes affect:
   - pipeline order
   - scenario names
   - file names
   - output locations
   - required commands
4. If changes are minor or not scientifically meaningful, do not edit `README.md`.

## How To Edit README
If an update is justified:

- Keep the current structure and tone unless it is clearly outdated.
- Make the smallest edit that restores accuracy.
- Preserve the distinction between:
  - SG90 outputs using `sg90_*`
  - all-cohort outputs using `allcohort_agebin_*`
  - beta-age outputs using `beta_age_*`
  - MICOM outputs in `Results/micom/`
- Prefer precise file and script references over vague prose.
- Do not invent results, interpretations, or completed analyses that are not visible in the repository.
- Do not overstate conclusions from sparse or exploratory outputs.

## Scientific Caution
Be conservative when describing results.

- Weighted community FBA may produce sparse optima; do not describe zero-growth taxa as a pipeline failure.
- Distinguish exploratory outputs from canonical results.
- Do not treat the presence of a file as proof of a robust biological conclusion without supporting context from the pipeline and summaries.
- If outputs conflict, prefer summaries and validated artifacts over ad hoc intermediate files.

## Ignore Noise
Usually ignore these when deciding whether to update the README:

- `.venv/`
- `.cache/`
- `.matplotlib/`
- `Library/`
- `Scripts/__pycache__/`
- `.DS_Store`

Also ignore accidental or malformed files unless they clearly matter to the documented workflow.

## Output Behavior
When acting as a maintenance agent:

- If there is no meaningful results or FBA change, leave `README.md` untouched.
- If you update `README.md`, include a short summary of:
  - what changed
  - why it mattered
  - which files or scripts justified the update

## Style Preferences
- Prefer concise, high-signal edits.
- Preserve existing terminology such as FBA, community FBA, SG90, all-cohort age-bin, and MICOM.
- Avoid rewriting large sections unless the current README is materially wrong.
