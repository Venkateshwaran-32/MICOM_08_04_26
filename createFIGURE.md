# createFIGURE.md

## Purpose
This file is a repo-local guide for creating, editing, and verifying presentation figures in this project.

It exists because figure work can fail in ways that are easy to miss during implementation:

- editing the wrong figure instead of the one the user asked for
- adding labels without checking the rendered result
- causing text overlaps or cramped headers
- fixing a single-species figure but forgetting the combined figure
- assuming a rerender finished before inspecting the final saved image
- prioritizing a new “better” visualization over the user’s explicitly requested edit

This file should be read before making plot or figure changes.

## Core Rule
When the user asks to change an existing figure, change that existing figure first unless they explicitly ask for a new figure.

Do not substitute a different figure type just because it is more informative.

Examples:

- If the user asks to improve `alistipes_shahii_crossfeeding_biomass_summary.png`, do that first.
- Do not switch to a mechanistic-chain figure unless the user specifically asks for a new mechanistic figure.

## Figure Workflow
Always follow this sequence:

1. Identify the exact target figure file the user cares about.
2. Identify the script that generates it.
3. Inspect the current rendered image before changing code.
4. Name the exact layout/content problems visible in the current image.
5. Edit the plotting code.
6. Rerender the figure.
7. Wait for the rerender to fully finish.
8. Inspect the final saved image again.
9. Check for overlaps, collisions, truncation, clipping, and misleading layout.
10. Only then report that the figure is done.

Do not skip step 8.

## Mandatory Visual QA
After every figure change, inspect the rendered PNG itself.

Do not rely only on reading the plotting code.
Do not assume that a layout change is correct just because it “should” work.
Do not describe a figure as fixed before visually checking the actual saved output.

For every updated figure, explicitly check:

- title does not overlap subtitle or panel headers
- panel headers do not overlap first row content
- column headers are visually distinct and have enough spacing
- long species names do not collide with values or adjacent columns
- pathway labels do not collide with bars, values, or titles
- row cards have enough vertical spacing
- right-aligned numbers have enough space from labels
- text is not clipped by figure bounds
- combined/multi-panel figures are also checked, not just single-panel outputs

## Common Mistakes To Avoid

### 1. Editing the wrong artifact
If the user asks for the summary figure, do not deliver only a different derived figure.

Bad:
- user asks to improve the summary PNG
- new mechanistic PNG is created
- original summary PNG still has readability issues

Good:
- update the requested summary PNG first
- create extra figures only if explicitly requested or clearly labeled as additional outputs

### 2. Adding labels without checking spacing
Adding column names can create a new overlap with the first data row.

If headers are added:

- give them their own header band or dedicated vertical space
- check whether they collide with the first row after rendering

### 3. Assuming one figure fix generalizes to all outputs
A single-species PNG may look fine while the combined figure is still cramped.

If a script generates:

- individual figures
- combined figures
- alternative layouts

then inspect all of them after the change.

### 4. Not accounting for different label lengths
One species may fit while another does not.

In this repo:

- `Alistipes shahii` is shorter
- `Faecalibacterium prausnitzii` is longer

A layout that works for the first may still overlap in the second.

Always test with the longest likely labels.

### 5. Looking at a stale image before the rerender finishes
If a plotting command is still running, the viewed image may still be the old version.

Always wait until the plotting command exits and confirms the saved file paths before final inspection.

### 6. Reporting “no overlaps” too early
Do not say there are no overlaps until:

- the rerender has completed
- the final image has been viewed
- long-label cases have been checked

## Layout Guidance

### Table-like left panels
If a panel is effectively a table, treat it like a table.

Use:

- fixed column anchors
- a dedicated header band
- consistent row heights
- consistent left/right alignment by column

Do not fake a table with loosely placed text elements if values must line up precisely.

Recommended columns for cross-feeding summaries:

- Producer
- Shared metabolite `[e]`
- Consumer
- Score
- Confidence

### Handling long labels
Use one of these approaches deliberately:

- shorten species names
- wrap long labels
- widen the figure/panel
- move value columns rightward

Do not leave long labels at full length if they run into numeric columns.

For row-based summary tables, a compact species label is acceptable if it preserves readability.

Example:

- `Faecalibacterium prausnitzii` -> `F. prausnitzii`

### Pathway bar panels
Pathway labels often need more space than expected.

Best practice:

- reserve a label region above each bar
- keep value labels aligned to a consistent right edge
- keep “Top reaction” on a separate line below the bar
- use conservative wrap widths

If the combined figure is cramped, increase overall figure size and panel spacing instead of squeezing text tighter.

## Content Discipline
Do not silently rename metrics or change meaning.

If the figure shows a metric such as `crossfeeding_support_score`, either:

- label it clearly, or
- intentionally shorten the label in the figure while preserving meaning

If a shortened header is used, ensure it is still understandable in context.

Examples:

- acceptable: `Score` if the panel context clearly indicates cross-feeding candidates
- better when space allows: `Support score`

## Communication Rules During Figure Work
When working on figures:

- say which exact figure you are changing
- say which exact visual issue you are fixing
- after rerendering, say which figures were checked
- if any overlap remains, say so directly and fix it

Do not claim completion before visual QA.

## Minimum Acceptance Checklist
Before closing a figure task, confirm all of the following:

- the requested figure was actually updated
- the plotting script was rerun successfully
- the final saved PNG was inspected
- no visible overlaps remain
- no key labels are ambiguous
- combined figures were inspected if generated by the same script
- any added labels still fit at presentation scale

## Project-Specific Note
For provenance figures in this repo, summary figures are often presentation-oriented and compact. That makes them especially prone to:

- header collisions
- long-species-name collisions
- cramped score/confidence columns
- pathway label crowding in the right panel

For figures under `Results/figures/micom/lysine_butyrate/`, preserve a consistent presentation footer unless the user asks otherwise:

- bottom-left interpretation box with 2 to 3 concise takeaways
- bottom-right source box naming the CSV inputs used for the figure
- enough bottom margin that both boxes remain readable at slide scale

If one figure in that lysine/butyrate set is updated for presentation clarity, check whether sibling figures generated by the same branch should keep the same footer style.

Bias toward slightly more whitespace, shorter row labels, and larger combined-figure canvases.

## If Unsure
If there is a tradeoff between density and readability, choose readability.

For slide figures, readable and slightly sparse is better than dense and technically complete but hard to parse.
