from pathlib import Path
import re
from typing import Optional, Set
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# Script 15: plot lysine-focused pathway fluxes from MICOM outputs
# Inputs:
# - Results/micom/pathway_flux/proper_age_bins/reaction_fluxes_long_by_agegroup_diet.csv
# Outputs:
# - Results/micom/lysine_butyrate/lysine_related_fluxes_long.csv
# - Results/micom/lysine_butyrate/lysine_biosynthesis_flux_summary.csv
# - Results/micom/lysine_butyrate/lysine_to_butyrate_flux_summary.csv
# - Results/micom/lysine_butyrate/lysine_to_butyrate_reaction_step_summary.csv
# - Results/figures/micom/lysine_butyrate/lysine_biosynthesis_flux_by_agegroup_diet.png
# - Results/figures/micom/lysine_butyrate/lysine_to_butyrate_flux_by_agegroup_diet.png
# - Results/figures/micom/lysine_butyrate/lysine_to_butyrate_reaction_steps_by_agegroup_diet.png
# Run:
# - .venv\Scripts\python Scripts\15_plot_lysine_pathways.py
# Expected runtime:
# - ~3 to 20 seconds
# -------------------------------------------------------------------

# This script extracts curated lysine-related rows from MICOM reaction flux outputs
# and generates summary tables + figures for:
# - L-lysine biosynthesis / diaminopimelate-linked lysine steps
# - L-lysine to butyrate candidate routes

PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_FLUX = (
    PROJECT_ROOT
    / "Results"
    / "micom"
    / "pathway_flux"
    / "proper_age_bins"
    / "reaction_fluxes_long_by_agegroup_diet.csv"
)

OUT_DIR = PROJECT_ROOT / "Results" / "micom" / "lysine_butyrate"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "lysine_butyrate"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_LYS = OUT_DIR / "lysine_related_fluxes_long.csv"
OUT_LYS_BIO = OUT_DIR / "lysine_biosynthesis_flux_summary.csv"
OUT_LYS_BUT = OUT_DIR / "lysine_to_butyrate_flux_summary.csv"
OUT_LYS_BUT_STEPS = OUT_DIR / "lysine_to_butyrate_reaction_step_summary.csv"

FIG_LYS_BIO = FIG_DIR / "lysine_biosynthesis_flux_by_agegroup_diet.png"
FIG_LYS_BUT = FIG_DIR / "lysine_to_butyrate_flux_by_agegroup_diet.png"
FIG_LYS_BUT_STEPS = FIG_DIR / "lysine_to_butyrate_reaction_steps_by_agegroup_diet.png"

FLUX_EPS = 1e-9

# Curated route-specific identifiers and names are preferable to broad "lys" text
# matching because generic fragments such as "lys" can hit unrelated words
# (for example "glycolysis").
CURATED_RULES = {
    "lysine_biosynthesis": {
        "reaction_ids": {
            "ACCOALYS2",
            "26DAPLLAT",
            "DAPDC",
            "DAPE",
            "DAPAT",
            "DAPDA",
            "DAPMDH",
            "SDPDS",
            "SDPTA",
            "PAPPT3",
            "UAAGDS",
        },
        "text_patterns": [
            r"\blysine biosynthesis\b",
            r"\bdiaminopimelate decarboxylase\b",
            r"\bdiaminopimelate epimerase\b",
            r"\bdiaminopimelate aminotransferase\b",
            r"\bmeso-2,6-diaminopimelate\b",
            r"\bdapdc\b",
            r"\bdape\b",
            r"\b26dapllat\b",
        ],
    },
    "lysine_to_butyrate_candidate": {
        "reaction_ids": {
            "LYSAM",
            "3ABUTCOAL",
            "BTCOADH",
            "PBUTT",
            "BUTKr",
            "BTCOAACCOAT",
            "BUTCTr",
        },
        "text_patterns": [
            r"\blysine 2,3-aminomutase\b",
            r"\b3-aminobutyryl-coa ammonia-lyase\b",
            r"\bcrotonoyl-coa\b",
            r"\bbutyryl-coa\b",
            r"\bbutyrate kinase\b",
            r"\bphosphate butyryltransferase\b",
            r"\bbutyryl-coa:acetate coa-transferase\b",
        ],
    },
}

REACTION_STEP_LABELS = {
    "LYSAM": "Lysine metabolism (LYSAM: lysine 2,3-aminomutase)",
    "3ABUTCOAL": "Lysine metabolism (3ABUTCOAL: 3-aminobutyryl-CoA ammonia-lyase)",
    "BTCOADH": "Butanoate metabolism (BTCOADH: crotonoyl-CoA to butyryl-CoA)",
    "PBUTT": "Butanoate metabolism (PBUTT: phosphate butyryltransferase)",
    "BUTKr": "Butanoate metabolism (BUTKr: butyrate kinase)",
    "BTCOAACCOAT": "Butanoate metabolism (BTCOAACCOAT: butyryl-CoA:acetate CoA-transferase)",
    "BUTCTr": "Butanoate metabolism (BUTCTr: acetyl-CoA:butyrate-CoA transferase)",
}

SUMMARY_BOX = {
    "biosynthesis": (
        "Summary\n"
        "1. Lysine metabolism remains high across all age groups.\n"
        "2. The strongest lysine biosynthesis signal appears in older bins,\n"
        "   especially 71_80 high_fiber.\n"
        "3. This figure shows pathway activity magnitude, not net export."
    ),
    "candidate": (
        "Summary\n"
        "1. Lysine-to-butyrate candidate flux is near zero in younger bins.\n"
        "2. Nonzero butanoate signal appears mainly in 71_80 and 81_plus.\n"
        "3. The candidate branch is much smaller than total lysine metabolism."
    ),
    "steps": (
        "Summary\n"
        "1. Only downstream butanoate steps are nonzero in these MICOM runs.\n"
        "2. BTCOADH and BTCOAACCOAT carry the observed older-bin signal.\n"
        "3. Upstream lysine-specific candidate steps remain effectively zero."
    ),
}

SOURCE_BOX = {
    "biosynthesis": (
        "Sources\n"
        "reaction_fluxes_long_\n"
        "by_agegroup_diet.csv\n"
        "lysine_biosynthesis_flux_summary.csv"
    ),
    "candidate": (
        "Sources\n"
        "reaction_fluxes_long_\n"
        "by_agegroup_diet.csv\n"
        "lysine_to_butyrate_flux_summary.csv"
    ),
    "steps": (
        "Sources\n"
        "reaction_fluxes_long_\n"
        "by_agegroup_diet.csv\n"
        "lysine_to_butyrate_reaction_step_summary.csv"
    ),
}


def text_join(row) -> str:
    # Merge searchable reaction text from multiple columns.
    vals = [row.get("reaction_id", ""), row.get("reaction_name", ""), row.get("pathway", "")]
    return " | ".join([str(v) for v in vals if pd.notna(v)]).lower()


def curated_match_mask(df: pd.DataFrame, reaction_ids: set[str], text_patterns: list[str]) -> pd.Series:
    rid = df["reaction_id"].fillna("").astype(str).str.upper()
    id_match = rid.isin(reaction_ids)
    if not text_patterns:
        return id_match

    text_match = pd.Series(False, index=df.index)
    for pat in text_patterns:
        text_match = text_match | df["text"].str.contains(pat, regex=True, na=False)
    return id_match | text_match


def classify_flux_rows(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["text"] = d.apply(text_join, axis=1)

    d["is_lysine_related"] = False
    d["is_biosyn_hint"] = False
    d["is_butyrate_hint"] = False
    d["pathway_class"] = "other"

    bio_mask = curated_match_mask(
        d,
        reaction_ids=CURATED_RULES["lysine_biosynthesis"]["reaction_ids"],
        text_patterns=CURATED_RULES["lysine_biosynthesis"]["text_patterns"],
    )
    but_mask = curated_match_mask(
        d,
        reaction_ids=CURATED_RULES["lysine_to_butyrate_candidate"]["reaction_ids"],
        text_patterns=CURATED_RULES["lysine_to_butyrate_candidate"]["text_patterns"],
    )

    d.loc[bio_mask, "is_lysine_related"] = True
    d.loc[bio_mask, "is_biosyn_hint"] = True
    d.loc[bio_mask, "pathway_class"] = "lysine_biosynthesis"

    d.loc[but_mask, "is_lysine_related"] = True
    d.loc[but_mask, "is_butyrate_hint"] = True
    d.loc[but_mask, "pathway_class"] = "lysine_to_butyrate_candidate"

    return d


def summarize(df: pd.DataFrame, klass: str, include_zero_flux: bool = False) -> pd.DataFrame:
    # Aggregate pathway signal for each age-group x diet block.
    sub = df[df["pathway_class"] == klass].copy()
    if not include_zero_flux:
        sub = sub[sub["abs_flux"] > FLUX_EPS].copy()
    if sub.empty:
        return pd.DataFrame(columns=["age_group", "diet", "pathway", "sum_abs_flux", "n_reactions"])

    out = (
        sub.groupby(["age_group", "diet", "pathway"], as_index=False)
        .agg(sum_abs_flux=("abs_flux", "sum"), n_reactions=("reaction_id", "count"))
        .sort_values(["age_group", "diet", "sum_abs_flux"], ascending=[True, True, False])
    )
    return out


def summarize_reaction_steps(df: pd.DataFrame, include_zero_flux: bool = True) -> pd.DataFrame:
    sub = df[df["pathway_class"] == "lysine_to_butyrate_candidate"].copy()
    sub["reaction_id_upper"] = sub["reaction_id"].fillna("").astype(str).str.upper()
    sub = sub[sub["reaction_id_upper"].isin(REACTION_STEP_LABELS)].copy()
    if sub.empty:
        return pd.DataFrame(
            columns=["age_group", "diet", "reaction_step", "sum_abs_flux", "n_rows"]
        )

    sub["reaction_step"] = sub["reaction_id_upper"].map(REACTION_STEP_LABELS)
    out = (
        sub.groupby(["age_group", "diet", "reaction_step"], as_index=False)
        .agg(sum_abs_flux=("abs_flux", "sum"), n_rows=("reaction_id", "count"))
        .sort_values(["age_group", "diet", "reaction_step"], ascending=[True, True, True])
    )
    if not include_zero_flux:
        active_steps = out.loc[out["sum_abs_flux"] > FLUX_EPS, "reaction_step"].unique().tolist()
        if not active_steps:
            return pd.DataFrame(
                columns=["age_group", "diet", "reaction_step", "sum_abs_flux", "n_rows"]
            )
        out = out[out["reaction_step"].isin(active_steps)].copy()
    return out


def format_flux_label(value: float) -> str:
    if abs(value) <= FLUX_EPS:
        return "0"
    if value >= 1000:
        return f"{value/1000:.1f}k"
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


def add_bar_labels(ax, min_fraction_of_max: float = 0.06, always_label_index: Optional[Set[str]] = None):
    heights = [bar.get_height() for container in ax.containers for bar in container]
    positive = [h for h in heights if h > FLUX_EPS]
    if not positive:
        return

    max_height = max(positive)
    min_height = max_height * min_fraction_of_max
    always_label_index = always_label_index or set()

    for container in ax.containers:
        labels = []
        for bar in container:
            h = bar.get_height()
            x_center = bar.get_x() + (bar.get_width() / 2)
            idx = None
            for tick, label in zip(ax.get_xticks(), ax.get_xticklabels()):
                if abs(x_center - tick) < 0.5:
                    idx = label.get_text()
                    break
            force_label = idx in always_label_index if idx is not None else False
            if not force_label and h <= max(FLUX_EPS, min_height):
                labels.append("")
            else:
                labels.append(format_flux_label(h))
        ax.bar_label(container, labels=labels, padding=2, fontsize=7)

    ax.set_ylim(top=max_height * 1.18)


def add_footer_boxes(fig, summary_text: str, source_text: str):
    fig.text(
        0.06,
        0.055,
        summary_text,
        fontsize=8.5,
        color="#173042",
        va="top",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#fff8e8", "edgecolor": "#d8c9a8", "linewidth": 0.8},
    )
    fig.text(
        0.94,
        0.055,
        source_text,
        fontsize=8,
        color="#5c7280",
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#f7f4ec", "edgecolor": "#d8d2c3", "linewidth": 0.8},
    )


def plot_grouped_bars(df: pd.DataFrame, title: str, out_path: Path):
    if df.empty:
        print(f"Skipping plot (no rows): {out_path}")
        return

    top_pathways = (
        # Keep only top pathways to make plots readable.
        df.groupby("pathway", as_index=False)["sum_abs_flux"]
        .sum()
        .sort_values("sum_abs_flux", ascending=False)
        .head(12)["pathway"]
        .tolist()
    )
    p = (
        df[df["pathway"].isin(top_pathways)]
        .groupby(["pathway", "age_group", "diet"], as_index=False)["sum_abs_flux"]
        .sum()
    )
    p["age_diet"] = p["age_group"].astype(str) + "|" + p["diet"].astype(str)

    pivot = p.pivot(index="pathway", columns="age_diet", values="sum_abs_flux").fillna(0.0)
    if pivot.empty:
        print(f"Skipping plot (empty pivot): {out_path}")
        return

    ax = pivot.plot(kind="bar", figsize=(13, 8.2))
    ax.set_title(title)
    ax.set_ylabel("sum abs flux")
    ax.set_xlabel("pathway")
    ax.tick_params(axis="x", rotation=40)
    always_label_index = None
    if "L-lysine Biosynthesis" in title:
        always_label_index = {"Cell wall biosynthesis"}
        summary_key = "biosynthesis"
    elif "L-lysine to Butyrate Candidate Pathway Flux" in title:
        always_label_index = {"Lysine metabolism"}
        summary_key = "candidate"
    else:
        summary_key = "biosynthesis"
    add_bar_labels(ax, min_fraction_of_max=0.08, always_label_index=always_label_index)
    fig = ax.get_figure()
    add_footer_boxes(fig, SUMMARY_BOX[summary_key], SOURCE_BOX[summary_key])
    plt.tight_layout(rect=[0.02, 0.12, 0.98, 0.98])
    plt.savefig(out_path, dpi=220)
    plt.close()
    print(f"Saved: {out_path}")


def plot_reaction_steps(df: pd.DataFrame, title: str, out_path: Path):
    if df.empty:
        print(f"Skipping plot (no rows): {out_path}")
        return

    order = [label for _, label in REACTION_STEP_LABELS.items()]
    p = df.copy()
    p["reaction_step"] = pd.Categorical(p["reaction_step"], categories=order, ordered=True)
    p = p.sort_values("reaction_step")
    p["age_diet"] = p["age_group"].astype(str) + "|" + p["diet"].astype(str)

    pivot = p.pivot(index="reaction_step", columns="age_diet", values="sum_abs_flux").fillna(0.0)
    if pivot.empty:
        print(f"Skipping plot (empty pivot): {out_path}")
        return

    ax = pivot.plot(kind="bar", figsize=(15, 9))
    ax.set_title(title)
    ax.set_ylabel("sum abs flux")
    ax.set_xlabel("curated reaction step")
    ax.tick_params(axis="x", rotation=0)
    add_bar_labels(ax, min_fraction_of_max=0.04)
    fig = ax.get_figure()
    add_footer_boxes(fig, SUMMARY_BOX["steps"], SOURCE_BOX["steps"])
    plt.tight_layout(rect=[0.02, 0.12, 0.98, 0.98])
    plt.savefig(out_path, dpi=220)
    plt.close()
    print(f"Saved: {out_path}")


def main():
    if not IN_FLUX.exists():
        raise FileNotFoundError(f"Missing input: {IN_FLUX}. Run Script 14 first.")

    df = pd.read_csv(IN_FLUX)
    required = {"age_group", "diet", "reaction_id", "reaction_name", "pathway", "abs_flux"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input flux table must contain columns: {sorted(required)}")

    c = classify_flux_rows(df)
    # Keep a detailed table for manual audit of curated matching.
    lys_only = c[c["is_lysine_related"] | (c["pathway_class"] == "lysine_to_butyrate_candidate")].copy()
    lys_only.to_csv(OUT_LYS, index=False)

    lys_bio = summarize(c, "lysine_biosynthesis", include_zero_flux=True)
    lys_but = summarize(c, "lysine_to_butyrate_candidate", include_zero_flux=True)
    lys_but_steps = summarize_reaction_steps(c, include_zero_flux=False)

    lys_bio.to_csv(OUT_LYS_BIO, index=False)
    lys_but.to_csv(OUT_LYS_BUT, index=False)
    lys_but_steps.to_csv(OUT_LYS_BUT_STEPS, index=False)

    plot_grouped_bars(
        lys_bio,
        "L-lysine Biosynthesis Pathway Flux by Age Group and Diet",
        FIG_LYS_BIO,
    )
    plot_grouped_bars(
        lys_but,
        "L-lysine to Butyrate Candidate Pathway Flux by Age Group and Diet",
        FIG_LYS_BUT,
    )
    plot_reaction_steps(
        lys_but_steps,
        "L-lysine to Butyrate Curated Reaction Steps by Age Group and Diet",
        FIG_LYS_BUT_STEPS,
    )

    print(f"Saved: {OUT_LYS}")
    print(f"Saved: {OUT_LYS_BIO}")
    print(f"Saved: {OUT_LYS_BUT}")
    print(f"Saved: {OUT_LYS_BUT_STEPS}")


if __name__ == "__main__":
    main()
