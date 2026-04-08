from pathlib import Path
import textwrap
import re
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROV_DIR = PROJECT_ROOT / "Results" / "fba" / "provenance"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "provenance"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SPECIES_CONFIG = [
    {
        "species_slug": "Alistipes_shahii_WAL_8301_AGORA1_03",
        "title": "Alistipes shahii",
        "diet": "high_fiber",
        "clean_csv": PROV_DIR / "high_fiber__Alistipes_shahii_WAL_8301_AGORA1_03__clean_crossfeeding_biomass_summary.csv",
        "pathway_csv": PROV_DIR / "high_fiber__Alistipes_shahii_WAL_8301_AGORA1_03__pathway_overview.csv",
        "out_png": FIG_DIR / "alistipes_shahii_crossfeeding_biomass_summary.png",
        "accent": "#1f5a89",
    },
    {
        "species_slug": "Faecalibacterium_prausnitzii_M21_2__AGORA1_03",
        "title": "Faecalibacterium prausnitzii",
        "diet": "western",
        "clean_csv": PROV_DIR / "western__Faecalibacterium_prausnitzii_M21_2__AGORA1_03__clean_crossfeeding_biomass_summary.csv",
        "pathway_csv": PROV_DIR / "western__Faecalibacterium_prausnitzii_M21_2__AGORA1_03__pathway_overview.csv",
        "out_png": FIG_DIR / "faecalibacterium_prausnitzii_crossfeeding_biomass_summary.png",
        "accent": "#7a4b1f",
    },
]

COMBINED_OUT = FIG_DIR / "crossfeeding_biomass_summary_two_species.png"
ALISTIPES_CHAIN_OUT = FIG_DIR / "alistipes_shahii_adenosine_mechanistic_chain.png"

BG = "#fcfbf7"
TEXT = "#173042"
MUTED = "#597182"
PANEL = "#f4efe6"
HIGH = "#1f8f6a"
MEDIUM = "#c5841f"
LOW = "#a85c59"


def short_species_name(value: str) -> str:
    text = str(value)
    parts = text.split("_")
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return text


def species_display_name(value: str) -> str:
    short = short_species_name(value)
    return short.replace("_", " ")


def compact_species_label(value: str) -> str:
    text = species_display_name(value)
    parts = text.split()
    if len(parts) >= 2 and len(text) > 22:
        return f"{parts[0][0]}. {parts[1]}"
    return text


def short_metabolite_label(met_id: str, met_name: str) -> str:
    name = str(met_name).strip()
    if name and name.lower() != "nan":
        return name
    return str(met_id)


def wrap(text: str, width: int) -> str:
    return "\n".join(textwrap.wrap(str(text), width=width)) if text else ""


def wrap_pathway_label(text: str, width: int) -> str:
    value = str(text).replace("/", "/ ")
    return "\n".join(textwrap.wrap(value, width=width)) if value else ""


def sanitize_for_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(value)).strip("_")


def confidence_color(level: str) -> str:
    if level == "high":
        return HIGH
    if level == "medium":
        return MEDIUM
    return LOW


def load_inputs(config: dict):
    clean_df = pd.read_csv(config["clean_csv"])
    pathway_df = pd.read_csv(config["pathway_csv"])
    clean_df = clean_df.sort_values(
        ["confidence", "crossfeeding_support_score", "biomass_association_score"],
        ascending=[True, False, False],
    ).copy()
    pathway_df = pathway_df.sort_values("biomass_association_score", ascending=False).copy()
    return clean_df, pathway_df


def draw_rounded_box(ax, x, y, w, h, fc, ec=None, lw=1.2, pad=0.012, rounding=0.015):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad={pad},rounding_size={rounding}",
        facecolor=fc,
        edgecolor=ec if ec else fc,
        linewidth=lw,
    )
    ax.add_patch(patch)
    return patch


def fit_text(text: str, max_len: int) -> str:
    value = str(text)
    if len(value) <= max_len:
        return value
    return value[: max_len - 1] + "…"


def wrap_id(value: str, width: int) -> str:
    text = str(value)
    chunks = textwrap.wrap(text, width=width, break_long_words=True, break_on_hyphens=False)
    return "\n".join(chunks)


def clip_with_ellipsis(text: str, max_len: int) -> str:
    value = str(text)
    if len(value) <= max_len:
        return value
    return value[: max_len - 3] + "..."


def display_source_name(path: Path) -> str:
    name = path.name
    for suffix in (
        "__clean_crossfeeding_biomass_summary.csv",
        "__pathway_overview.csv",
        "__provenance_edges.csv",
        "__reaction_steps.csv",
    ):
        if suffix in name:
            return suffix.replace("__", "", 1)
    return name


def source_box_text(derived_paths: list[Path], upstream_paths: Optional[list[Path]] = None, width: int = 58) -> str:
    lines = ["Sources", "Derived"]
    lines.extend(f"- {display_source_name(fp)}" for fp in derived_paths)
    if upstream_paths:
        lines.append("Upstream inputs")
        lines.extend(f"- {display_source_name(fp)}" for fp in upstream_paths)
    return "\n".join(wrap(line, width) for line in lines)


def load_alistipes_chain_inputs():
    config = SPECIES_CONFIG[0]
    clean_df, pathway_df = load_inputs(config)
    top = clean_df.iloc[0]

    met_slug = sanitize_for_filename(top["shared_metabolite_id"])
    edges_csv = PROV_DIR / f"{config['diet']}__{config['species_slug']}__{met_slug}__provenance_edges.csv"
    steps_csv = PROV_DIR / f"{config['diet']}__{config['species_slug']}__Nucleotide_interconversion__reaction_steps.csv"
    edges = pd.read_csv(edges_csv)
    steps = pd.read_csv(steps_csv)

    top_producer = (
        edges[edges["role"] == "shared_pool_secretion"]
        .sort_values("abs_flux", ascending=False)
        .iloc[0]
    )
    uptake = (
        edges[(edges["role"] == "shared_pool_uptake") & (edges["species_to"] == config["species_slug"])]
        .sort_values("abs_flux", ascending=False)
        .iloc[0]
    )
    internal = (
        edges[(edges["role"] == "internal_processing") & (edges["species_to"] == config["species_slug"])]
        .sort_values("abs_flux", ascending=False)
    )
    selected_ids = ["PUNP1", "ADK1", "NDPK3", "ADNK1"]
    selected_steps = steps[steps["reaction_id"].isin(selected_ids)].copy()
    selected_steps["sort_key"] = selected_steps["reaction_id"].map({rid: i for i, rid in enumerate(selected_ids)})
    selected_steps = selected_steps.sort_values(["sort_key", "step"]).drop(columns=["sort_key"])

    pathway_row = pathway_df[pathway_df["pathway"] == top["consumer_top_pathway"]].iloc[0]
    other_producers = (
        edges[edges["role"] == "shared_pool_secretion"]["species_from"]
        .dropna()
        .astype(str)
        .nunique()
        - 1
    )
    other_consumers = (
        edges[edges["role"] == "shared_pool_uptake"]["species_to"]
        .dropna()
        .astype(str)
        .nunique()
        - 1
    )

    return {
        "config": config,
        "top": top,
        "top_producer": top_producer,
        "uptake": uptake,
        "internal": internal,
        "steps": selected_steps,
        "pathway_row": pathway_row,
        "other_producers": max(int(other_producers), 0),
        "other_consumers": max(int(other_consumers), 0),
        "edges_csv": edges_csv,
        "steps_csv": steps_csv,
    }


def draw_species_panel(ax, config: dict):
    clean_df, pathway_df = load_inputs(config)
    focus_rows = clean_df.head(5).copy()
    pathway_rows = pathway_df.head(4).copy()

    ax.set_facecolor(BG)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    draw_rounded_box(ax, 0.02, 0.03, 0.96, 0.94, PANEL, ec="#dfd4c3", lw=1.4)
    ax.text(0.05, 0.94, config["title"], fontsize=18, weight="bold", color=TEXT, va="top")
    ax.text(0.05, 0.90, f"Community FBA provenance summary | diet = {config['diet']}", fontsize=10.5, color=MUTED, va="top")

    left_x = 0.05
    left_w = 0.58
    right_x = 0.68
    right_w = 0.27

    ax.text(left_x, 0.83, "Cross-feeding candidates retained", fontsize=12.5, weight="bold", color=TEXT)
    draw_rounded_box(ax, left_x, 0.774, left_w, 0.034, fc="#efe8dc", ec="#e3d8c7", lw=0.8)
    col_prod = left_x + 0.015
    col_arrow1 = left_x + 0.165
    col_met = left_x + 0.19
    col_arrow2 = left_x + 0.345
    col_cons = left_x + 0.37
    col_score = left_x + 0.505
    col_conf = left_x + 0.575
    header_y = 0.790
    ax.text(col_prod, header_y, "Producer", fontsize=8.0, color=MUTED, va="center", ha="left", weight="bold")
    ax.text(col_met, header_y, "Shared metabolite [e]", fontsize=8.0, color=MUTED, va="center", ha="left", weight="bold")
    ax.text(col_cons, header_y, "Consumer", fontsize=8.0, color=MUTED, va="center", ha="left", weight="bold")
    ax.text(col_score, header_y, "Score", fontsize=8.0, color=MUTED, va="center", ha="right", weight="bold")
    ax.text(col_conf, header_y, "Confidence", fontsize=8.0, color=MUTED, va="center", ha="right", weight="bold")

    y = 0.736
    step = 0.116
    for row in focus_rows.itertuples(index=False):
        color = confidence_color(row.confidence)
        draw_rounded_box(ax, left_x, y - 0.048, left_w, 0.076, fc="white", ec="#d7d1c7", lw=1.0)
        ax.text(col_prod, y, wrap(short_species_name(row.producer_species), 21), fontsize=8.8, color=TEXT, va="center", ha="left")
        ax.text(col_arrow1, y, "\u2192", fontsize=10.5, color=MUTED, va="center", ha="center")
        ax.text(col_met, y, wrap(short_metabolite_label(row.shared_metabolite_id, row.shared_metabolite_name), 18), fontsize=8.8, color=TEXT, va="center", ha="left")
        ax.text(col_arrow2, y, "\u2192", fontsize=10.5, color=MUTED, va="center", ha="center")
        ax.text(col_cons, y, compact_species_label(config["species_slug"]), fontsize=8.8, color=TEXT, va="center", ha="left")
        ax.text(col_score, y, f"{row.crossfeeding_support_score:.1f}", fontsize=8.9, color=color, va="center", ha="right", weight="bold")
        ax.text(col_conf, y, row.confidence.upper(), fontsize=8.0, color=color, va="center", ha="right", weight="bold")
        y -= step

    ax.text(right_x, 0.83, "Top biomass-adjacent pathways", fontsize=12.5, weight="bold", color=TEXT)
    y = 0.765
    bar_max = pathway_rows["biomass_association_score"].max() if not pathway_rows.empty else 1.0
    for row in pathway_rows.itertuples(index=False):
        ax.text(right_x, y + 0.032, wrap_pathway_label(row.pathway, 17), fontsize=8.4, color=TEXT, va="bottom", ha="left")
        draw_rounded_box(ax, right_x, y - 0.002, right_w - 0.03, 0.015, fc="#ece5d9", ec="#ece5d9", lw=0.5)
        width = (right_w - 0.03) * (row.biomass_association_score / bar_max if bar_max else 0)
        draw_rounded_box(ax, right_x, y - 0.002, width, 0.015, fc=config["accent"], ec=config["accent"], lw=0.5)
        ax.text(right_x + right_w - 0.01, y + 0.005, f"{row.biomass_association_score:.2f}", fontsize=8.0, color=TEXT, ha="right", va="center")
        ax.text(right_x, y - 0.024, f"Top reaction: {row.top_reaction_id}", fontsize=7.4, color=MUTED, va="center", ha="left")
        y -= 0.135

    ax.text(0.05, 0.190, "Score meaning", fontsize=12.5, weight="bold", color=TEXT)
    score_note = (
        "In this table, Score = biomass_association_score for the consumer top pathway. "
        "Formula: biomass_association_score = pathway_flux_share x consumer_growth, "
        "where pathway_flux_share = consumer_top_pathway_abs_flux / consumer_total_internal_abs_flux."
    )
    ax.text(0.05, 0.180, wrap(score_note, 92), fontsize=8.6, color=MUTED, va="top", ha="left")

    ax.text(0.05, 0.085, "Interpretation", fontsize=12.5, weight="bold", color=TEXT)
    if not clean_df.empty:
        top = clean_df.iloc[0]
        interp = (
            f"Strongest cross-feeding signal: {short_species_name(top['producer_species'])} supplies "
            f"{short_metabolite_label(top['shared_metabolite_id'], top['shared_metabolite_name'])}. "
            f"Growth-relevant pathways are dominated by central metabolism, with "
            f"{pathway_rows.iloc[min(3, len(pathway_rows)-1)]['pathway'] if not pathway_rows.empty else 'internal processing'} "
            f"appearing as an additional support layer."
        )
    else:
        interp = "No retained cross-feeding candidates passed the strict summary filter."
    ax.text(0.05, 0.0750, wrap(interp, 88), fontsize=8.4, color=TEXT, va="top", ha="left")
    upstream_paths = [
        PROV_DIR.parent / "community_species_connector_fluxes_by_diet.csv",
        PROV_DIR.parent / "community_biomass_associated_pathways_by_diet.csv",
        PROV_DIR.parent / "community_reaction_fluxes_by_diet.csv",
    ]
    ax.text(
        0.94,
        0.13,
        source_box_text([config["clean_csv"], config["pathway_csv"]], upstream_paths=upstream_paths),
        fontsize=8.0,
        color=MUTED,
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#f7f4ec", "edgecolor": "#d8d2c3", "linewidth": 0.8},
    )


def draw_metric_pill(ax, x, y, w, h, title, value, fill):
    draw_rounded_box(ax, x, y, w, h, fc=fill, ec=fill, lw=0.8)
    ax.text(x + 0.012, y + h * 0.62, title, fontsize=8.3, color=MUTED, va="center", ha="left")
    ax.text(x + 0.012, y + h * 0.28, value, fontsize=10.0, color=TEXT, va="center", ha="left", weight="bold")


def draw_alistipes_mechanistic_chain(ax):
    data = load_alistipes_chain_inputs()
    config = data["config"]
    top = data["top"]
    top_producer = data["top_producer"]
    uptake = data["uptake"]
    pathway_row = data["pathway_row"]
    steps = data["steps"]

    ax.set_facecolor(BG)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    draw_rounded_box(ax, 0.02, 0.03, 0.96, 0.94, PANEL, ec="#dfd4c3", lw=1.4)
    ax.text(0.05, 0.95, "Alistipes shahii: adenosine provenance chain", fontsize=18, weight="bold", color=TEXT, va="top")
    ax.text(
        0.05,
        0.915,
        "Original COBRApy community FBA | high_fiber | one highlighted route from shared metabolite uptake to biomass-supporting pathway context",
        fontsize=10.2,
        color=MUTED,
        va="top",
    )

    ax.text(0.05, 0.865, "Mechanistic chain", fontsize=13.0, weight="bold", color=TEXT, va="top")

    x_positions = [0.05, 0.245, 0.435, 0.64, 0.855]
    widths = [0.145, 0.135, 0.155, 0.165, 0.11]
    y_box = 0.64
    h_box = 0.16
    titles = [
        "Producer species",
        "Shared metabolite [e]",
        "Consumer uptake [e -> c]",
        "Key internal reactions [c]",
        "Biomass-support context",
    ]
    fills = ["#ffffff", "#f8f4ec", "#ffffff", "#ffffff", "#f8f4ec"]
    for x, w, title, fill in zip(x_positions, widths, titles, fills):
        draw_rounded_box(ax, x, y_box, w, h_box, fc=fill, ec="#d7d1c7", lw=1.0, pad=0.006, rounding=0.014)
        ax.text(x + 0.012, y_box + h_box - 0.02, title, fontsize=9.2, weight="bold", color=TEXT, va="top", ha="left")

    ax.text(
        x_positions[0] + 0.012,
        y_box + 0.092,
        wrap(species_display_name(top_producer["species_from"]), 18),
        fontsize=10.0,
        color=TEXT,
        va="center",
        ha="left",
    )
    ax.text(
        x_positions[0] + 0.012,
        y_box + 0.04,
        f"Top producer flux\n{float(top_producer['abs_flux']):.2f}",
        fontsize=8.8,
        color=MUTED,
        va="center",
        ha="left",
    )

    ax.text(
        x_positions[1] + 0.012,
        y_box + 0.1,
        fit_text(f"{top['shared_metabolite_id']}  {top['shared_metabolite_name']}", 23),
        fontsize=10.0,
        color=TEXT,
        va="center",
        ha="left",
    )
    ax.text(
        x_positions[1] + 0.012,
        y_box + 0.04,
        "Extracellular shared-pool metabolite",
        fontsize=8.5,
        color=MUTED,
        va="center",
        ha="left",
    )

    ax.text(
        x_positions[2] + 0.012,
        y_box + 0.1,
        wrap(config["title"], 18),
        fontsize=10.0,
        color=TEXT,
        va="center",
        ha="left",
    )
    ax.text(
        x_positions[2] + 0.012,
        y_box + 0.045,
        f"{wrap_id(uptake['reaction_id'], 18)}\nADNCNT3tc transport\nto adn[c]",
        fontsize=7.0,
        color=MUTED,
        va="center",
        ha="left",
    )

    reaction_lines = []
    for row in steps.itertuples(index=False):
        reaction_lines.append(f"{row.reaction_id}  {row.flux:.2f}")
        reaction_lines.append(clip_with_ellipsis(row.equation_in_flux_direction, 29))
    ax.text(
        x_positions[3] + 0.012,
        y_box + 0.102,
        "\n".join(reaction_lines),
        fontsize=7.0,
        color=TEXT,
        va="top",
        ha="left",
    )

    ax.text(
        x_positions[4] + 0.012,
        y_box + 0.102,
        wrap(pathway_row["pathway"], 14),
        fontsize=9.4,
        color=TEXT,
        va="center",
        ha="left",
        weight="bold",
    )
    ax.text(
        x_positions[4] + 0.012,
        y_box + 0.04,
        wrap("Supports nucleotide precursor supply for growth", 20),
        fontsize=7.8,
        color=MUTED,
        va="center",
        ha="left",
    )

    for idx in range(4):
        x0 = x_positions[idx] + widths[idx]
        x1 = x_positions[idx + 1]
        y_mid = y_box + h_box * 0.5
        ax.annotate(
            "",
            xy=(x1 - 0.01, y_mid),
            xytext=(x0 + 0.01, y_mid),
            arrowprops=dict(arrowstyle="->", lw=1.8, color=config["accent"]),
        )

    ax.text(0.05, 0.56, "Labeled flux metrics", fontsize=12.5, weight="bold", color=TEXT, va="top")
    draw_metric_pill(ax, 0.05, 0.46, 0.18, 0.08, "Producer flux", f"{float(top['producer_flux']):.2f}", "#ffffff")
    draw_metric_pill(ax, 0.25, 0.46, 0.18, 0.08, "Uptake flux", f"{float(top['uptake_flux']):.2f}", "#ffffff")
    draw_metric_pill(ax, 0.45, 0.46, 0.22, 0.08, "Crossfeeding support score", f"{float(top['crossfeeding_support_score']):.2f}", "#ffffff")
    draw_metric_pill(ax, 0.69, 0.46, 0.22, 0.08, "Biomass association score", f"{float(pathway_row['biomass_association_score']):.3f}", "#ffffff")

    ax.text(0.05, 0.405, "Interpretation", fontsize=12.5, weight="bold", color=TEXT, va="top")
    interp = (
        f"Adenosine in the shared extracellular pool ({top['shared_metabolite_id']}) is supplied mainly by "
        f"{species_display_name(top_producer['species_from'])}. Alistipes shahii takes it up through "
        f"{uptake['reaction_id']} and then routes adn[c] through active nucleotide interconversion reactions "
        f"including PUNP1, ADK1, NDPK3, and ADNK1. This highlighted route supports a biomass-relevant nucleotide "
        f"precursor pool rather than directly quantifying biomass contribution."
    )
    ax.text(0.05, 0.365, wrap(interp, 120), fontsize=9.1, color=TEXT, va="top", ha="left")

    note = (
        f"Shared-pool caveat: this is not a unique one-to-one exchange. "
        f"{data['other_producers']} other producer(s) and {data['other_consumers']} other consumer(s) also interact with adn[e] "
        f"in the saved provenance edges. Crossfeeding support score = min(producer flux, uptake flux)."
    )
    ax.text(0.05, 0.12, wrap(note, 135), fontsize=8.6, color=MUTED, va="top", ha="left")
    upstream_paths = [
        PROV_DIR.parent / "community_exchange_fluxes_by_diet.csv",
        PROV_DIR.parent / "community_species_connector_fluxes_by_diet.csv",
        PROV_DIR.parent / "community_reaction_fluxes_by_diet.csv",
        PROV_DIR.parent / "community_biomass_associated_pathways_by_diet.csv",
    ]
    ax.text(
        0.94,
        0.12,
        source_box_text(
            [config["clean_csv"], config["pathway_csv"], data["edges_csv"], data["steps_csv"]],
            upstream_paths=upstream_paths,
            width=54,
        ),
        fontsize=8.0,
        color=MUTED,
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#f7f4ec", "edgecolor": "#d8d2c3", "linewidth": 0.8},
    )


def main():
    for config in SPECIES_CONFIG:
        fig, ax = plt.subplots(figsize=(14, 9.1))
        fig.patch.set_facecolor(BG)
        draw_species_panel(ax, config)
        fig.tight_layout()
        fig.savefig(config["out_png"], dpi=240, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {config['out_png']}")

    fig, axes = plt.subplots(2, 1, figsize=(16, 17))
    fig.patch.set_facecolor(BG)
    for ax, config in zip(axes, SPECIES_CONFIG):
        draw_species_panel(ax, config)
    fig.tight_layout(h_pad=1.8)
    fig.savefig(COMBINED_OUT, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {COMBINED_OUT}")

    fig, ax = plt.subplots(figsize=(16, 9))
    fig.patch.set_facecolor(BG)
    draw_alistipes_mechanistic_chain(ax)
    fig.tight_layout()
    fig.savefig(ALISTIPES_CHAIN_OUT, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {ALISTIPES_CHAIN_OUT}")


if __name__ == "__main__":
    main()
