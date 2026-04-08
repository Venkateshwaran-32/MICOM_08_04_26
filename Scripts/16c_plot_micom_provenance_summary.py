from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_SUMMARY = (
    PROJECT_ROOT
    / "Results"
    / "micom"
    / "provenance"
    / "proper_age_bins"
    / "micom_provenance_case_overview.csv"
)
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "provenance"
FIG_DIR.mkdir(parents=True, exist_ok=True)

AGE_ORDER = ["61_70", "71_80", "81_plus"]
DIET_ORDER = ["high_fiber", "western"]
CASE_COLUMNS = [f"{age}|{diet}" for age in AGE_ORDER for diet in DIET_ORDER]
METRIC_ROWS = [
    ("top_pathway_score", "Top pathway score"),
    ("anchor_uptake_flux", "Anchor uptake by focal species"),
    ("crossfeeding_support_score", "Producer-supported uptake"),
]

BG = "#fcfbf7"
TEXT = "#173042"
SUBTLE = "#5c7280"
GRID = "#ddd6c8"
HF_COLOR = "#3e7d5e"
WESTERN_COLOR = "#b45f04"


def short_case_label(case: str) -> str:
    age, diet = case.split("|", 1)
    diet_short = "HF" if diet == "high_fiber" else "W"
    return f"{age}\n{diet_short}"


def display_species_name(species: str) -> str:
    base = str(species).replace("__AGORA1_03", "")
    parts = base.split("_")
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return base.replace("_", " ")


def compact_text(value: str, width: int = 18) -> str:
    text = "" if pd.isna(value) else str(value).replace("_", " ").strip()
    if not text:
        return "none"
    words = text.split()
    lines = []
    current = []
    length = 0
    for word in words:
        add_len = len(word) + (1 if current else 0)
        if current and length + add_len > width:
            lines.append(" ".join(current))
            current = [word]
            length = len(word)
        else:
            current.append(word)
            length += add_len
    if current:
        lines.append(" ".join(current))
    return "\n".join(lines)


def build_metric_heatmap(sub: pd.DataFrame) -> pd.DataFrame:
    work = sub.copy()
    work["case"] = work["age_group"].astype(str) + "|" + work["diet"].astype(str)
    rows = []
    for metric, label in METRIC_ROWS:
        row = work.set_index("case")[metric].reindex(CASE_COLUMNS).fillna(0.0)
        rows.append(row.rename(label))
    return pd.DataFrame(rows)


def build_case_table(sub: pd.DataFrame) -> pd.DataFrame:
    work = sub.copy()
    work["case"] = work["age_group"].astype(str) + "|" + work["diet"].astype(str)
    work = work.set_index("case").reindex(CASE_COLUMNS).reset_index()
    work["Case"] = work["case"].map(short_case_label)
    work["Anchor"] = work["anchor_metabolite_name"].fillna(work["anchor_metabolite_id"]).map(lambda v: compact_text(v, 14))
    work["Origin"] = work["likely_origin"].map(lambda v: compact_text(v, 18))
    work["Producer"] = work["top_producer_species"].map(lambda v: compact_text(v, 20))
    return work[["Case", "Anchor", "Origin", "Producer"]]


def summarize_species(sub: pd.DataFrame) -> str:
    top_row = sub.sort_values("top_pathway_score", ascending=False).iloc[0]
    nonzero_cf = int((sub["crossfeeding_support_score"] > 1e-9).sum())
    anchors = ", ".join(sorted(set(sub["anchor_metabolite_name"].fillna(sub["anchor_metabolite_id"]).astype(str))))
    return (
        "Summary\n"
        f"1. Top pathway stays {top_row['top_pathway']} across all cases.\n"
        f"2. Strongest pathway score: {top_row['age_group']} {top_row['diet']} ({top_row['top_pathway_score']:.2f}).\n"
        f"3. Anchor metabolites used here: {anchors}; producer-supported uptake is nonzero in {nonzero_cf}/6 cases."
    )


def draw_heatmap(ax, heat: pd.DataFrame):
    im = ax.imshow(heat.values, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(range(heat.shape[1]))
    ax.set_xticklabels([short_case_label(c) for c in heat.columns], fontsize=9)
    ax.set_yticks(range(heat.shape[0]))
    ax.set_yticklabels(heat.index, fontsize=10)
    ax.set_title("Quantitative provenance summary", loc="left", fontsize=12, color=TEXT, weight="bold", pad=10)
    ax.tick_params(length=0)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            val = heat.iloc[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=9, color="#c62828")
    for spine in ax.spines.values():
        spine.set_visible(False)
    return im


def draw_case_table(ax, case_table: pd.DataFrame):
    ax.axis("off")
    ax.set_title("Likely source of anchor metabolite by case", loc="left", fontsize=12, color=TEXT, weight="bold", pad=10)
    tbl = ax.table(
        cellText=case_table.values,
        colLabels=case_table.columns,
        cellLoc="left",
        colLoc="left",
        bbox=[0.0, 0.05, 1.0, 0.9],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#d8d2c3")
        cell.set_linewidth(0.6)
        if row == 0:
            cell.set_facecolor("#efe6d3")
            cell.set_text_props(weight="bold", color=TEXT)
        else:
            cell.set_facecolor("#fffdf8")
            if col == 0:
                cell.set_text_props(color=TEXT, weight="bold")
    tbl.scale(1.0, 1.35)


def plot_species(sub: pd.DataFrame, species: str):
    out_path = FIG_DIR / f"{species}__micom_provenance_summary.png"
    heat = build_metric_heatmap(sub)
    case_table = build_case_table(sub)
    summary_text = summarize_species(sub)

    fig = plt.figure(figsize=(14, 8.6), facecolor=BG)
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.45, 1.15],
        height_ratios=[1.0, 0.32],
        left=0.05,
        right=0.97,
        top=0.88,
        bottom=0.08,
        wspace=0.18,
        hspace=0.18,
    )

    ax_heat = fig.add_subplot(gs[0, 0])
    ax_table = fig.add_subplot(gs[0, 1])
    ax_footer = fig.add_subplot(gs[1, :])
    ax_footer.axis("off")

    im = draw_heatmap(ax_heat, heat)
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.035, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("Metric value", fontsize=9, color=TEXT)

    draw_case_table(ax_table, case_table)

    fig.suptitle(
        f"MICOM provenance dashboard: {display_species_name(species)}",
        fontsize=18,
        color=TEXT,
        weight="bold",
        y=0.965,
    )
    fig.text(
        0.05,
        0.925,
        "Rows summarize the focal species pathway score, anchor-metabolite uptake, and uptake plausibly supported by another species. "
        "Right table names the anchor metabolite and its likely source in each case.",
        fontsize=10,
        color=SUBTLE,
    )

    ax_footer.text(
        0.0,
        0.88,
        summary_text,
        fontsize=9,
        color=TEXT,
        va="top",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#fff8e8", "edgecolor": "#d8c9a8", "linewidth": 0.8},
    )
    ax_footer.text(
        0.995,
        0.88,
        "Source\nmicom_provenance_case_overview.csv\nDerived by 16b_map_micom_flux_provenance.py",
        fontsize=8.5,
        color=SUBTLE,
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#f7f4ec", "edgecolor": "#d8d2c3", "linewidth": 0.8},
    )
    ax_footer.text(
        0.0,
        0.18,
        "Producer-supported uptake is the portion of anchor-metabolite uptake that the script can plausibly match to secretion "
        "from another species. Cases remain discrete age-group and diet snapshots, not a continuous trajectory.",
        fontsize=8.5,
        color=SUBTLE,
        va="bottom",
    )

    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    if not IN_SUMMARY.exists():
        raise FileNotFoundError(f"Missing MICOM provenance summary input: {IN_SUMMARY}")

    df = pd.read_csv(IN_SUMMARY)
    required = {
        "age_group",
        "diet",
        "species",
        "top_pathway",
        "top_pathway_score",
        "anchor_metabolite_id",
        "anchor_metabolite_name",
        "anchor_uptake_flux",
        "crossfeeding_support_score",
        "likely_origin",
        "top_producer_species",
    }
    if not required.issubset(df.columns):
        raise ValueError(f"MICOM provenance summary table must contain columns: {sorted(required)}")

    df["age_group"] = pd.Categorical(df["age_group"], categories=AGE_ORDER, ordered=True)
    df["diet"] = pd.Categorical(df["diet"], categories=DIET_ORDER, ordered=True)
    df = df.sort_values(["species", "age_group", "diet"]).reset_index(drop=True)

    for species, sub in df.groupby("species", sort=False):
        plot_species(sub.copy(), species)


if __name__ == "__main__":
    main()
