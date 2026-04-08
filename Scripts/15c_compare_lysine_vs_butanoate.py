from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LYS_BIO_SUMMARY = PROJECT_ROOT / "Results" / "micom" / "lysine_butyrate" / "lysine_biosynthesis_flux_summary.csv"
LYS_BUT_SUMMARY = PROJECT_ROOT / "Results" / "micom" / "lysine_butyrate" / "lysine_to_butyrate_flux_summary.csv"
STEP_SUMMARY = PROJECT_ROOT / "Results" / "micom" / "lysine_butyrate" / "lysine_to_butyrate_reaction_step_summary.csv"

OUT_DIR = PROJECT_ROOT / "Results" / "micom" / "lysine_butyrate"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "lysine_butyrate"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_COMPARE = OUT_DIR / "lysine_vs_butanoate_comparison_by_agegroup_diet.csv"
OUT_HEAT = FIG_DIR / "lysine_vs_butanoate_heatmap_by_agegroup_diet.png"
OUT_STEP_HEAT = FIG_DIR / "butanoate_downstream_steps_heatmap_by_agegroup_diet.png"

AGE_ORDER = ["21_40", "41_60", "61_70", "71_80", "81_plus"]
DIET_ORDER = ["western", "high_fiber"]
BTCOADH_LABEL = "Butanoate metabolism (BTCOADH: crotonoyl-CoA to butyryl-CoA)"
BTCOAACCOAT_LABEL = "Butanoate metabolism (BTCOAACCOAT: butyryl-CoA:acetate CoA-transferase)"
COMPARISON_SUMMARY = (
    "Summary\n"
    "1. Lysine metabolism stays high across all age groups.\n"
    "2. Older bins, especially 71_80, retain strong lysine flux.\n"
    "3. Butanoate signal appears mainly in older groups, making\n"
    "   71_80 the clearest lysine-to-butanoate contrast."
)
COMPARISON_SOURCES = (
    "Sources\n"
    "lysine_biosynthesis_flux_summary.csv\n"
    "lysine_to_butyrate_flux_summary.csv\n"
    "lysine_to_butyrate_reaction_step_summary.csv"
)
STEP_SUMMARY_TEXT = (
    "Summary\n"
    "1. Only downstream butanoate steps are nonzero here.\n"
    "2. BTCOADH and BTCOAACCOAT carry the full visible signal.\n"
    "3. The nonzero step signal is concentrated in older age bins."
)
STEP_SOURCES = (
    "Sources\n"
    "lysine_to_butyrate_reaction_step_summary.csv"
)


def ordered_columns():
    return [f"{age}|{diet}" for age in AGE_ORDER for diet in DIET_ORDER]


def build_comparison_table(lys_bio_df: pd.DataFrame, lys_but_df: pd.DataFrame, step_df: pd.DataFrame) -> pd.DataFrame:
    lys = lys_bio_df[lys_bio_df["pathway"] == "Lysine metabolism"].copy()
    but = lys_but_df[lys_but_df["pathway"] == "Butanoate metabolism"].copy()

    merged = (
        lys[["age_group", "diet", "sum_abs_flux"]]
        .rename(columns={"sum_abs_flux": "lysine_metabolism_flux"})
        .merge(
            but[["age_group", "diet", "sum_abs_flux"]].rename(columns={"sum_abs_flux": "butanoate_metabolism_flux"}),
            on=["age_group", "diet"],
            how="outer",
        )
        .fillna(0.0)
    )

    btcoadh = (
        step_df[step_df["reaction_step"] == BTCOADH_LABEL][["age_group", "diet", "sum_abs_flux"]]
        .rename(columns={"sum_abs_flux": "btcoadh_flux"})
    )
    btcoaaccoat = (
        step_df[step_df["reaction_step"] == BTCOAACCOAT_LABEL][["age_group", "diet", "sum_abs_flux"]]
        .rename(columns={"sum_abs_flux": "btcoaaccoat_flux"})
    )

    merged = (
        merged.merge(btcoadh, on=["age_group", "diet"], how="left")
        .merge(btcoaaccoat, on=["age_group", "diet"], how="left")
        .fillna(0.0)
    )
    merged["delta_lysine_minus_butanoate"] = merged["lysine_metabolism_flux"] - merged["butanoate_metabolism_flux"]
    merged["downstream_butanoate_step_flux"] = merged["btcoadh_flux"] + merged["btcoaaccoat_flux"]
    merged["butanoate_nonzero_flag"] = merged["butanoate_metabolism_flux"] > 1e-9

    merged["age_group"] = pd.Categorical(merged["age_group"], categories=AGE_ORDER, ordered=True)
    merged["diet"] = pd.Categorical(merged["diet"], categories=DIET_ORDER, ordered=True)
    merged = merged.sort_values(["age_group", "diet"]).reset_index(drop=True)
    return merged


def draw_heatmap(ax, data: pd.DataFrame, title: str, cmap: str, fmt: str = ".1f"):
    im = ax.imshow(data.values, aspect="auto", cmap=cmap)
    ax.set_title(title, fontsize=12, loc="left", weight="bold", color="#173042", pad=10)
    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels(data.columns, rotation=40, ha="right", fontsize=9)
    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels(data.index, fontsize=10)
    ax.tick_params(axis="both", length=0)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data.iloc[i, j]
            text = format(val, fmt)
            ax.text(j, i, text, ha="center", va="center", fontsize=8, color="#c62828")

    for spine in ax.spines.values():
        spine.set_visible(False)
    return im


def plot_comparison_heatmap(compare_df: pd.DataFrame, out_path: Path):
    cols = ordered_columns()
    work = compare_df.copy()
    work["age_diet"] = work["age_group"].astype(str) + "|" + work["diet"].astype(str)

    pivot_rows = []
    for metric_name, label in [
        ("lysine_metabolism_flux", "Lysine metabolism"),
        ("butanoate_metabolism_flux", "Butanoate metabolism"),
        ("delta_lysine_minus_butanoate", "Lysine - Butanoate delta"),
    ]:
        row = (
            work.set_index("age_diet")[metric_name]
            .reindex(cols)
            .fillna(0.0)
        )
        pivot_rows.append(row.rename(label))

    heat = pd.DataFrame(pivot_rows)

    fig, ax = plt.subplots(figsize=(14, 4.8))
    fig.patch.set_facecolor("#fcfbf7")
    im = draw_heatmap(
        ax,
        heat,
        "Lysine vs Butanoate Metabolism by Age Group and Diet",
        cmap="YlGnBu",
        fmt=".1f",
    )
    fig.text(0.07, 0.94, "Pathway-level sum abs flux; columns are age_group|diet", fontsize=10, color="#5c7280")
    fig.text(
        0.07,
        0.09,
        COMPARISON_SUMMARY,
        fontsize=9,
        color="#173042",
        va="top",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#fff8e8", "edgecolor": "#d8c9a8", "linewidth": 0.8},
    )
    fig.text(
        0.93,
        0.09,
        COMPARISON_SOURCES,
        fontsize=8,
        color="#5c7280",
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#f7f4ec", "edgecolor": "#d8d2c3", "linewidth": 0.8},
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    plt.tight_layout(rect=[0.03, 0.16, 0.98, 0.92])
    plt.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_step_heatmap(compare_df: pd.DataFrame, out_path: Path):
    cols = ordered_columns()
    work = compare_df.copy()
    work["age_diet"] = work["age_group"].astype(str) + "|" + work["diet"].astype(str)

    heat = pd.DataFrame(
        [
            work.set_index("age_diet")["btcoadh_flux"].reindex(cols).fillna(0.0).rename("BTCOADH"),
            work.set_index("age_diet")["btcoaaccoat_flux"].reindex(cols).fillna(0.0).rename("BTCOAACCOAT"),
        ]
    )

    fig, ax = plt.subplots(figsize=(14, 3.8))
    fig.patch.set_facecolor("#fcfbf7")
    im = draw_heatmap(
        ax,
        heat,
        "Downstream Butanoate Steps by Age Group and Diet",
        cmap="Oranges",
        fmt=".1f",
    )
    fig.text(0.07, 0.92, "Only nonzero downstream butanoate steps from the lysine-to-butyrate candidate branch", fontsize=10, color="#5c7280")
    fig.text(
        0.07,
        0.09,
        STEP_SUMMARY_TEXT,
        fontsize=9,
        color="#173042",
        va="top",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#fff8e8", "edgecolor": "#d8c9a8", "linewidth": 0.8},
    )
    fig.text(
        0.93,
        0.09,
        STEP_SOURCES,
        fontsize=8,
        color="#5c7280",
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#f7f4ec", "edgecolor": "#d8d2c3", "linewidth": 0.8},
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    plt.tight_layout(rect=[0.03, 0.17, 0.98, 0.89])
    plt.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main():
    if not LYS_BUT_SUMMARY.exists():
        raise FileNotFoundError(f"Missing input: {LYS_BUT_SUMMARY}")
    if not LYS_BIO_SUMMARY.exists():
        raise FileNotFoundError(f"Missing input: {LYS_BIO_SUMMARY}")
    if not STEP_SUMMARY.exists():
        raise FileNotFoundError(f"Missing input: {STEP_SUMMARY}")

    lys_bio_df = pd.read_csv(LYS_BIO_SUMMARY)
    lys_but_df = pd.read_csv(LYS_BUT_SUMMARY)
    step_df = pd.read_csv(STEP_SUMMARY)

    compare_df = build_comparison_table(lys_bio_df, lys_but_df, step_df)
    compare_df.to_csv(OUT_COMPARE, index=False)

    plot_comparison_heatmap(compare_df, OUT_HEAT)
    plot_step_heatmap(compare_df, OUT_STEP_HEAT)

    print(f"Saved comparison table: {OUT_COMPARE}")
    print(f"Saved heatmap: {OUT_HEAT}")
    print(f"Saved step heatmap: {OUT_STEP_HEAT}")


if __name__ == "__main__":
    main()
