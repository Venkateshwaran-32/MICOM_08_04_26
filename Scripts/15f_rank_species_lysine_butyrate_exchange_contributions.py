from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PATHWAY_DIR = PROJECT_ROOT / "Results" / "micom" / "pathway_flux" / "proper_age_bins"
GROWTH_DIR = PROJECT_ROOT / "Results" / "micom" / "growth" / "proper_age_bins"
PROC_MICOM = PROJECT_ROOT / "Data" / "Processed" / "micom"
OUT_DIR = PROJECT_ROOT / "Results" / "micom" / "lysine_butyrate"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "lysine_butyrate"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

IN_FLUX = PATHWAY_DIR / "reaction_fluxes_long_by_agegroup_diet.csv"
IN_GROWTH = GROWTH_DIR / "organism_growth_rates_by_agegroup_diet.csv"
IN_ABUND = PROC_MICOM / "agegroup_median_taxonomy_for_micom.csv"

OUT_CANONICAL = OUT_DIR / "species_lysine_butyrate_exchange_contributions_by_agegroup_diet.csv"
OUT_GROWERS = OUT_DIR / "species_lysine_butyrate_exchange_contributions_ranked_growers_only.csv"
OUT_FIG = FIG_DIR / "species_lysine_butyrate_exchange_contributions_by_agegroup_diet.png"

AGE_ORDER = ["21_40", "41_60", "61_70", "71_80", "81_plus"]
DIET_ORDER = ["high_fiber", "western"]
MIN_GROWTH = 1e-6

LYSINE_EXCHANGE_RXN = "EX_lys_L(e)"
BUTYRATE_EXCHANGE_RXN = "EX_but(e)"

BG = "#fcfbf7"
TEXT = "#173042"
GRID = "#ddd6c8"
LYSINE_COLOR = "#1f5a89"
BUTYRATE_COLOR = "#c76b29"
SUMMARY_TEXT = (
    "Summary\n"
    "1. Positive values export to the shared environment; negative values indicate uptake.\n"
    "2. Older bins show the clearest butyrate exchange signal.\n"
    "3. Panels show only species with positive MICOM growth at TF = 0.5."
)
SOURCE_TEXT = (
    "Sources\n"
    "reaction_fluxes_long_\n"
    "by_agegroup_diet.csv\n"
    "organism_growth_rates_\n"
    "by_agegroup_diet.csv\n"
    "species_lysine_butyrate_exchange_\n"
    "contributions_by_agegroup_diet.csv"
)


def pretty_age_group(value: str) -> str:
    mapping = {
        "21_40": "Ages 21-40",
        "41_60": "Ages 41-60",
        "61_70": "Ages 61-70",
        "71_80": "Ages 71-80",
        "81_plus": "Ages 81+",
    }
    return mapping.get(str(value), str(value))


def species_display_name(row: pd.Series) -> str:
    if pd.notna(row.get("table2_taxon")) and str(row["table2_taxon"]).strip():
        return str(row["table2_taxon"])
    if pd.notna(row.get("profile_taxon")) and str(row["profile_taxon"]).strip():
        return str(row["profile_taxon"]).replace("s__", "").replace("_", " ")
    return str(row["id"])


def load_abundance() -> pd.DataFrame:
    abund = pd.read_csv(IN_ABUND).copy()
    required = {"sample_id", "id", "abundance"}
    if not required.issubset(abund.columns):
        raise ValueError(f"MICOM abundance file must contain columns: {sorted(required)}")

    abund = abund.rename(columns={"sample_id": "age_group"})
    keep = ["age_group", "id", "abundance", "table2_taxon", "profile_taxon", "model_species_id"]
    keep = [c for c in keep if c in abund.columns]
    abund = abund[keep].drop_duplicates(["age_group", "id"]).copy()
    abund["species_label"] = abund.apply(species_display_name, axis=1)
    return abund


def load_growth() -> pd.DataFrame:
    growth = pd.read_csv(IN_GROWTH).copy()
    required = {"age_group", "diet", "id", "growth_rate"}
    if not required.issubset(growth.columns):
        raise ValueError(f"MICOM growth file must contain columns: {sorted(required)}")
    growth["grows_in_diet"] = growth["growth_rate"] > MIN_GROWTH
    return growth


def build_exchange_summary() -> pd.DataFrame:
    flux = pd.read_csv(IN_FLUX).copy()
    required = {"id", "age_group", "diet", "reaction_id", "flux"}
    if not required.issubset(flux.columns):
        raise ValueError(f"MICOM flux file must contain columns: {sorted(required)}")

    flux["flux"] = pd.to_numeric(flux["flux"], errors="coerce").fillna(0.0)

    ex = (
        flux[flux["reaction_id"].isin([LYSINE_EXCHANGE_RXN, BUTYRATE_EXCHANGE_RXN])]
        .groupby(["age_group", "diet", "id", "reaction_id"], as_index=False)
        .agg(net_exchange_flux=("flux", "sum"))
    )
    pivot = (
        ex.pivot(index=["age_group", "diet", "id"], columns="reaction_id", values="net_exchange_flux")
        .reset_index()
        .fillna(0.0)
    )
    pivot.columns.name = None
    if LYSINE_EXCHANGE_RXN not in pivot.columns:
        pivot[LYSINE_EXCHANGE_RXN] = 0.0
    if BUTYRATE_EXCHANGE_RXN not in pivot.columns:
        pivot[BUTYRATE_EXCHANGE_RXN] = 0.0

    pivot = pivot.rename(
        columns={
            LYSINE_EXCHANGE_RXN: "net_lysine_exchange_flux_to_environment",
            BUTYRATE_EXCHANGE_RXN: "net_butyrate_exchange_flux_to_environment",
        }
    )
    pivot["lysine_export_flux_to_environment"] = pivot["net_lysine_exchange_flux_to_environment"].clip(lower=0.0)
    pivot["lysine_uptake_flux_from_environment"] = (-pivot["net_lysine_exchange_flux_to_environment"]).clip(lower=0.0)
    pivot["butyrate_export_flux_to_environment"] = pivot["net_butyrate_exchange_flux_to_environment"].clip(lower=0.0)
    pivot["butyrate_uptake_flux_from_environment"] = (-pivot["net_butyrate_exchange_flux_to_environment"]).clip(lower=0.0)
    pivot["exchange_salience"] = pivot[
        ["net_lysine_exchange_flux_to_environment", "net_butyrate_exchange_flux_to_environment"]
    ].abs().max(axis=1)
    return pivot


def build_canonical_table(abund: pd.DataFrame, growth: pd.DataFrame, exchange: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for age_group in sorted(abund["age_group"].astype(str).unique(), key=lambda x: AGE_ORDER.index(x)):
        age_abund = abund[abund["age_group"].astype(str) == age_group].copy()
        for diet in DIET_ORDER:
            block = age_abund.copy()
            block["diet"] = diet
            rows.append(block)
    base = pd.concat(rows, ignore_index=True)

    out = (
        base.merge(growth, on=["age_group", "diet", "id"], how="left", suffixes=("", "_growth"))
        .merge(exchange, on=["age_group", "diet", "id"], how="left")
    )

    out["growth_rate"] = out["growth_rate"].fillna(0.0)
    out["grows_in_diet"] = out["grows_in_diet"].where(out["grows_in_diet"].notna(), False).astype(bool)
    out["shown_in_growth_figure"] = out["growth_rate"] > MIN_GROWTH

    numeric_fill = [
        "net_lysine_exchange_flux_to_environment",
        "net_butyrate_exchange_flux_to_environment",
        "lysine_export_flux_to_environment",
        "lysine_uptake_flux_from_environment",
        "butyrate_export_flux_to_environment",
        "butyrate_uptake_flux_from_environment",
        "exchange_salience",
    ]
    for col in numeric_fill:
        out[col] = out[col].fillna(0.0)

    out["rank_exchange_salience_within_panel"] = (
        out.groupby(["age_group", "diet"])["exchange_salience"]
        .rank(method="dense", ascending=False)
        .astype(int)
    )

    out["age_group"] = pd.Categorical(out["age_group"], categories=AGE_ORDER, ordered=True)
    out["diet"] = pd.Categorical(out["diet"], categories=DIET_ORDER, ordered=True)
    out = out.sort_values(
        ["age_group", "diet", "shown_in_growth_figure", "rank_exchange_salience_within_panel", "species_label"],
        ascending=[True, True, False, True, True],
    ).reset_index(drop=True)
    return out


def build_growers_only_table(canonical: pd.DataFrame) -> pd.DataFrame:
    growers = canonical[canonical["shown_in_growth_figure"]].copy()
    growers = growers.sort_values(
        ["age_group", "diet", "rank_exchange_salience_within_panel", "species_label"]
    ).reset_index(drop=True)
    return growers


def format_value(value: float) -> str:
    value = abs(float(value))
    if value >= 1000:
        return f"{value / 1000:.1f}k"
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}"
    if value >= 1:
        return f"{value:.2f}"
    if value > 0:
        return f"{value:.3f}"
    return "0"


def add_value_labels(ax, bars, values):
    labels = []
    for val in values:
        if abs(val) <= 1e-12:
            labels.append("")
        else:
            labels.append(format_value(val))
    ax.bar_label(bars, labels=labels, padding=2, fontsize=7, color=TEXT)


def plot_ranked_panels(growers: pd.DataFrame, canonical: pd.DataFrame):
    fig, axes = plt.subplots(len(AGE_ORDER), len(DIET_ORDER), figsize=(16, 2.9 * len(AGE_ORDER)))
    fig.patch.set_facecolor(BG)

    for i, age_group in enumerate(AGE_ORDER):
        for j, diet in enumerate(DIET_ORDER):
            ax = axes[i, j]
            ax.set_facecolor(BG)
            sub = growers[(growers["age_group"].astype(str) == age_group) & (growers["diet"].astype(str) == diet)].copy()

            if sub.empty:
                ax.axis("off")
                continue

            sub = sub.sort_values(["exchange_salience", "species_label"], ascending=[True, True])
            y = np.arange(len(sub))
            h = 0.36

            lys_bars = ax.barh(
                y + h / 2,
                sub["net_lysine_exchange_flux_to_environment"],
                height=h,
                color=LYSINE_COLOR,
                label="Lysine",
            )
            but_bars = ax.barh(
                y - h / 2,
                sub["net_butyrate_exchange_flux_to_environment"],
                height=h,
                color=BUTYRATE_COLOR,
                label="Butyrate",
            )

            ax.set_yticks(y)
            ax.set_yticklabels(sub["species_label"], fontsize=8, color=TEXT)
            ax.set_title(f"{pretty_age_group(age_group)} | {diet} | TF = 0.5", fontsize=11, color=TEXT)
            ax.set_xlabel("Net exchange with shared environment", fontsize=9, color=TEXT)
            ax.tick_params(axis="x", labelsize=8, colors=TEXT)
            ax.grid(axis="x", color=GRID, linewidth=0.8, alpha=0.7)
            ax.set_axisbelow(True)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_color("#bfb6a5")
            ax.spines["bottom"].set_color("#bfb6a5")
            ax.axvline(0.0, color="#8a8172", linewidth=1.0, alpha=0.9)

            panel_max = max(
                float(sub["net_lysine_exchange_flux_to_environment"].abs().max()),
                float(sub["net_butyrate_exchange_flux_to_environment"].abs().max()),
                1.0,
            )
            ax.set_xlim(-panel_max * 1.18, panel_max * 1.18)

            add_value_labels(ax, lys_bars, sub["net_lysine_exchange_flux_to_environment"].tolist())
            add_value_labels(ax, but_bars, sub["net_butyrate_exchange_flux_to_environment"].tolist())

            nongrowing_active = canonical[
                (canonical["age_group"].astype(str) == age_group)
                & (canonical["diet"].astype(str) == diet)
                & (~canonical["shown_in_growth_figure"])
                & (canonical["exchange_salience"] > 0)
            ].copy()
            if not nongrowing_active.empty:
                note_lines = ["Non-growing with exchange:"]
                for row in nongrowing_active.sort_values("exchange_salience", ascending=False).head(4).itertuples(index=False):
                    note_lines.append(
                        f"{row.species_label}: Lys {format_value(row.net_lysine_exchange_flux_to_environment)}, "
                        f"But {format_value(row.net_butyrate_exchange_flux_to_environment)}"
                    )
                ax.text(
                    0.985,
                    0.03,
                    "\n".join(note_lines),
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=7,
                    bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.88, "edgecolor": "#888888"},
                )

    handles = [
        plt.Line2D([0], [0], color=LYSINE_COLOR, lw=8),
        plt.Line2D([0], [0], color=BUTYRATE_COLOR, lw=8),
    ]
    labels = ["Net lysine exchange", "Net butyrate exchange"]
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.995))
    fig.suptitle(
        "Species-Level MICOM Net Lysine And Butyrate Exchange With The Shared Environment",
        fontsize=15,
        color=TEXT,
        weight="bold",
        y=1.02,
    )
    fig.text(
        0.5,
        0.072,
        "Positive = export to shared environment. Negative = uptake from shared environment. "
        "Panels show only species with positive MICOM growth in the corresponding TF = 0.5 solution.",
        fontsize=9,
        color="#5c7280",
        ha="center",
    )
    fig.text(
        0.02,
        0.01,
        SUMMARY_TEXT,
        fontsize=8.5,
        color=TEXT,
        va="bottom",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#fff8e8", "edgecolor": "#d8c9a8", "linewidth": 0.8},
    )
    fig.text(
        0.98,
        0.01,
        SOURCE_TEXT,
        fontsize=8,
        color="#5c7280",
        ha="right",
        va="bottom",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#f7f4ec", "edgecolor": "#d8d2c3", "linewidth": 0.8},
    )
    plt.tight_layout(rect=[0, 0.14, 1, 0.97])
    plt.savefig(OUT_FIG, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main():
    for path in [IN_FLUX, IN_GROWTH, IN_ABUND]:
        if not path.exists():
            raise FileNotFoundError(f"Missing required input: {path}")

    abund = load_abundance()
    growth = load_growth()
    exchange = build_exchange_summary()

    canonical = build_canonical_table(abund, growth, exchange)
    growers = build_growers_only_table(canonical)

    canonical.to_csv(OUT_CANONICAL, index=False)
    growers.to_csv(OUT_GROWERS, index=False)
    plot_ranked_panels(growers, canonical)

    print(f"Saved: {OUT_CANONICAL}")
    print(f"Saved: {OUT_GROWERS}")
    print(f"Saved: {OUT_FIG}")


if __name__ == "__main__":
    main()
