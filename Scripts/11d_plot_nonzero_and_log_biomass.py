from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCENARIO_DIR = PROJECT_ROOT / "Results" / "fba" / "scenarios"
FIG_DIR = PROJECT_ROOT / "Results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_COUNT = FIG_DIR / "allcohort_agebin_nonzero_species_count_by_diet.png"
OUT_LOG = FIG_DIR / "allcohort_agebin_species_biomass_log10_by_diet.png"
OUT_CSV = PROJECT_ROOT / "Results" / "fba" / "allcohort_agebin_nonzero_species_count_by_diet.csv"

AGE_ORDER = ["21_40", "41_60", "61_70", "71_80", "81_90"]
DIET_ORDER = ["western", "high_fiber"]
EPS = 1e-9


def short_name(full_name: str) -> str:
    base = str(full_name).replace("__", "_")
    parts = [p for p in base.split("_") if p]
    if not parts:
        return str(full_name)
    genus = parts[0]
    species = ""
    for token in parts[1:]:
        if token.isalpha() and token.lower() == token:
            species = token
            break
    return f"{genus} {species}".strip()


def short_abbrev(name: str) -> str:
    parts = str(name).split()
    if len(parts) >= 2:
        return f"{parts[0][0]}. {parts[1]}"
    return str(name)


def main():
    files = sorted(SCENARIO_DIR.glob("allcohort_agebin_*_species_biomass_by_diet.csv"))
    if not files:
        raise RuntimeError("No allcohort age-bin species files found. Run 09b first.")

    long_rows = []
    for fp in files:
        age = fp.name.replace("allcohort_agebin_", "").replace("_species_biomass_by_diet.csv", "")
        df = pd.read_csv(fp)
        df["age_group"] = age
        df["species_short"] = df["species"].apply(short_name)
        df["log10_biomass_flux_plus_eps"] = (df["biomass_flux"].astype(float) + EPS).map(lambda x: np.log10(x))
        df["nonzero"] = df["biomass_flux"].astype(float) > EPS
        long_rows.append(df)

    d = pd.concat(long_rows, ignore_index=True)
    d["age_group"] = pd.Categorical(d["age_group"], categories=AGE_ORDER, ordered=True)

    count_df = (
        d.groupby(["age_group", "diet"], as_index=False)["nonzero"]
        .sum()
        .rename(columns={"nonzero": "n_nonzero_species"})
        .sort_values(["age_group", "diet"])
    )
    count_df.to_csv(OUT_CSV, index=False)

    # Build species labels so each bar shows exactly which species are nonzero.
    nonzero_species = (
        d[d["nonzero"]]
        .groupby(["age_group", "diet"], as_index=False)["species_short"]
        .agg(lambda s: ", ".join(sorted(set(map(str, s)))))
    )
    nonzero_species["species_short"] = nonzero_species["species_short"].apply(
        lambda txt: ", ".join(short_abbrev(x.strip()) for x in str(txt).split(",") if x.strip())
    )
    species_label_map = {
        (str(r["age_group"]), str(r["diet"])): str(r["species_short"])
        for _, r in nonzero_species.iterrows()
    }

    # Figure 1: nonzero species count by age bin and diet.
    # Keep bar chart clean; show species labels in a separate table panel.
    piv_count = count_df.pivot(index="age_group", columns="diet", values="n_nonzero_species").reindex(AGE_ORDER).fillna(0.0)
    x = np.arange(len(AGE_ORDER))
    width = 0.35
    fig = plt.figure(figsize=(12, 7))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[3, 1])
    ax = fig.add_subplot(gs[0, 0])

    bars_w = ax.bar(x - width / 2, piv_count["western"].values, width=width, label="western", color="#C62828")
    bars_h = ax.bar(x + width / 2, piv_count["high_fiber"].values, width=width, label="high_fiber", color="#2E7D32")

    for bars, diet in [(bars_w, "western"), (bars_h, "high_fiber")]:
        for i, b in enumerate(bars):
            y = float(b.get_height())
            ax.text(
                b.get_x() + b.get_width() / 2,
                y + 0.05,
                f"{int(round(y))}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    ax.set_title("Species with Positive Growth by Age Bin and Diet")
    ax.set_xlabel("age bin")
    ax.set_ylabel("nonzero species count")
    ax.set_xticks(x)
    ax.set_xticklabels(AGE_ORDER)
    ax.legend(frameon=False)
    ax.set_ylim(0, max(1.0, float(piv_count.values.max()) + 1.0))

    ax_tbl = fig.add_subplot(gs[1, 0])
    ax_tbl.axis("off")
    table_rows = []
    for age in AGE_ORDER:
        table_rows.append(
            [
                age,
                species_label_map.get((age, "western"), "none"),
                species_label_map.get((age, "high_fiber"), "none"),
            ]
        )
    tbl = ax_tbl.table(
        cellText=table_rows,
        colLabels=["age bin", "western nonzero species", "high_fiber nonzero species"],
        cellLoc="left",
        loc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 1.2)

    fig.tight_layout()
    fig.savefig(OUT_COUNT, dpi=220)
    plt.close(fig)

    # Figure 2: log10 biomass flux heatmap-like matrix via pivot + imshow for readability.
    for diet in DIET_ORDER:
        sub = d[d["diet"] == diet].copy()
        if sub.empty:
            continue
        heat = (
            sub.pivot_table(index="species_short", columns="age_group", values="log10_biomass_flux_plus_eps", aggfunc="mean")
            .reindex(columns=AGE_ORDER)
            .fillna(np.log10(EPS))
        )
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(heat.values, aspect="auto", cmap="viridis")
        ax.set_yticks(range(len(heat.index)))
        ax.set_yticklabels(heat.index)
        ax.set_xticks(range(len(heat.columns)))
        ax.set_xticklabels(heat.columns)
        ax.set_title(f"Species Biomass Flux (log10) by Age Bin | {diet}")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("log10(biomass_flux + 1e-9)")
        fig.tight_layout()
        out = OUT_LOG.with_name(f"allcohort_agebin_species_biomass_log10_{diet}.png")
        fig.savefig(out, dpi=220)
        plt.close(fig)

    print(f"Saved: {OUT_CSV}")
    print(f"Saved: {OUT_COUNT}")
    print(f"Saved: {OUT_LOG.with_name('allcohort_agebin_species_biomass_log10_western.png')}")
    print(f"Saved: {OUT_LOG.with_name('allcohort_agebin_species_biomass_log10_high_fiber.png')}")


if __name__ == "__main__":
    main()
