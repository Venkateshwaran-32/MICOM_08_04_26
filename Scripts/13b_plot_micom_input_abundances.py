from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_CSV = PROJECT_ROOT / "Data" / "Processed" / "micom" / "agegroup_median_taxonomy_for_micom.csv"
FIG_DIR = PROJECT_ROOT / "Results" / "figures" / "micom" / "growth" / "proper_age_bins"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_FIG = FIG_DIR / "micom_input_abundances_by_agegroup_proper_age_bins.png"
MIN_ABUNDANCE = 1e-12
NOTE_BOX_X = 0.9825
NOTE_BOX_Y = 0.05


def pretty_age_group(value: str) -> str:
    mapping = {
        "21_40": "Ages 21-40",
        "41_60": "Ages 41-60",
        "61_70": "Ages 61-70",
        "71_80": "Ages 71-80",
        "81_plus": "Ages 81+",
    }
    return mapping.get(str(value), str(value))


def build_palette(labels: list[str]) -> dict[str, tuple]:
    cmap = plt.get_cmap("tab10")
    return {label: cmap(i % 10) for i, label in enumerate(sorted(labels))}


def main():
    if not IN_CSV.exists():
        raise FileNotFoundError(f"Missing MICOM abundance table: {IN_CSV}. Run Script 13 first.")

    df = pd.read_csv(IN_CSV)
    required = {"sample_id", "table2_taxon", "abundance"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input table missing columns: {sorted(required)}")

    df["sample_id"] = df["sample_id"].astype(str)
    df["table2_taxon"] = df["table2_taxon"].astype(str)
    df["abundance"] = df["abundance"].astype(float)

    species_order = sorted(df["table2_taxon"].drop_duplicates().tolist())
    palette = build_palette(species_order)
    age_groups = sorted(df["sample_id"].drop_duplicates().tolist())

    fig, axes = plt.subplots(len(age_groups), 1, figsize=(14, max(12, 2.4 * len(age_groups))), sharex=True)
    if len(age_groups) == 1:
        axes = [axes]

    for ax, age_group in zip(axes, age_groups):
        sub = df[df["sample_id"] == age_group].copy()
        if sub.empty:
            ax.axis("off")
            continue

        sub = sub.sort_values("abundance", ascending=True)
        colors = [palette[label] for label in sub["table2_taxon"]]
        bars = ax.barh(sub["table2_taxon"], sub["abundance"], color=colors)

        n_nonzero = int((sub["abundance"] > MIN_ABUNDANCE).sum())
        n_total = len(sub)
        ax.set_title(f"{pretty_age_group(age_group)} | {n_nonzero}/{n_total} nonzero abundance")
        ax.set_ylabel("")
        ax.grid(axis="x", alpha=0.25)

        for bar, value in zip(bars, sub["abundance"]):
            x = float(bar.get_width())
            y = bar.get_y() + bar.get_height() / 2.0
            label = f"{value:.5f}"
            if value > MIN_ABUNDANCE:
                ax.text(x + 0.01, y, label, va="center", ha="left", fontsize=8)
            else:
                ax.text(0.005, y, label, va="center", ha="left", fontsize=8, color="#666666")

        zero_species = sub.loc[sub["abundance"] <= MIN_ABUNDANCE, "table2_taxon"].tolist()
        if zero_species:
            note = "Zero abundance:\n" + "\n".join(f"{idx}. {name}" for idx, name in enumerate(zero_species, start=1))
        else:
            note = "Zero abundance:\nNone"
        ax.text(
            NOTE_BOX_X,
            NOTE_BOX_Y,
            note,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            multialignment="left",
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.88, "edgecolor": "#888888"},
        )

    axes[-1].set_xlabel("Normalized median abundance")
    axes[-1].set_xlim(0.0, max(0.82, float(df["abundance"].max()) * 1.18))

    handles = [
        plt.Line2D([0], [0], color=palette[label], linewidth=6, label=label)
        for label in sorted(palette)
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.95),
        title="Species",
        frameon=True,
        fontsize=8,
        title_fontsize=9,
        ncol=3,
        columnspacing=1.2,
        handlelength=1.8,
    )
    fig.suptitle("MICOM Input Abundances by Age Group", y=0.99, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.84])
    fig.savefig(OUT_FIG, dpi=220)
    plt.close(fig)
    print(f"Saved: {OUT_FIG}")


if __name__ == "__main__":
    main()
