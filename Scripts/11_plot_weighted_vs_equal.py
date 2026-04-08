from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FBA_DIR = PROJECT_ROOT / "Results" / "fba"
SCENARIO_DIR = FBA_DIR / "scenarios"
FIG_DIR = PROJECT_ROOT / "Results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

IN_COMPARE = SCENARIO_DIR / "sg90_median_vs_equal_exchange_flux_comparison_by_diet.csv"

OUT_TOPDELTA = FIG_DIR / "sg90_vs_equal_flux_delta_top25.png"
OUT_AGE_HEAT = FIG_DIR / "agegroup_vs_equal_flux_delta_heatmap.png"


def make_top_delta_plot(df: pd.DataFrame):
    import matplotlib.pyplot as plt

    sub = df[df["scenario"] == "sg90_median"].copy()
    if sub.empty:
        print("Skipping SG90 top-delta plot: no sg90_median rows.")
        return
    top = sub.sort_values("abs_delta_flux_vs_equal", ascending=False).head(25).copy()
    top["label"] = top["diet"] + " | " + top["exchange_id"]
    top = top.sort_values("delta_flux_vs_equal")

    fig, ax = plt.subplots(figsize=(11, 8))
    colors = ["#2E7D32" if v > 0 else "#C62828" for v in top["delta_flux_vs_equal"]]
    ax.barh(top["label"], top["delta_flux_vs_equal"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_title("Top Community Exchange Flux Shifts: SG90 Median vs Equal")
    ax.set_xlabel("delta flux (scenario - equal)")
    fig.tight_layout()
    fig.savefig(OUT_TOPDELTA, dpi=220)
    plt.close(fig)
    print(f"Saved: {OUT_TOPDELTA}")


def make_agegroup_heatmap(df: pd.DataFrame):
    import matplotlib.pyplot as plt
    import numpy as np

    if "scenario" not in df.columns:
        print("Skipping age-group heatmap: no scenario column in comparison table.")
        return
    sub = df[df["scenario"].astype(str).str.startswith("sg90_age_")].copy()
    if sub.empty:
        print("Skipping age-group heatmap: no sg90_age_* scenarios found.")
        return

    # Top exchanges by aggregate absolute delta to keep heatmap readable.
    top_ex = (
        sub.groupby("exchange_id", as_index=False)["abs_delta_flux_vs_equal"]
        .sum()
        .sort_values("abs_delta_flux_vs_equal", ascending=False)
        .head(30)["exchange_id"]
        .tolist()
    )
    heat = (
        sub[sub["exchange_id"].isin(top_ex)]
        .groupby(["scenario", "exchange_id"], as_index=False)["delta_flux_vs_equal"]
        .mean()
        .pivot(index="exchange_id", columns="scenario", values="delta_flux_vs_equal")
        .fillna(0.0)
    )
    if heat.empty:
        print("Skipping age-group heatmap: empty after filtering.")
        return

    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(heat.values, aspect="auto", cmap="coolwarm")
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=45, ha="right")
    ax.set_title("Age-group vs Equal: Exchange Flux Delta Heatmap")
    fig.colorbar(im, ax=ax, label="delta flux (scenario - equal)")
    fig.tight_layout()
    fig.savefig(OUT_AGE_HEAT, dpi=220)
    plt.close(fig)
    print(f"Saved: {OUT_AGE_HEAT}")


def main():
    if not IN_COMPARE.exists():
        raise FileNotFoundError(f"Missing input: {IN_COMPARE}. Run script 10 first.")

    df = pd.read_csv(IN_COMPARE)
    make_top_delta_plot(df)
    make_agegroup_heatmap(df)


if __name__ == "__main__":
    main()
