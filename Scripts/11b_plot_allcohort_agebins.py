from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_COMPARE = PROJECT_ROOT / "Results" / "fba" / "allcohort_agebin_vs_equal_exchange_comparisons.csv"
FIG_DIR = PROJECT_ROOT / "Results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_TOPDELTA = FIG_DIR / "allcohort_agebin_vs_equal_flux_delta_top30.png"
OUT_HEAT = FIG_DIR / "allcohort_agebin_vs_equal_flux_delta_heatmap_top40.png"


def make_top_delta_plot(df: pd.DataFrame):
    top = df.sort_values("abs_delta_flux_vs_equal", ascending=False).head(30).copy()
    if top.empty:
        print("Skipping top-delta plot: no rows.")
        return

    top["label"] = top["scenario"].astype(str) + " | " + top["diet"].astype(str) + " | " + top["exchange_id"].astype(str)
    top = top.sort_values("delta_flux_vs_equal")

    fig, ax = plt.subplots(figsize=(12, 9))
    colors = ["#2E7D32" if v > 0 else "#C62828" for v in top["delta_flux_vs_equal"]]
    ax.barh(top["label"], top["delta_flux_vs_equal"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_title("Top Exchange Flux Shifts: All-Cohort Age Bins vs Equal")
    ax.set_xlabel("delta flux (age-bin scenario - allcohort equal)")
    fig.tight_layout()
    fig.savefig(OUT_TOPDELTA, dpi=220)
    plt.close(fig)
    print(f"Saved: {OUT_TOPDELTA}")


def make_heatmap(df: pd.DataFrame):
    if df.empty:
        print("Skipping heatmap: no rows.")
        return

    top_ex = (
        df.groupby("exchange_id", as_index=False)["abs_delta_flux_vs_equal"]
        .sum()
        .sort_values("abs_delta_flux_vs_equal", ascending=False)
        .head(40)["exchange_id"]
        .tolist()
    )

    sub = df[df["exchange_id"].isin(top_ex)].copy()
    sub["scenario_diet"] = sub["scenario"].astype(str) + "|" + sub["diet"].astype(str)
    heat = (
        sub.groupby(["exchange_id", "scenario_diet"], as_index=False)["delta_flux_vs_equal"]
        .mean()
        .pivot(index="exchange_id", columns="scenario_diet", values="delta_flux_vs_equal")
        .fillna(0.0)
    )

    if heat.empty:
        print("Skipping heatmap: empty after filtering.")
        return

    fig, ax = plt.subplots(figsize=(13, 11))
    im = ax.imshow(heat.values, aspect="auto", cmap="coolwarm")
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=45, ha="right")
    ax.set_title("All-Cohort Age-Bin vs Equal: Exchange Flux Delta Heatmap")
    fig.colorbar(im, ax=ax, label="delta flux")
    fig.tight_layout()
    fig.savefig(OUT_HEAT, dpi=220)
    plt.close(fig)
    print(f"Saved: {OUT_HEAT}")


def main():
    if not IN_COMPARE.exists():
        raise FileNotFoundError(f"Missing input: {IN_COMPARE}. Run script 10b first.")

    df = pd.read_csv(IN_COMPARE)
    need = {"scenario", "diet", "exchange_id", "delta_flux_vs_equal", "abs_delta_flux_vs_equal"}
    if not need.issubset(df.columns):
        raise ValueError(f"Input comparison table missing required columns: {sorted(need)}")

    make_top_delta_plot(df)
    make_heatmap(df)


if __name__ == "__main__":
    main()
