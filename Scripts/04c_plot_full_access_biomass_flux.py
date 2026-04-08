import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors


# ---------- PATHS ----------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IN_CSV = PROJECT_ROOT / "Results" / "fba" / "community_species_biomass_flux_full_access.csv"
OUT_DIR = PROJECT_ROOT / "Results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FLUX_EPS = 1e-9


def prettify_species_name(text: str) -> str:
    base = str(text).replace("__", "_").replace("_AGORA1_03", "")
    parts = [p for p in base.split("_") if p]
    if len(parts) >= 2:
        genus = parts[0]
        species = parts[1]
        if genus == "Parabacteroides":
            return f"Parabact\nerium\n{species}"
        return f"{genus}\n{species}"
    return base.replace("_", " ")


def get_input_csv() -> Path:
    raw = os.environ.get("FULL_ACCESS_BIOMASS_CSV")
    return Path(raw) if raw else DEFAULT_IN_CSV


def output_path_for_input(input_csv: Path) -> Path:
    stem = input_csv.stem.replace("community_species_biomass_flux_", "")
    return OUT_DIR / f"community_species_biomass_flux_{stem}_top_growers_by_diet.png"


def plot_title(df: pd.DataFrame) -> str:
    if "uptake_bound_per_exchange" in df.columns:
        vals = df["uptake_bound_per_exchange"].dropna().unique()
        if len(vals) == 1:
            val = float(vals[0])
            label = str(int(val)) if val.is_integer() else str(val)
            return f"Top Growers Under Full-Access Community FBA (both diets; uptake bound = {label})"
    return "Top Growers Under Full-Access Community FBA"


def light_palette(n: int) -> list[str]:
    base = list(plt.get_cmap("Pastel1").colors) + list(plt.get_cmap("Pastel2").colors) + list(plt.get_cmap("Set3").colors)
    colors = []
    for i in range(n):
        rgb = base[i % len(base)]
        colors.append(mcolors.to_hex(rgb))
    return colors


def identical_diet_profiles(df: pd.DataFrame, flux_eps: float = FLUX_EPS) -> bool:
    diets = list(dict.fromkeys(df["diet"].astype(str)))
    if len(diets) != 2:
        return False
    left = (
        df[df["diet"] == diets[0]][["species", "biomass_flux"]]
        .copy()
        .sort_values("species")
        .reset_index(drop=True)
    )
    right = (
        df[df["diet"] == diets[1]][["species", "biomass_flux"]]
        .copy()
        .sort_values("species")
        .reset_index(drop=True)
    )
    if list(left["species"]) != list(right["species"]):
        return False
    return ((left["biomass_flux"] - right["biomass_flux"]).abs() <= flux_eps).all()


def main():
    in_csv = get_input_csv()
    out_png = output_path_for_input(in_csv)

    df = pd.read_csv(in_csv)
    if not {"diet", "species", "biomass_flux"}.issubset(df.columns):
        raise ValueError("Input CSV must contain diet, species, and biomass_flux columns.")

    df = df.dropna(subset=["biomass_flux"]).copy()
    if df.empty:
        raise ValueError("No biomass flux values found in the input CSV.")

    diets = list(dict.fromkeys(df["diet"].astype(str)))
    combine_diets = identical_diet_profiles(df)
    panel_diets = [diets[0]] if combine_diets else diets
    fig, axes = plt.subplots(1, len(panel_diets), figsize=(6 * len(panel_diets), 5), sharey=True)
    if len(panel_diets) == 1:
        axes = [axes]

    for ax, diet_name in zip(axes, panel_diets):
        plot_df = (
            df[df["diet"] == diet_name]
            .copy()
            .sort_values("biomass_flux", ascending=False)
        )
        plot_df = plot_df[plot_df["biomass_flux"] > FLUX_EPS].copy()

        if plot_df.empty:
            ax.set_title(f"{diet_name}\n(no positive growers)")
            ax.axis("off")
            continue

        plot_df["species_label"] = plot_df["species"].map(prettify_species_name)
        bars = ax.bar(
            plot_df["species_label"],
            plot_df["biomass_flux"],
            color=light_palette(len(plot_df)),
            edgecolor="#6b6b6b",
            linewidth=0.8,
        )
        ax.set_title("western_full_access and high_fiber_full_access" if combine_diets else diet_name)
        ax.set_xlabel("Species")
        ax.tick_params(axis="x", rotation=0)
        total_flux = float(df[df["diet"] == diet_name]["biomass_flux"].fillna(0.0).sum())
        ax.text(
            0.98,
            0.98,
            f"Total biomass flux: {total_flux:.2f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            color="black",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "#8a8a8a", "boxstyle": "round,pad=0.25"},
        )
        for bar, value in zip(bars, plot_df["biomass_flux"]):
            y = bar.get_height() * 0.55
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                y,
                f"{value:.2f}",
                ha="center",
                va="center",
                color="black",
                fontsize=8,
            )

    axes[0].set_ylabel("Biomass flux")
    fig.suptitle(plot_title(df), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Input: {in_csv}")
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
