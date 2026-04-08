from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FBA_DIR = PROJECT_ROOT / "Results" / "fba"
OUT_DIR = FBA_DIR / "scenarios"
OUT_DIR.mkdir(parents=True, exist_ok=True)

IN_SUMMARY = FBA_DIR / "community_weighting_scenarios_summary.csv"
IN_SPECIES = FBA_DIR / "community_weighting_scenarios_species_biomass.csv"
IN_EXCHANGE = FBA_DIR / "community_weighting_scenarios_exchange_fluxes.csv"
IN_COMPARE = FBA_DIR / "community_weighting_vs_equal_exchange_comparison.csv"


def save_if_nonempty(df: pd.DataFrame, path: Path):
    if df.empty:
        return
    df.to_csv(path, index=False)
    print(f"Saved: {path}")


def main():
    for p in [IN_SUMMARY, IN_SPECIES, IN_EXCHANGE, IN_COMPARE]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required file: {p}")

    df_summary = pd.read_csv(IN_SUMMARY)
    df_species = pd.read_csv(IN_SPECIES)
    df_exchange = pd.read_csv(IN_EXCHANGE)
    df_compare = pd.read_csv(IN_COMPARE)

    # SG90 median (main requested scenario)
    save_if_nonempty(
        df_summary[df_summary["scenario"] == "sg90_median"],
        OUT_DIR / "sg90_median_community_summary_by_diet.csv",
    )
    save_if_nonempty(
        df_species[df_species["scenario"] == "sg90_median"],
        OUT_DIR / "sg90_median_species_biomass_by_diet.csv",
    )
    save_if_nonempty(
        df_exchange[df_exchange["scenario"] == "sg90_median"],
        OUT_DIR / "sg90_median_community_exchange_fluxes_by_diet.csv",
    )

    # Equal abundance baseline
    save_if_nonempty(
        df_summary[df_summary["scenario"] == "equal_abundance"],
        OUT_DIR / "equal_abundance_community_summary_by_diet.csv",
    )
    save_if_nonempty(
        df_species[df_species["scenario"] == "equal_abundance"],
        OUT_DIR / "equal_abundance_species_biomass_by_diet.csv",
    )
    save_if_nonempty(
        df_exchange[df_exchange["scenario"] == "equal_abundance"],
        OUT_DIR / "equal_abundance_community_exchange_fluxes_by_diet.csv",
    )

    # SG90 age-group scenarios (if present)
    age_sc = sorted([s for s in df_summary["scenario"].astype(str).unique() if s.startswith("sg90_age_")])
    for scenario in age_sc:
        tag = scenario.replace("sg90_age_", "")
        save_if_nonempty(
            df_summary[df_summary["scenario"] == scenario],
            OUT_DIR / f"sg90_age_{tag}_community_summary_by_diet.csv",
        )
        save_if_nonempty(
            df_species[df_species["scenario"] == scenario],
            OUT_DIR / f"sg90_age_{tag}_species_biomass_by_diet.csv",
        )
        save_if_nonempty(
            df_exchange[df_exchange["scenario"] == scenario],
            OUT_DIR / f"sg90_age_{tag}_community_exchange_fluxes_by_diet.csv",
        )

    # Direct comparison tables
    save_if_nonempty(
        df_compare[df_compare["scenario"] == "sg90_median"],
        OUT_DIR / "sg90_median_vs_equal_exchange_flux_comparison_by_diet.csv",
    )
    for scenario in age_sc:
        tag = scenario.replace("sg90_age_", "")
        save_if_nonempty(
            df_compare[df_compare["scenario"] == scenario],
            OUT_DIR / f"sg90_age_{tag}_vs_equal_exchange_flux_comparison_by_diet.csv",
        )

    print(f"\nScenario-specific files are in: {OUT_DIR}")


if __name__ == "__main__":
    main()
