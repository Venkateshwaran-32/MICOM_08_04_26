from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = PROJECT_ROOT / "Results" / "fba" / "community_species_biomass_flux_full_access.csv"


def main() -> None:
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"CSV not found: {CSV_PATH}")

    df = pd.read_csv(CSV_PATH)

    print(f"Loaded: {CSV_PATH}")
    print(f"Shape: {df.shape}")
    print()
    print(df.head())
    print()
    print(df.dtypes)


if __name__ == "__main__":
    main()
