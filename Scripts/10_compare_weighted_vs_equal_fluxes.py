from pathlib import Path
import pandas as pd
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FBA_DIR = PROJECT_ROOT / "Results" / "fba"
SCENARIO_DIR = FBA_DIR / "scenarios"

IN_EXCHANGE = SCENARIO_DIR / "sg90_median_community_exchange_fluxes_by_diet.csv"
IN_EQUAL = SCENARIO_DIR / "equal_abundance_community_exchange_fluxes_by_diet.csv"
OUT_COMPARE = SCENARIO_DIR / "sg90_median_vs_equal_exchange_flux_comparison_by_diet.csv"


def main():
    if not IN_EXCHANGE.exists() or not IN_EQUAL.exists():
        raise FileNotFoundError(
            f"Missing scenario inputs: {IN_EXCHANGE} and/or {IN_EQUAL}. Run script 09 first."
        )

    wt = pd.read_csv(IN_EXCHANGE).copy()
    eq = pd.read_csv(IN_EQUAL)[["diet", "exchange_id", "flux"]].rename(columns={"flux": "flux_equal"}).copy()

    out = wt.merge(eq, on=["diet", "exchange_id"], how="left")
    out["scenario"] = "sg90_median"
    out["delta_flux_vs_equal"] = out["flux"] - out["flux_equal"]
    out["abs_delta_flux_vs_equal"] = out["delta_flux_vs_equal"].abs()
    out["log2_fc_abs_flux_vs_equal"] = np.log2((out["abs_flux"] + 1e-9) / (out["flux_equal"].abs() + 1e-9))
    out = out.sort_values(["scenario", "diet", "abs_delta_flux_vs_equal"], ascending=[True, True, False])
    out.to_csv(OUT_COMPARE, index=False)

    print(f"Saved: {OUT_COMPARE} ({len(out)} rows)")


if __name__ == "__main__":
    main()
