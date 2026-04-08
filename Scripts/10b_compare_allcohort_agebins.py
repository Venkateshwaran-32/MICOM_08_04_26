from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCENARIO_DIR = PROJECT_ROOT / "Results" / "fba" / "scenarios"
OUT_DIR = PROJECT_ROOT / "Results" / "fba"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EQUAL_FILE = SCENARIO_DIR / "allcohort_equal_abundance_community_exchange_fluxes_by_diet.csv"
OUT_MASTER = OUT_DIR / "allcohort_agebin_vs_equal_exchange_comparisons.csv"


def main():
    if not EQUAL_FILE.exists():
        raise FileNotFoundError(f"Missing baseline file: {EQUAL_FILE}. Run script 09b first.")

    eq = pd.read_csv(EQUAL_FILE)[["diet", "exchange_id", "flux"]].rename(columns={"flux": "flux_equal"})

    age_files = sorted(SCENARIO_DIR.glob("allcohort_agebin_*_community_exchange_fluxes_by_diet.csv"))
    if not age_files:
        raise RuntimeError("No allcohort age-bin scenario exchange files found. Run script 09b first.")

    all_rows = []
    for fp in age_files:
        sc_name = fp.name.replace("_community_exchange_fluxes_by_diet.csv", "")
        wt = pd.read_csv(fp).copy()

        out = wt.merge(eq, on=["diet", "exchange_id"], how="left")
        out["scenario"] = sc_name
        out["delta_flux_vs_equal"] = out["flux"] - out["flux_equal"]
        out["abs_delta_flux_vs_equal"] = out["delta_flux_vs_equal"].abs()
        out["log2_fc_abs_flux_vs_equal"] = np.log2((out["abs_flux"] + 1e-9) / (out["flux_equal"].abs() + 1e-9))
        out = out.sort_values(["diet", "abs_delta_flux_vs_equal"], ascending=[True, False])

        out_file = SCENARIO_DIR / f"{sc_name}_vs_allcohort_equal_exchange_flux_comparison_by_diet.csv"
        out.to_csv(out_file, index=False)
        print(f"Saved: {out_file} ({len(out)} rows)")

        all_rows.append(out)

    master = pd.concat(all_rows, ignore_index=True)
    master.to_csv(OUT_MASTER, index=False)
    print(f"Saved: {OUT_MASTER} ({len(master)} rows)")


if __name__ == "__main__":
    main()
