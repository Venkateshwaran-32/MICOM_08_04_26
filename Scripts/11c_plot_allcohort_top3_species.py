from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCENARIO_DIR = PROJECT_ROOT / "Results" / "fba" / "scenarios"
OUT_CSV = PROJECT_ROOT / "Results" / "fba" / "allcohort_agebin_top3_species_by_diet.csv"
OUT_FIG = PROJECT_ROOT / "Results" / "figures" / "allcohort_agebin_top3_species_by_diet.png"

AGE_ORDER = ["21_40", "41_60", "61_70", "71_80", "81_90"]
DIET_ORDER = ["western", "high_fiber"]
FLUX_EPS = 1e-9


def age_from_scenario(name: str) -> str:
    m = re.search(r"allcohort_agebin_(\d+_\d+)", name)
    return m.group(1) if m else "unknown"


def short_species_name(full_name: str) -> str:
    # Convert AGORA model IDs to short labels like "Escherichia coli".
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


def main():
    files = sorted(SCENARIO_DIR.glob("allcohort_agebin_*_species_biomass_by_diet.csv"))
    if not files:
        raise RuntimeError("No allcohort age-bin species biomass files found. Run 09b first.")

    rows = []
    for fp in files:
        scenario = fp.name.replace("_species_biomass_by_diet.csv", "")
        age_group = age_from_scenario(scenario)
        df = pd.read_csv(fp)
        for diet in DIET_ORDER:
            sub = df[df["diet"] == diet].copy()
            if sub.empty:
                continue
            top = sub.sort_values("biomass_flux", ascending=False).head(3).copy()
            top["rank"] = range(1, len(top) + 1)
            top["scenario"] = scenario
            top["age_group"] = age_group
            top["species_short"] = top["species"].apply(short_species_name)
            rows.append(top[["scenario", "age_group", "diet", "rank", "species", "species_short", "biomass_flux", "weight"]])

    out = pd.concat(rows, ignore_index=True)
    out["age_group"] = pd.Categorical(out["age_group"], categories=AGE_ORDER, ordered=True)
    out = out.sort_values(["age_group", "diet", "rank"])
    out.to_csv(OUT_CSV, index=False)

    # Also save strict rank table (including zeros) so "top 3" is explicit.
    out_rank = out.copy()
    out_rank_file = OUT_CSV.with_name("allcohort_agebin_top3_species_by_diet_ranked_including_zeros.csv")
    out_rank.to_csv(out_rank_file, index=False)

    # Cleaner view: one panel per diet, stacked by short species names.
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5), sharey=True)
    for j, diet in enumerate(DIET_ORDER):
        ax = axes[j]
        sub = out[(out["diet"] == diet) & (out["biomass_flux"] > FLUX_EPS)].copy()
        if sub.empty:
            ax.set_title(f"{diet} (no positive growth rows)")
            ax.axis("off")
            continue

        pivot = (
            sub.groupby(["age_group", "species_short"], as_index=False)["biomass_flux"]
            .sum()
            .pivot(index="age_group", columns="species_short", values="biomass_flux")
            .fillna(0.0)
            .reindex(AGE_ORDER)
            .fillna(0.0)
        )
        pivot.plot(kind="bar", stacked=True, ax=ax, width=0.75)
        ax.set_title(diet)
        ax.set_xlabel("age bin")
        ax.set_ylabel("biomass_flux")
        ax.tick_params(axis="x", rotation=0)
        ax.legend(title="species", fontsize=8, title_fontsize=9, loc="upper right")

    fig.suptitle("Dominant Positive-Growth Species by Age Bin and Diet", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_FIG, dpi=220)
    plt.close(fig)

    print(f"Saved: {OUT_CSV} ({len(out)} rows)")
    print(f"Saved: {out_rank_file} ({len(out_rank)} rows)")
    print(f"Saved: {OUT_FIG}")


if __name__ == "__main__":
    main()
