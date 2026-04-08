from pathlib import Path          # lets us work with file/folder paths safely (Windows-friendly)
import pandas as pd               # for making tables (DataFrames) and saving CSV
import cobra                      # COBRApy: loads SBML models and runs FBA

# Resolve paths from the project root so the script works across machines/OSes.
PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Folder where your downloaded VMH SBML (.xml) files are stored
MODELS_DIR = PROJECT_ROOT / "Models" / "vmh_agora_sbml"

# Where we will save the QC summary CSV
OUT_CSV = PROJECT_ROOT / "Results" / "qc" / "model_qc_summary.csv"

# --- Basic debugging prints so you can see what's happening ---
print("MODELS_DIR =", MODELS_DIR)            # shows the folder path your script is looking at
print("Exists?    =", MODELS_DIR.exists())   # checks if that folder actually exists

# Find all .xml files in the models folder (your SBML files)
files = sorted(MODELS_DIR.glob("*.xml"))     # glob("*.xml") = match all XML files
print("XML files found:", len(files))        # how many XML files were found
print("First few:", [f.name for f in files[:5]])  # show first 5 filenames (quick sanity check)

# Make sure the output folder Results/qc exists (creates it if it doesn't)
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

rows = []  # we will store one "row" of QC output per model here

# --- Loop through each SBML model file and try to load it ---
for fp in files:   
    try:
        m = cobra.io.read_sbml_model(str(fp))      # load the SBML model into COBRApy

        # Collect simple QC metrics ("model size") + objective info
        rows.append({
            "file": fp.name,                       # filename of the SBML
            "model_id": m.id,                      # internal model ID inside SBML
            "reactions": len(m.reactions),         # number of reactions
            "metabolites": len(m.metabolites),     # number of metabolites
            "genes": len(m.genes),                 # number of genes in GPR rules
            "objective": str(m.objective.expression),  # what the model is optimizing (usually biomass)
        })

        print("OK  ", fp.name)                     # tells you this model loaded successfully

    except Exception as e:
        # If loading fails, store the error instead (so you know which model broke)
        rows.append({
            "file": fp.name,
            "error": repr(e)                           # error message from Python/COBRApy
        })

        print("BAD ", fp.name, repr(e))            # tells you this model failed to load

# Convert all collected rows into a table
df = pd.DataFrame(rows)

# Save the QC table to a CSV so you can send to PI / keep records
df.to_csv(OUT_CSV, index=False)

print("\nSaved QC CSV to:", OUT_CSV)               # confirm where it saved
print(df.head(20))                                 # print the first 20 rows in terminal
