import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit.Chem import PandasTools

PandasTools.molRepresentation = "SVG"

INPUT_FN = "20230307_MAC1_IC50_data.csv"
OUTPUT_FN = "20230307_MAC1_IC50_data.html"
POTENCY_COLS = [
    "MAC1_IC50_GMean_MOD",
    "MAC1_IC50_GMean_NUM (µM)",
    "MAC1_IC50_GMean_STDDEV (×/÷)",
    "MAC1_IC50_GMean_COUNT",
]


def addEnhancedStereoAnnotations(m):
    gpLookup = {
        Chem.StereoGroupType.STEREO_OR: "or",
        Chem.StereoGroupType.STEREO_AND: "&",
        Chem.StereoGroupType.STEREO_ABSOLUTE: "abs",
    }
    sgs = m.GetStereoGroups()
    for i, sg in enumerate(sgs):
        typ = gpLookup[sg.GetGroupType()]
        for at in sg.GetAtoms():
            nt = ""
            if at.HasProp("atomNote"):
                nt += at.GetProp("atomNote") + ","
            nt += f"{typ}{i+1}"
            at.SetProp("atomNote", nt)


### LOAD DATA ###
df = pd.read_csv(INPUT_FN)

### ADD MOL WITH ENHANCED STEREO INFO ###
PandasTools.AddMoleculeColumnToFrame(df, "parent_CXSMILES", "parent_molecule")
for _, row in df.iterrows():
    addEnhancedStereoAnnotations(row["parent_molecule"])

### CLEAN UP DF FOR EXPORT ###
df = df[
    [
        "parent_molecule",
        "ASAP_molecule_ID",
        "ASAP_batch_ID",
        "ASAP_VCID",
        "vendor_catalog_ID",
        "vendor_batch_ID",
        "vendor",
        "parent_CXSMILES",
        "salt_name",
        "salt_ratio",
    ]
    + POTENCY_COLS
]
df.fillna("", inplace=True)

### WRITE OUT HTML FILE ###
html_content = df.to_html()
with open(OUTPUT_FN, "w") as f:
    f.write(html_content)
