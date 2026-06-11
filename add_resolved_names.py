#!/usr/bin/env python3
"""
Enrich the annotated / shortlist tables with three columns that make the GNPS
'Compound_Name' field usable for the genome-pairing step:

  resolved_name                 - a cleaned, human-readable name: library prefixes
                                  (Massbank:/ReSpect:/'Spectral Match to ... from'),
                                  collision-energy/eV suffixes, trailing adduct tags
                                  and 'ID_formula_Name' wrappers are stripped; a small
                                  curated map resolves confirmed registry IDs.
  name_is_catalog_id            - True when the name is still a bare registry/catalog
                                  code (CAS / Enamine Z / NCGC / vendor ST / ZINC /
                                  CHEMBL). For these, trust molecular_formula / InChIKey
                                  / Smiles, NOT the name.
  likely_synthetic_or_contaminant - True for matches to synthetic screening / drug
                                  libraries (Enamine, ZINC, CHEMBL, MCE-DRUG) or to
                                  well-known LC-MS contaminants/drugs. Heuristic flag:
                                  treat such hits skeptically as bacterial metabolites.

Confirmed by mass+formula+InChIKey: CAS 83-88-5 = riboflavin (C17H20N4O6, [M+H]+ 377.15).
Columns are inserted right after Compound_Name; tables are rewritten in place.
"""
import os, re
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANN = os.path.join(ROOT, "e1fe0f44e01f4fd0a8deef0718c20c3d-merged_results_with_gnps.tsv")
TAB = os.path.join(ROOT, "results", "tables")

TABLES = ["annotated_features_all_stats.csv", "annotated_significant_features.csv",
          "annotated_significant_glucose.csv", "annotated_significant_starch.csv",
          "confident_shortlist.csv", "confident_shortlist_glucose.csv",
          "confident_shortlist_starch.csv"]

CURATED_ID = {"83-88-5": "Riboflavin (vitamin B2)"}

BARE_ID = re.compile(r"^(Z\d{6,}|NCGC[\dA-Za-z\-]*|ST\d{5,}|ZINC\d+|S?CHEMBL\d+|\d{2,7}-\d{2}-\d)$")
# synthetic SCREENING-library catalog IDs. NOTE: vendor stock codes like 'ST…'
# (TimTec) are intentionally NOT here — vendors also sell natural products, so an
# ST id alone is not evidence of a synthetic compound (e.g. ST50331689 = an aurone).
SYN_VENDOR = re.compile(r"^(Z\d{6,}|ZINC\d+|S?CHEMBL\d+)")
DRUG_LIBS = {"MCE-DRUG.mgf"}
CONTAM_KW = ["sulfadimethoxine", "sulfamethazine", "sulfamethoxazole", "sulfadiazine",
             "sulfaquinoxaline", "amitriptyline", "polidocanol", "lauryldiethanolamine",
             "diethanolamine", "palmitamide", "oleamide", "erucamide", "stearamide",
             "phthalate", "polyethylene glycol", "tributyl", "triton", "tween"]


def resolve_name(name):
    s = str(name).strip().strip('"').strip()
    m = re.match(r"^Spectral Match to (.+?) from .+$", s)
    if m:
        s = m.group(1)
    s = re.sub(r"^(?:Massbank|MassBank|ReSpect|MoNA):\S*\s*", "", s).strip()   # accession prefix
    if "!" in s:                                                              # ID!name
        s = s.split("!", 1)[1]
    s = re.sub(r"\s*\[IIN-based:[^\]]*\]", "", s)
    m = re.match(r"^[A-Za-z0-9\-]+_C\d[\dA-Za-z]*_(.+)$", s)                   # ID_formula_Name
    if m:
        s = m.group(1)
    s = re.split(r"\s+CollisionEnergy", s)[0]
    s = re.split(r"\s+-\s+[\d.]+\s*eV", s)[0]
    s = re.sub(r"\s*\[M[\+\-][^\]]*\]\+?$", "", s)                            # trailing adduct
    s = s.split("|")[0].strip().strip('"').strip()
    key = re.split(r"\s", s)[0] if s else s
    if key in CURATED_ID:
        s = CURATED_ID[key]
    return s if s else str(name)


def flags(orig_name, resolved, libname):
    is_id = bool(BARE_ID.match(resolved.strip()))
    n = str(orig_name).lower()
    syn = (bool(SYN_VENDOR.match(str(orig_name).strip()))
           or (str(libname) in DRUG_LIBS)
           or any(kw in n for kw in CONTAM_KW))
    return is_id, syn


def main():
    ann = pd.read_csv(ANN, sep="\t")
    lib = ann.set_index("#Scan#")["LibraryName"].to_dict()

    for fn in TABLES:
        path = os.path.join(TAB, fn)
        if not os.path.exists(path):
            print("skip (missing):", fn); continue
        df = pd.read_csv(path)
        df["resolved_name"] = df["Compound_Name"].map(resolve_name)
        res = df.apply(lambda r: flags(r["Compound_Name"], r["resolved_name"],
                                       lib.get(int(r["feature_id"]))), axis=1)
        df["name_is_catalog_id"] = [a for a, _ in res]
        df["likely_synthetic_or_contaminant"] = [b for _, b in res]
        # reorder: put the 3 new cols right after Compound_Name
        cols = list(df.columns)
        for c in ["resolved_name", "name_is_catalog_id", "likely_synthetic_or_contaminant"]:
            cols.remove(c)
        i = cols.index("Compound_Name") + 1
        cols[i:i] = ["resolved_name", "name_is_catalog_id", "likely_synthetic_or_contaminant"]
        df = df[cols]
        df.to_csv(path, index=False)
        print(f"{fn:42}  n={len(df):3}  catalog_id={int(df.name_is_catalog_id.sum()):2}  "
              f"synthetic/contam={int(df.likely_synthetic_or_contaminant.sum()):2}")


if __name__ == "__main__":
    main()
