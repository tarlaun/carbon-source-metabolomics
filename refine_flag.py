#!/usr/bin/env python3
"""
Recompute only the `likely_synthetic_or_contaminant` column in the enriched tables,
using the refined rule (vendor 'ST…' stock codes no longer trigger it — they are not
evidence of a synthetic compound). resolved_name / name_is_catalog_id / pubchem_* are
left untouched, so the verified PubChem resolutions are preserved (no re-query).
"""
import os, re
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
ANN = os.path.join(ROOT, "e1fe0f44e01f4fd0a8deef0718c20c3d-merged_results_with_gnps.tsv")
TABLES = ["annotated_features_all_stats.csv", "annotated_significant_features.csv",
          "annotated_significant_glucose.csv", "annotated_significant_starch.csv",
          "confident_shortlist.csv", "confident_shortlist_glucose.csv",
          "confident_shortlist_starch.csv"]

SYN_VENDOR = re.compile(r"^(Z\d{6,}|ZINC\d+|S?CHEMBL\d+)")   # ST… removed
DRUG_LIBS = {"MCE-DRUG.mgf"}
CONTAM_KW = ["sulfadimethoxine", "sulfamethazine", "sulfamethoxazole", "sulfadiazine",
             "sulfaquinoxaline", "amitriptyline", "polidocanol", "lauryldiethanolamine",
             "diethanolamine", "palmitamide", "oleamide", "erucamide", "stearamide",
             "phthalate", "polyethylene glycol", "tributyl", "triton", "tween"]


def main():
    lib = pd.read_csv(ANN, sep="\t").set_index("#Scan#")["LibraryName"].to_dict()

    def flag(name, fid):
        n = str(name).lower()
        return (bool(SYN_VENDOR.match(str(name).strip()))
                or (str(lib.get(int(fid))) in DRUG_LIBS)
                or any(k in n for k in CONTAM_KW))

    for fn in TABLES:
        p = os.path.join(TAB, fn)
        if not os.path.exists(p):
            print("skip:", fn); continue
        df = pd.read_csv(p)
        df["likely_synthetic_or_contaminant"] = df.apply(
            lambda r: flag(r["Compound_Name"], r["feature_id"]), axis=1)
        df.to_csv(p, index=False)
        print(f"{fn:42}  flagged={int(df.likely_synthetic_or_contaminant.sum()):2}/{len(df)}")


if __name__ == "__main__":
    main()
