#!/usr/bin/env python3
"""
Enrich the annotated/shortlist tables with verified PubChem cross-references.

For every annotated feature (looked up by exact InChIKey):
  pubchem_cid            - PubChem Compound ID
  pubchem_url            - clickable structure page
  pubchem_title          - PubChem's curated compound title
  pubchem_formula_match  - True iff PubChem's MolecularFormula == the annotation's
                           molecular_formula  (ACCURACY GUARD: a name is only trusted
                           when this is True)

For the bare catalog-ID features (CAS / Enamine Z / NCGC / vendor ST) whose PubChem
formula matches, 'resolved_name' is upgraded to the PubChem title and
'name_is_catalog_id' is set False. Nothing is invented: every value comes from
PubChem keyed on the InChIKey, and is discarded if the formula does not match.

Network: PubChem PUG REST via curl (python ssl fails behind the proxy).
A standalone audit file 'pubchem_resolution.csv' documents the catalog-ID rescues.
"""
import os, re, json, subprocess, time
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
TABLES = ["annotated_features_all_stats.csv", "annotated_significant_features.csv",
          "annotated_significant_glucose.csv", "annotated_significant_starch.csv",
          "confident_shortlist.csv", "confident_shortlist_glucose.csv",
          "confident_shortlist_starch.csv"]


def get(url):
    r = subprocess.run(["curl", "-s", "--max-time", "30", url], capture_output=True, text=True)
    try:
        return json.loads(r.stdout)
    except Exception:
        return {}


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def clean_title(t):
    t = str(t).strip()
    m = re.match(r"^[A-Za-z0-9\-]+_C\d[\dA-Za-z]*_(.+)$", t)   # NCGC..._formula_Name
    if m:
        t = m.group(1)
    return t.strip()


def main():
    base = pd.read_csv(os.path.join(TAB, "annotated_features_all_stats.csv"))
    keys = sorted(base.InChIKey.dropna().unique().tolist())
    print(f"querying PubChem for {len(keys)} unique InChIKeys ...")

    ik2cid, ik2pf = {}, {}
    for ch in chunks(keys, 25):
        r = get(f"{BASE}/compound/inchikey/{','.join(ch)}/property/MolecularFormula,InChIKey/JSON")
        for p in r.get("PropertyTable", {}).get("Properties", []):
            ik = p.get("InChIKey")
            if ik:
                ik2cid[ik] = p.get("CID"); ik2pf[ik] = p.get("MolecularFormula")
        time.sleep(0.35)
    print(f"  resolved CIDs for {len(ik2cid)}/{len(keys)} keys")

    cids = sorted({c for c in ik2cid.values() if c})
    cid2title = {}
    for ch in chunks(cids, 25):
        r = get(f"{BASE}/compound/cid/{','.join(map(str, ch))}/description/JSON")
        for info in r.get("InformationList", {}).get("Information", []):
            if info.get("Title") and info.get("CID") not in cid2title:
                cid2title[info["CID"]] = info["Title"]
        time.sleep(0.35)
    # individual fallback for CIDs the batch did not return a Title for
    missing = [c for c in cids if c not in cid2title]
    for c in missing:
        r = get(f"{BASE}/compound/cid/{c}/description/JSON")
        for info in r.get("InformationList", {}).get("Information", []):
            if info.get("Title"):
                cid2title[c] = info["Title"]; break
        time.sleep(0.25)
    print(f"  got titles for {len(cid2title)}/{len(cids)} CIDs ({len(missing)} via fallback)")

    def enrich(df):
        cid = df.InChIKey.map(ik2cid)
        pf = df.InChIKey.map(ik2pf)
        title = cid.map(cid2title)
        match = (pf.notna() & df.molecular_formula.notna()
                 & (pf.astype(str) == df.molecular_formula.astype(str)))
        df = df.copy()
        df["pubchem_cid"] = cid.apply(lambda x: int(x) if pd.notna(x) else pd.NA).astype("Int64")
        df["pubchem_url"] = cid.apply(lambda x: f"https://pubchem.ncbi.nlm.nih.gov/compound/{int(x)}"
                                      if pd.notna(x) else "")
        df["pubchem_title"] = title.fillna("")
        df["pubchem_formula_match"] = match
        # upgrade resolved_name for catalog-ID rows that pass the formula guard
        up = df.name_is_catalog_id & match & title.notna()
        df.loc[up, "resolved_name"] = title[up].map(clean_title)
        df.loc[up, "name_is_catalog_id"] = False
        return df

    for fn in TABLES:
        p = os.path.join(TAB, fn)
        if not os.path.exists(p):
            print("skip:", fn); continue
        df = enrich(pd.read_csv(p))
        df.to_csv(p, index=False)
        print(f"{fn:42}  pubchem_cid={int(df.pubchem_cid.notna().sum()):3}/{len(df):3}  "
              f"formula_match={int(df.pubchem_formula_match.sum()):3}")

    # standalone audit of the catalog-ID rescues
    allt = pd.read_csv(os.path.join(TAB, "annotated_features_all_stats.csv"))
    audit = allt[allt.pubchem_title.astype(str).str.len() > 0][
        ["feature_id", "Compound_Name", "resolved_name", "molecular_formula",
         "pubchem_cid", "pubchem_formula_match", "pubchem_title", "pubchem_url"]]
    # focus the audit on the originally-cryptic ones
    cryptic = allt["Compound_Name"].str.contains(
        r"^(?:Z\d{6,}|NCGC|ST\d{5,}|\d{2,7}-\d{2}-\d)", regex=True, na=False)
    allt[cryptic][["feature_id", "Compound_Name", "resolved_name", "molecular_formula",
                   "pubchem_cid", "pubchem_formula_match", "pubchem_title", "pubchem_url"]
                  ].drop_duplicates("Compound_Name").to_csv(
        os.path.join(TAB, "pubchem_resolution.csv"), index=False)
    print("\nwrote pubchem_resolution.csv (catalog-ID audit)")


if __name__ == "__main__":
    main()
