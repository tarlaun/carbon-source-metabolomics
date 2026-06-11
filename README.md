# Phase 2 — Annotated differential metabolites (Starch vs Glucose)

Untargeted LC-MS/MS metabolomics of a bacterium grown on **glucose** vs **starch**.
This bundle contains the GNPS-annotated metabolites that differ between the two media,
with statistics and figures, ready for the genome-pairing step.

## ⚠ Read first — how to interpret these numbers
- **Test:** Student's t-test (equal variance) on log10 relative abundances; Benjamini–Hochberg FDR.
  *Significant* = **FDR < 0.05 AND |log2FC| ≥ 1**.
- **Direction:** `log2FC_Glucose_vs_Starch` > 0 → **higher in glucose**; < 0 → **higher in starch**.
- **Replicates are TECHNICAL** (4 per group = repeated extractions/injections of one culture per condition,
  so biological n = 1 vs 1). The p-values/FDR therefore describe *these specific extracts* and are
  **descriptive / hypothesis-generating, not a biological population inference.**
- **Annotations are putative** GNPS MS² library matches (level 2/3). Every row carries match quality and flags:
  - `MQScore` (0–1, higher = better), `MZErrorPPM` — match confidence.
  - `resolved_name` — cleaned, human-readable name (cryptic catalog/CAS IDs resolved via PubChem).
  - `likely_synthetic_or_contaminant` — **TRUE** = match to a synthetic/drug library or a known LC-MS
    contaminant (surfactant, sulfa drug). Treat these skeptically as bacterial metabolites.
  - `pubchem_cid` / `pubchem_url` / `pubchem_title` / `pubchem_formula_match` — verified PubChem cross-reference
    (only trust the PubChem name when `pubchem_formula_match = TRUE`; all included ones are TRUE).

## ➜ Recommended input for genome pairing
Use **`tables/confident_shortlist.csv`** (MQScore ≥ 0.85, adduct/isotope ions collapsed to one feature per
metabolite). Dropping the rows flagged `likely_synthetic_or_contaminant` leaves **~41 clean, biologically
plausible metabolites (24 higher in glucose / 17 higher in starch)**.

## Files

**tables/**
| file | rows | what |
|---|---|---|
| `annotated_significant_features.csv` | 74 | all annotated & significant features + stats + annotation (`direction` column) |
| `annotated_significant_glucose.csv` | 44 | …split: higher in glucose |
| `annotated_significant_starch.csv` | 30 | …split: higher in starch |
| `confident_shortlist.csv` | 45 | **recommended** — MQScore≥0.85, adducts collapsed, resolved names |
| `confident_shortlist_glucose.csv` | 26 | …split: higher in glucose |
| `confident_shortlist_starch.csv` | 19 | …split: higher in starch |
| `pubchem_resolution.csv` | — | audit: cryptic catalog-ID names → PubChem structures (formula-checked) |

**figures/**
| file | what |
|---|---|
| `Q4_boxplots_grid.png` | **recommended for the paper** — faceted grid, one panel per metabolite, Starch vs Glucose, replicate points, numeric FDR |
| `Q4_boxplots_top10.png` | top-10 per direction (literal: raw annotated, ranked by Student's t) |
| `Q4_boxplots_shortlist.png` | compact horizontal overview of the confident shortlist |

**REPORT.md** — full methods and results (Phase 1 unsupervised analysis + Phase 2 annotation).

## Key columns
`feature_id` (= GNPS `#Scan#`), `mz`, `rt`, `resolved_name`, `Compound_Name`, `direction`,
`log2FC_Glucose_vs_Starch`, `mean_Starch`, `mean_Glucose`, `student_t/p/FDR`, `welch_p/FDR`,
`MQScore`, `MZErrorPPM`, `molecular_formula`, `InChIKey`, `Smiles`, `superclass`/`class`/`subclass`,
`npclassifier_*`, `pubchem_cid/url/title/formula_match`, `likely_synthetic_or_contaminant`.
Shortlist also has `n_collapsed_adducts`, `collapsed_feature_ids`, `isomers_same_annotation`.
