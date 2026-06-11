# Phase 1 — Metabolite profiling: Starch vs Glucose

**Driving question:** Does the bacterium's metabolite profile differ when grown on **starch** (`EXT_*`)
versus **glucose** (`aG_*`)? If so, are the data normal (→ parametric vs non-parametric), and which
features differ?

**Input:** `2026-06-06T21-12_export.csv` — already blank-subtracted, imputed, TIC-normalised.
**Design:** 8 samples — 4 Starch + 4 Glucose — × **7 032 features** (`id_mz_rt&NA`; the leading `id`
is the unique feature identifier for later molecular look-up). Every row sums to **1.0** → values are
**relative abundances (compositional)**; no zeros/NaN/negatives (imputed floor ≈ 1.5×10⁻¹⁰; only **9 of
7 032 features** ever touch it).

> Every headline number below was **independently re-derived from the raw CSV by separate agents writing
> their own code** (PCA, exact PERMANOVA, PCoA, BH-FDR, Mann–Whitney) and reproduced to full precision.
> The methodology was also adversarially reviewed; the resulting caveats are folded in here.

---

## TL;DR

| Question | Answer |
|---|---|
| **Q1 — Do the profiles differ?** | **Yes, decisively** (descriptively). PCA PC1 = **67 %** and fully separates the media; PCoA1 (Bray–Curtis) = **90 %**; PERMANOVA F = **52.1**, R² = **0.90** (Aitchison/CLR R² = **0.68**), p = **0.029** — the *smallest value the 4-vs-4 design can return*. PERMDISP is **n.s.** (separation is a location shift, not a spread artifact); silhouette = **0.75**; clustering recovers the groups perfectly; dropping the within-group outlier leaves R² = 0.93. |
| **Q2 — Are the data normal?** | **Untestable at n = 4/group, so don't depend on it.** Shapiro is essentially powerless at n=4 (it "passes" ~80 % of features by default); log10 reduces skew but the bulk distribution is bimodal. **Decision:** work on log10 and use a **variance-moderated (empirical-Bayes) test** as primary, cross-checked against ordinary Welch t — *not* a per-feature normality screen. |
| **Q3 — Which features differ?** | A **majority of the metabolome shifts** — moderated-t FDR<0.05: **4 644**; Welch+BH: **4 289**; Welch+BY (dependence-robust): **2 469**. That ~60 % rate reflects the strong *global* glucose-vs-starch shift, so it is **not an actionable shortlist**. A stricter, robust shortlist (FDR<0.05 **and** \|log2FC\|≥1 **and** not floor-driven **and** detected **and** low within-group CV) gives **538** features. See `Q3_volcano.png`. |

> ⚠️ **Before trusting any p-value, read [Caveat 1](#1-biologicalvs-technical-replication--most-important).**

---

## Methods (one consistent preprocessing for every step)

* **`R`** = relative abundance (table as given, rows sum to 1). **`L = log10(R)`** for symmetric-scale work.
* **PCA**: `L` → **Pareto scaling** (mean-centre, ÷√SD), the metabolomics convention.
* **Distances**: **Bray–Curtis** on `R` (standard for relative abundances) for PCoA/PERMANOVA;
  **Aitchison** (Euclidean on the centred-log-ratio, CLR) as the compositional cross-check; Euclidean(log10) also checked.
* **PERMANOVA / PERMDISP**: exact tests over all C(8,4)=70 label permutations.
* **Differential**: per-feature fold change = log2 of the **geometric-mean ratio** (Glucose/Starch),
  i.e. the difference of mean-log10s — consistent with testing on logs. Tests: **moderated t**
  (limma-style empirical-Bayes variance shrinkage, implemented by hand) as primary; ordinary **Welch t**
  for comparison; **Mann–Whitney U** as a separation flag. Multiple testing: **Benjamini–Hochberg** (primary)
  and **Benjamini–Yekutieli** (dependence-robust sensitivity).
* Pipelines: `analysis/metabolite_analysis.py` (primary) + `analysis/robustness_checks.py` (sensitivity).
  PERMANOVA/PERMDISP/PCoA/BH/BY/moderated-t coded from scratch (`scikit-bio`/`statsmodels` unavailable).

---

## Q1 — The two profiles are clearly different

| Evidence | Result |
|---|---|
| **PCA** (`Q1_PCA.png`) | PC1 = **67.1 %**, PC2 = 8.6 %; 4 Starch at one PC1 extreme, 4 Glucose at the other — complete separation. CLR-based PCA gives an essentially identical PC1 (67.1 %), so closure is not driving it. |
| **PCoA / Bray–Curtis** (`Q1_PCoA_BrayCurtis.png`) | PCoA1 = **89.7 %**; clean split. |
| **PERMANOVA (exact)** | Bray–Curtis: F=**52.1**, R²=**0.897**; Aitchison/CLR: F=**12.5**, R²=**0.676**; both p=**0.0286**. |
| **PERMDISP (exact)** | F=0.28, **p=0.51 (n.s.)** → within-group dispersions equal (0.10 vs 0.12); the signal is a centroid shift, not a spread difference. |
| **Silhouette (Bray–Curtis)** | **0.75** (strong, well-separated clusters). |
| **Hierarchical clustering** (`Q1_dendrograms.png`, `Q1_heatmap_top50.png`) | Ward/Euclidean and average/Bray–Curtis both split perfectly into the two media. |
| **Leave-one-out** (drop outlier `EXT_wo_ami_1`) | PERMANOVA F=**68.0**, R²=**0.93**, p=0.0286 — separation is not driven by the outlier. |

**On the p-value:** with 4 vs 4 there are only C(8,4)=70 label permutations, and a partition and its complement
give the identical pseudo-F, so the **smallest attainable p is 2/70 = 0.0286** — and the true labelling is the
single most-separated of all 70 (next-best F ≈ 2). The test is therefore **saturated**: the evidence for Q1 is
the **effect size** (R²=0.90 Bray–Curtis / 0.68 Aitchison) and the visual separation, *not* a precise p-value.
Reporting the more conservative **Aitchison R²=0.68** avoids overselling, since Bray–Curtis is weighted toward a
few high-abundance ions.

**→ Yes, the profiles differ. Proceed to Q2/Q3 (subject to Caveat 1).**

## Q2 — Normality: do not gate on it (n = 4)

* Per-feature normality is **statistically untestable at n = 4**: Shapiro–Wilk has almost no power (it "passes"
  ~80 % of features here, and in simulation rejects a true exponential only ~13 % of the time), and standardised
  residuals from 4 points are mechanically bounded — the flat QQ tails in `Q2_normality.png` are that artifact, not data.
* The bulk log10-intensity distribution is **bimodal** (near-floor mode + signal mode); log10 reduces per-feature
  skew (median 0.14 → −0.07) but a single "normal/log-normal" label oversimplifies.
* **Decision rule (stated explicitly):** *n = 4/group makes a formal normality screen uninformative, so we do not
  use one to choose the test.* We instead work on the log10 scale and adopt a **variance-moderated empirical-Bayes
  test** (the standard fix for small-n omics, which borrows variance information across the 7 032 features to
  stabilise the ~3-df per-feature estimates), cross-checked against ordinary Welch t. Results are concordant.

## Q3 — Differential features + volcano

Direction: **log2FC > 0 ⇒ higher in Glucose**; **< 0 ⇒ higher in Starch** (geometric-mean ratio).

| Test | FDR<0.05 | Note |
|---|---|---|
| **Moderated t (empirical-Bayes)** — *primary ranking* | **4 644** | prior df = 1.1, posterior df = 7.1; stabilises Welch's unstable ~3-df variances |
| Welch t + BH | 4 289 | with \|log2FC\|≥1: 3 999 (↑Glucose 2 374 / ↑Starch 1 625) |
| Welch t + **BY** (dependence-robust) | 2 469 | features are correlated under closure; BY is the conservative bound |
| Mann–Whitney + BH | 4 720 | **degenerate** — see below |

* **The two parametric tests agree** (and agree with the CLR transform: significant-set Jaccard **0.87**, p-value
  Spearman **0.90**), so the differential signal is robust to transform and to the variance model.
* **Mann–Whitney is degenerate at n=4 and is *not* independent corroboration.** Its smallest possible two-sided p
  is 2/70 = 0.0286, reached only at complete group separation. All 4 720 "significant" features are tied at exactly
  that value, so the test provides **zero ability to rank** them — treat it only as a binary "completely separated" flag.
* **Why ~60 % of features are flagged:** this reflects the strong *global* metabolome shift between the two media
  (amplified by compositional closure), not 4 000+ independent biological discoveries. **Do not treat the FDR list
  as a shortlist.** Rank by effect size and reproducibility: the **robust shortlist** (`shortlist_robust_features.csv`,
  **538 features** — FDR<0.05 & \|log2FC\|≥1 & not floor-driven & mean ≥ 1e-5 & within-group CV ≤ 0.5) is the
  actionable set.
* **Volcano** (`Q3_volcano.png`): classic two-wing shape, gap at FC=0, top-12 features labelled by `feature_id`.

Outputs: `differential_features.csv` (Welch + MW, sorted by p) and the richer
**`differential_features_annotated.csv`** (adds moderated-t, BY-FDR, within-group CV, `near_imputation_floor`
flag — use this one for feature look-up).

---

## Critical caveats

### 1. Replication is TECHNICAL — **confirmed by the user**
The four samples per group are **technical/analytical replicates of a single culture per condition** (Starch
`EXT_wo_ami_1, _1, _2, _3`; Glucose `aG_LLE, _1, _2, _3`), confirmed by the user. This was also evident in the
data: a nested Bray–Curtis structure (each base sample farther from its own `_1/_2/_3` triplet than the triplet
members are from each other) and high within-group correlation (r ≈ 0.6–0.8). **Consequence: the biological n is
1 vs 1.** Every p-value/FDR therefore describes whether *these two specific extracts* differ — it quantifies the
reproducibility of the separation, **not** a generalizable carbon-source effect on the organism. All results here
are **descriptive / hypothesis-generating**; the effect sizes (fold change) and ordination are the meaningful
deliverables. To make a population-level claim, repeat with independent biological replicates (≥5/group).

### 2. Tiny-n statistical resolution
PERMANOVA p is pinned at its 0.0286 floor (so effect size, not p, carries Q1); per-feature Welch t has only ~3 df
and produces parametric p-values (e.g. 1e-12) far below what an exact 4-vs-4 test could justify — hence the
moderated-t and the "rank by effect size" guidance above. With proper biological replication, repeat with ≥5/group.

### 3. Compositional constraint (checked, empirically minor here)
Values sum to 1, so a rise in some features forces apparent drops in others; raw-proportion fold-change conflates a
feature's own change with this closure effect. We verified it is **not materially distorting conclusions** here
(matrix is near-complete; CLR reproduces PC1, PERMANOVA structure, and 87 % of the differential hits), but
fold-change *directions* for a minority of features can flip under CLR — read up/down calls as relative-to-total.

### 4. Feature identities are unannotated
`feature_id` is an **unannotated m/z@RT feature** (annotation level 4). One metabolite can appear as several IDs
(adducts/isotopes/in-source fragments). Before counting "molecules", de-duplicate/feature-group, and treat any
m/z-based identity as **putative** pending MS² / authentic-standard confirmation.

---

---

# Phase 2 — Annotation cross-referencing (GNPS) + boxplots

**Goal:** intersect the GNPS library annotations with the significantly different features (split by which
medium they are more prevalent in), build an annotated stats table, and box-plot the top hits.

**Test used (per user's stated assumption of normality + equal variance):** **Student's t-test** (`equal_var=True`)
on log10 relative abundances, BH-FDR. *Significant* = FDR<0.05 & |log2FC|≥1. This agrees with the Phase-1
Welch-based volcano (significant-set Jaccard = **0.94**), so the equal-variance assumption does not change the set.

**Merge:** annotation `#Scan#` ↔ leading feature `id`. Validated chemically — all 174 matched features agree on
m/z to within 0.0005 Da.

| | count |
|---|---|
| Annotated features (of 640 GNPS hits) present in the cleaned dataset | **174** |
| Annotated **and** significant | **74** |
| → higher in **Glucose** (log2FC > 0) | **44** |
| → higher in **Starch** (log2FC < 0) | **30** |

**Top hits** (ranked by Student's t):
* *Higher in glucose:* dodecanoic (lauric) acid, isoflavones (prunetin, biochanin A, afrormosin), a triterpenoid
  saponin, oleanolic acid, indole-3-acetyl-tryptophan.
* *Higher in starch:* diketopiperazines **cyclo(Leu-Pro) / cyclo(Pro-Leu)** (classic bacterial metabolites),
  3-hydroxypropionyl tryptamine, genistein, an Ile-Tyr dipeptide.

### Confident biological shortlist (MQScore filter + adduct de-duplication)

Starting from the 74 annotated + significant features:
1. **Confidence filter** — keep GNPS `MQScore ≥ 0.85` → 52 features.
2. **Collapse adducts/isotopes** — features sharing an InChIKey skeleton **and co-eluting** (|ΔRT| ≤ 0.2 min) are
   one metabolite; keep the best-MQScore representative (e.g. biliverdin's two ions, 2-methylthiopurine's three).
3. **Keep isomers separate, flagged** — same annotation at a *different* RT means a different isomer / degenerate
   library match (frequent for diketopiperazines: cyclo(Leu-Pro)/(Pro-Leu)/(Phe-Pro)), so these are retained as
   distinct rows and marked `*`, not merged.

**Result: 45 confident metabolites — 26 higher in glucose, 19 higher in starch** (11 carry the isomeric-ambiguity
flag; 7 annotations appear in both directions across isomers). See `confident_shortlist.csv`
(+`_glucose`/`_starch` splits) with `n_collapsed_adducts`, `collapsed_feature_ids`, and the ambiguity flags.

**Name resolution + provenance flags** (added to every annotated/shortlist table): `resolved_name` gives a clean
human-readable name (library prefixes, collision-energy tags and `ID_formula_Name` wrappers stripped).
`name_is_catalog_id` marked names that were still a bare registry code; `likely_synthetic_or_contaminant` flags
matches to synthetic/drug libraries (Enamine, MCE-DRUG…) or known LC-MS contaminants (surfactants, sulfa drugs).
**Dropping all flagged rows leaves 41 clean, biologically plausible metabolites (24 glucose / 17 starch)** — the
recommended input for genome pairing.

**Verified PubChem cross-references** (`analysis/pubchem_enrich.py`): every feature was looked up in PubChem by its
exact `InChIKey`, adding `pubchem_cid`, `pubchem_url` (clickable structure page), `pubchem_title`, and
`pubchem_formula_match`. **Accuracy guard:** a PubChem name is only trusted when its molecular formula matches the
annotation's — 126/174 features got a CID and **all 126 matched (0 mismatches)**; the other 48 InChIKeys are simply
absent from PubChem and were left unresolved (nothing invented). This resolved **all** bare catalog IDs to real
structures — e.g. CAS `83-88-5` → riboflavin; `5654-86-4` → cyclo(Leu-Pro) (3-isobutyl-hexahydropyrrolo[1,2-a]
pyrazine-1,4-dione); `82597-82-8` → cyclo(Trp-Phe); `ST50331689` → a hydroxy-aurone; `Z445877824` → an adamantyl
urea (synthetic, flagged). The catalog-ID rescues are itemised in `pubchem_resolution.csv`.

Notable confident hits — *glucose:* dodecanoic acid, isoflavones (prunetin/biochanin A/afrormosin), a triterpenoid
saponin, biliverdin, gramine, the bacterial peptide **fMLF** (N-formyl-Met-Leu-Phe); *starch:* 3-hydroxypropionyl
tryptamine, diketopiperazines, khasianine, 2-methylthiopurine-6,8-diol.

**Boxplots** (all use `resolved_name`; Starch vs Glucose with the 4 replicate points overlaid; each panel labelled
with the **numeric BH-FDR** — not significance stars, since at n=4 technical replicates the p-values are descriptive,
not biological inference, and an exact 4-vs-4 test cannot resolve p below ≈0.03). Generated by
`analysis/make_boxplots.py` (reads the enriched tables, does not rebuild them):
* `Q4_boxplots_grid.png` — **publication-standard format**: a faceted grid, one panel per metabolite (groups on
  x, abundance on log-y, points overlaid), top-10 per direction from the confident, non-flagged set. This is the
  conventional metabolomics figure (MetaboAnalyst-style). **Use this for the paper.**
* `Q4_boxplots_shortlist.png` — compact horizontal overview of the confident shortlist (⚠ = likely contaminant/synthetic).
* `Q4_boxplots_top10.png` — literal raw top-10 per direction by Student's t (no filtering), matching the original wording.

> **Annotation caveats:** (a) some hits are likely **library contaminants / false matches** (e.g.,
> Lauryldiethanolamine, Polidocanol — surfactants; sulfa drugs; amitriptyline) — `MQScore`/`MZErrorPPM` are carried
> in every table. (b) `*`-flagged rows have **degenerate annotations** (one name, several distinct features) — treat
> the specific identity cautiously. (c) All matches are **putative** (MS²-library level) pending authentic-standard
> confirmation, and remain **descriptive** (technical replicates, n=1 biological/group).

---

## Output files

```
analysis/metabolite_analysis.py            Phase-1 primary pipeline
analysis/robustness_checks.py              PERMDISP, CLR/Aitchison, moderated-t, LOO, QC, shortlist
analysis/annotation_merge.py               Phase-2 annotation merge + Student-t
analysis/confident_shortlist.py            Phase-2b MQScore filter + adduct dedup + boxplots
analysis/add_resolved_names.py             resolved_name + catalog-id / synthetic flags
analysis/pubchem_enrich.py                 verified PubChem CID/URL/title (InChIKey, formula-checked)
analysis/boxplots_top10_literal.py         raw top-10 boxplots (matches collaborator's wording)
analysis/make_boxplots.py                  regenerates all boxplots from enriched tables (resolved_name)
results/figures/Q4_boxplots_grid.png       publication-standard faceted grid (recommended for paper)
results/tables/pubchem_resolution.csv      audit of catalog-ID -> PubChem structure rescues
results/figures/Q4_boxplots_top10.png      raw top-10 per direction (literal request)
results/figures/Q4_boxplots_shortlist.png  confident shortlist boxplots (refined)
results/tables/annotated_features_all_stats.csv      174 annotated features + stats + chem class
results/tables/annotated_significant_features.csv    74 annotated & significant (+ direction)
results/tables/annotated_significant_glucose.csv     44 higher in glucose
results/tables/annotated_significant_starch.csv      30 higher in starch
results/tables/confident_shortlist.csv               45 confident metabolites (MQ>=0.85, deduped)
results/tables/confident_shortlist_glucose.csv       26 higher in glucose
results/tables/confident_shortlist_starch.csv        19 higher in starch
results/figures/Q1_PCA.png                 PCA scores + scree
results/figures/Q1_PCoA_BrayCurtis.png     PCoA + PERMANOVA
results/figures/Q1_dendrograms.png         Ward & Bray-Curtis dendrograms
results/figures/Q1_heatmap_top50.png       clustered heatmap (top-50 variance features)
results/figures/Q2_normality.png           histograms, QQ, skewness
results/figures/Q3_volcano.png             volcano, all 7032 features
results/figures/Q3b_robustness.png         CLR concordance, PERMDISP, moderated-t
results/tables/differential_features.csv          Welch + Mann-Whitney, sorted by p
results/tables/differential_features_annotated.csv  + moderated-t, BY-FDR, CV, floor flag  <-- use for lookup
results/tables/shortlist_robust_features.csv      538 actionable features
results/tables/top25_features.csv                 preview
results/tables/qc_per_sample.csv                  per-sample QC
results/tables/distance_braycurtis.csv            8x8 distance matrix
results/tables/summary.json / summary_robustness.json   all numeric results
```

## Suggested next step
Resolve the `feature_id`s of the **robust shortlist** (`shortlist_robust_features.csv`) against your spectral/
molecular database (after the de-duplication caveat above) to identify the metabolites driving the
glucose-vs-starch difference — **once the replicate structure (Caveat 1) is confirmed.**
