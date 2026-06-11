#!/usr/bin/env python3
"""
Phase 2 — Cross-reference GNPS annotations with the statistically significant
features (Starch vs Glucose), split by direction, and visualise the top hits.

Per the user's stated assumption (normal, equal variance) the differential test
here is **Student's t-test (equal_var=True)** on log10 relative abundances, with
Benjamini-Hochberg FDR. Welch results are carried alongside for transparency.

NOTE: the 4 samples per group are TECHNICAL replicates (confirmed by the user),
so all p-values are descriptive of these specific extracts, not a biological
population inference.

Merge key: dataset feature column 'ID_mz_rt&annot' -> leading ID == annotation '#Scan#'.

Outputs (results/tables, results/figures):
  annotated_features_all_stats.csv         all 174 annotated features + stats
  annotated_significant_features.csv       significant annotated features (+direction)
  annotated_significant_glucose.csv        ... higher in glucose
  annotated_significant_starch.csv         ... higher in starch
  Q4_boxplots_glucose_top10.png            top-10 glucose-prevalent annotated features
  Q4_boxplots_starch_top10.png             top-10 starch-prevalent annotated features
"""
import os, re, json
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
CSV = os.path.join(ROOT, "2026-06-06T21-12_export.csv")
ANN = os.path.join(ROOT, "e1fe0f44e01f4fd0a8deef0718c20c3d-merged_results_with_gnps.tsv")
FIG = os.path.join(ROOT, "results", "figures"); TAB = os.path.join(ROOT, "results", "tables")
GC = {"Starch": "#D55E00", "Glucose": "#0072B2"}
FC_THR, FDR_THR = 1.0, 0.05

ANN_COLS = ["Compound_Name", "Adduct", "molecular_formula", "ExactMass", "MQScore",
            "MZErrorPPM", "SharedPeaks", "Smiles", "InChIKey", "superclass", "class",
            "subclass", "npclassifier_superclass", "npclassifier_class", "npclassifier_pathway"]


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1); r = np.minimum.accumulate(r[::-1])[::-1]
    q = np.empty(n); q[o] = np.clip(r, 0, 1); return q


def clean_name(s):
    s = str(s).strip().strip('"').strip()
    m = re.match(r"Spectral Match to (.+?) from ", s)
    if m: s = m.group(1)
    s = re.split(r"\s+CollisionEnergy", s)[0]
    s = re.split(r"\s+-\s+[\d.]+\s*eV", s)[0]
    return s.strip().strip('"').strip()


def load_dataset():
    df = pd.read_csv(CSV).rename(columns={pd.read_csv(CSV, nrows=0).columns[0]: "sample"})
    meta = ["sample", "ATTRIBUTE_Sample", "ATTRIBUTE_Blank", "ATTRIBUTE_Blank_Ctrl",
            "ATTRIBUTE_CarbonSource"]
    feat = [c for c in df.columns if c not in meta]
    grp = df["ATTRIBUTE_CarbonSource"].replace({"Startch": "Starch"}).to_numpy()
    R = df[feat].apply(pd.to_numeric, errors="coerce"); R.index = df["sample"].to_numpy()
    info = []
    for c in feat:
        m = re.match(r"^(\d+)_([\d.]+)_([\d.]+)", c)
        info.append((c, int(m.group(1)), float(m.group(2)), float(m.group(3))))
    info = pd.DataFrame(info, columns=["col", "feature_id", "mz", "rt"])
    return R, grp, info


def differential(R, grp, info):
    L = np.log10(R); m0, m1 = grp == "Starch", grp == "Glucose"
    log2fc = (L.values[m1].mean(0) - L.values[m0].mean(0)) / np.log10(2)
    # Student's t (equal variance) -- the user's chosen test
    t_s, p_s = stats.ttest_ind(L.values[m1], L.values[m0], axis=0, equal_var=True)
    # Welch (reference)
    t_w, p_w = stats.ttest_ind(L.values[m1], L.values[m0], axis=0, equal_var=False)
    res = info.copy()
    res["mean_Starch"] = R.values[m0].mean(0)
    res["mean_Glucose"] = R.values[m1].mean(0)
    res["log2FC_Glucose_vs_Starch"] = log2fc
    res["student_t"] = t_s; res["student_p"] = p_s; res["student_FDR"] = bh_fdr(p_s)
    res["welch_p"] = p_w; res["welch_FDR"] = bh_fdr(p_w)
    res["significant"] = (res.student_FDR < FDR_THR) & (res.log2FC_Glucose_vs_Starch.abs() >= FC_THR)
    res["direction"] = np.where(log2fc > 0, "higher_in_Glucose", "higher_in_Starch")
    return res, L


def boxplot_panel(features, R, grp, title, path):
    n = len(features); ncol = 5; nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow), squeeze=False)
    m0, m1 = grp == "Starch", grp == "Glucose"
    for i, (_, r) in enumerate(features.iterrows()):
        ax = axes[i // ncol][i % ncol]
        col = r["col"]
        vals = {"Starch": R[col].values[m0], "Glucose": R[col].values[m1]}
        bp = ax.boxplot([vals["Starch"], vals["Glucose"]], labels=["Starch", "Glucose"],
                        widths=.55, patch_artist=True, showfliers=False, zorder=1)
        for patch, g in zip(bp["boxes"], ["Starch", "Glucose"]):
            patch.set_facecolor(GC[g]); patch.set_alpha(.35)
        for med in bp["medians"]: med.set_color("k")
        for j, g in enumerate(["Starch", "Glucose"], start=1):
            x = np.random.normal(j, 0.06, len(vals[g]))
            ax.scatter(x, vals[g], color=GC[g], edgecolor="k", s=45, zorder=3)
        nm = clean_name(r["Compound_Name"])
        ax.set_title(f"{nm[:34]}\n(id {int(r.feature_id)}, m/z {r.mz:.3f})", fontsize=9)
        ax.set_ylabel("relative abundance", fontsize=8)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.text(0.5, 0.97, f"log2FC={r.log2FC_Glucose_vs_Starch:+.1f}  "
                f"p={r.student_p:.1e}\nFDR={r.student_FDR:.1e}  MQ={r.MQScore:.2f}",
                transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
                bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=.8))
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle(title, fontsize=14, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def main():
    np.random.seed(42)
    R, grp, info = load_dataset()
    res, L = differential(R, grp, info)
    ann = pd.read_csv(ANN, sep="\t")
    ann_sub = ann[["#Scan#"] + ANN_COLS].rename(columns={"#Scan#": "feature_id"})

    # merge annotations onto the per-feature stats (inner = annotated features only)
    merged = res.merge(ann_sub, on="feature_id", how="inner")
    print(f"annotated features present in dataset: {len(merged)}")

    order_cols = (["feature_id", "mz", "rt", "Compound_Name", "Adduct", "molecular_formula",
                   "ExactMass", "direction", "mean_Starch", "mean_Glucose",
                   "log2FC_Glucose_vs_Starch", "student_t", "student_p", "student_FDR",
                   "welch_p", "welch_FDR", "significant", "MQScore", "MZErrorPPM",
                   "SharedPeaks", "superclass", "class", "subclass", "npclassifier_pathway",
                   "npclassifier_superclass", "npclassifier_class", "Smiles", "InChIKey"])
    allt = merged[order_cols].sort_values("student_p").reset_index(drop=True)
    allt.to_csv(f"{TAB}/annotated_features_all_stats.csv", index=False)

    sig = allt[allt.significant].copy()
    sig.to_csv(f"{TAB}/annotated_significant_features.csv", index=False)
    glu = sig[sig.direction == "higher_in_Glucose"].sort_values("student_p")
    sta = sig[sig.direction == "higher_in_Starch"].sort_values("student_p")
    glu.to_csv(f"{TAB}/annotated_significant_glucose.csv", index=False)
    sta.to_csv(f"{TAB}/annotated_significant_starch.csv", index=False)

    # agreement Student vs Welch significance
    welch_sig = set(res.loc[(res.welch_FDR < FDR_THR) & (res.log2FC_Glucose_vs_Starch.abs() >= FC_THR), "feature_id"])
    stud_sig = set(res.loc[res.significant, "feature_id"])
    jac = len(welch_sig & stud_sig) / len(welch_sig | stud_sig)

    summary = {
        "annotated_in_dataset": int(len(merged)),
        "annotated_AND_significant": int(len(sig)),
        "annotated_sig_higher_in_glucose": int(len(glu)),
        "annotated_sig_higher_in_starch": int(len(sta)),
        "student_vs_welch_sigset_jaccard": round(jac, 3),
        "test": "Student t-test (equal_var=True) on log10 relative abundance; BH-FDR; "
                "significant = FDR<0.05 & |log2FC|>=1",
    }
    print(json.dumps(summary, indent=2))

    # need 'col' for plotting -> re-attach from res
    colmap = res.set_index("feature_id")["col"]
    for tab in (glu, sta):
        tab["col"] = tab["feature_id"].map(colmap)

    # boxplots: top 10 per direction by Student p
    n_g = min(10, len(glu)); n_s = min(10, len(sta))
    boxplot_panel(glu.head(n_g), R, grp,
                  f"Top {n_g} annotated features higher in GLUCOSE (Student's t-test)",
                  f"{FIG}/Q4_boxplots_glucose_top10.png")
    boxplot_panel(sta.head(n_s), R, grp,
                  f"Top {n_s} annotated features higher in STARCH (Student's t-test)",
                  f"{FIG}/Q4_boxplots_starch_top10.png")
    print(f"boxplots written: {n_g} glucose, {n_s} starch")

    with open(f"{TAB}/annotation_merge_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== TOP 10 higher in GLUCOSE ===")
    print(glu.head(10)[["feature_id", "Compound_Name", "log2FC_Glucose_vs_Starch",
                        "student_p", "student_FDR", "MQScore"]].to_string(index=False))
    print("\n=== TOP 10 higher in STARCH ===")
    print(sta.head(10)[["feature_id", "Compound_Name", "log2FC_Glucose_vs_Starch",
                        "student_p", "student_FDR", "MQScore"]].to_string(index=False))


if __name__ == "__main__":
    main()
