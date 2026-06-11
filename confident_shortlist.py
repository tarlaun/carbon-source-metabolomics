#!/usr/bin/env python3
"""
Phase 2b — Confident biological shortlist + redesigned boxplots.

From the 74 annotated & significant features (Student's t, FDR<0.05 & |log2FC|>=1):
  1. Keep only confident library matches: MQScore >= MQ_MIN (default 0.85).
  2. Collapse adduct/isotope/charge-state redundancy: features that share the same
     InChIKey skeleton AND co-elute (|dRT| <= RT_TOL) are one metabolite -> keep the
     best-MQScore representative (record the collapsed ids).
     Same-skeleton features at DIFFERENT retention times are isomers / degenerate
     library matches (common for diketopiperazines) -> kept separate and FLAGGED.
  3. Flag any molecule whose features disagree on direction (annotation ambiguity).

Redesigned figure: a single plot with two stacked panels (higher-in-glucose,
higher-in-starch), horizontal grouped boxplots on a shared log abundance axis with
the 4 replicate points overlaid.
"""
import os, re, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
CSV = os.path.join(ROOT, "2026-06-06T21-12_export.csv")
ANN = os.path.join(ROOT, "e1fe0f44e01f4fd0a8deef0718c20c3d-merged_results_with_gnps.tsv")
FIG = os.path.join(ROOT, "results", "figures"); TAB = os.path.join(ROOT, "results", "tables")
GC = {"Starch": "#D55E00", "Glucose": "#0072B2"}
MQ_MIN = 0.85
RT_TOL = 0.20
TOP_N = 12


def clean_name(s):
    s = str(s).strip().strip('"').strip()
    for pat in [r"^Spectral Match to (.+?) from ", r"^Massbank:\S+\s+(.+)$",
                r"^ReSpect:\S+\s+(.+)$"]:
        m = re.match(pat, s)
        if m: s = m.group(1); break
    s = re.split(r"\s+CollisionEnergy", s)[0]
    s = re.split(r"\s+-\s+[\d.]+\s*eV", s)[0]
    s = s.split("|")[0].strip().strip('"').strip()
    return s if s else "(unnamed)"


def load_R():
    df = pd.read_csv(CSV).rename(columns={pd.read_csv(CSV, nrows=0).columns[0]: "sample"})
    meta = ["sample", "ATTRIBUTE_Sample", "ATTRIBUTE_Blank", "ATTRIBUTE_Blank_Ctrl",
            "ATTRIBUTE_CarbonSource"]
    feat = [c for c in df.columns if c not in meta]
    grp = df["ATTRIBUTE_CarbonSource"].replace({"Startch": "Starch"}).to_numpy()
    R = df[feat].apply(pd.to_numeric, errors="coerce"); R.index = df["sample"].to_numpy()
    fid2col = {int(c.split("_")[0]): c for c in feat}
    return R, grp, fid2col


def build_shortlist():
    sig = pd.read_csv(f"{TAB}/annotated_significant_features.csv")
    ann = pd.read_csv(ANN, sep="\t")
    sig = sig.merge(ann[["#Scan#", "InChIKey-Planar"]].rename(
        columns={"#Scan#": "feature_id"}), on="feature_id", how="left")

    def skel(r):
        for c in ["InChIKey-Planar", "InChIKey"]:
            v = r.get(c)
            if isinstance(v, str) and len(v) > 3:
                return v[:14]
        return "NAME:" + clean_name(r.Compound_Name).lower()
    sig["skel"] = sig.apply(skel, axis=1)
    sig["clean_name"] = sig.Compound_Name.map(clean_name)

    # confidence filter
    conf = sig[sig.MQScore >= MQ_MIN].copy()

    # cluster by skeleton + co-elution (RT gap <= RT_TOL)
    conf = conf.sort_values(["skel", "rt"])
    cluster_id = []
    last_skel, last_rt, cid = None, None, 0
    for _, r in conf.iterrows():
        if r.skel != last_skel or (r.rt - last_rt) > RT_TOL:
            cid += 1
        cluster_id.append(cid); last_skel, last_rt = r.skel, r.rt
    conf["cluster"] = cluster_id

    # per-skeleton annotation-ambiguity flags (computed on confident set)
    skel_dirs = conf.groupby("skel")["direction"].nunique()
    skel_nclu = conf.groupby("skel")["cluster"].nunique()

    reps = []
    for cl, sub in conf.groupby("cluster"):
        rep = sub.sort_values("MQScore", ascending=False).iloc[0].copy()
        rep["n_collapsed_adducts"] = len(sub)
        rep["collapsed_feature_ids"] = ";".join(map(str, sorted(sub.feature_id.astype(int))))
        rep["isomers_same_annotation"] = bool(skel_nclu[rep.skel] > 1)
        rep["direction_conflict_in_annotation"] = bool(skel_dirs[rep.skel] > 1)
        reps.append(rep)
    short = pd.DataFrame(reps)

    cols = ["feature_id", "clean_name", "Compound_Name", "mz", "rt", "direction",
            "log2FC_Glucose_vs_Starch", "mean_Starch", "mean_Glucose",
            "student_t", "student_p", "student_FDR", "MQScore", "MZErrorPPM",
            "n_collapsed_adducts", "collapsed_feature_ids", "isomers_same_annotation",
            "direction_conflict_in_annotation", "molecular_formula", "superclass",
            "class", "npclassifier_pathway", "Smiles", "InChIKey"]
    short = short[[c for c in cols if c in short.columns]].sort_values("student_p").reset_index(drop=True)
    return sig, short


def boxplots(short, R, grp, fid2col, path):
    glu = short[short.direction == "higher_in_Glucose"].head(TOP_N).iloc[::-1]
    sta = short[short.direction == "higher_in_Starch"].head(TOP_N).iloc[::-1]
    m = {g: grp == g for g in ("Starch", "Glucose")}

    fig, axes = plt.subplots(2, 1, figsize=(12, 13),
                             gridspec_kw={"height_ratios": [len(glu), len(sta)]})
    for ax, tab, title in [(axes[0], glu, f"Higher in GLUCOSE  (top {len(glu)})"),
                           (axes[1], sta, f"Higher in STARCH  (top {len(sta)})")]:
        ax.set_xscale("log")
        for i, (_, r) in enumerate(tab.iterrows()):
            col = fid2col[int(r.feature_id)]
            for g, off in [("Starch", -0.20), ("Glucose", +0.20)]:
                v = R[col].values[m[g]]
                bp = ax.boxplot([v], positions=[i + off], widths=0.34, vert=False,
                                patch_artist=True, showfliers=False, manage_ticks=False)
                bp["boxes"][0].set_facecolor(GC[g]); bp["boxes"][0].set_alpha(.35)
                bp["medians"][0].set_color("k")
                ax.scatter(v, np.full(len(v), i + off) + np.random.normal(0, .03, len(v)),
                           color=GC[g], edgecolor="k", s=28, zorder=3)
        labels = [f"{r.clean_name[:30]}  (id {int(r.feature_id)})"
                  + ("  *" if r.isomers_same_annotation else "")
                  for _, r in tab.iterrows()]
        ax.set_yticks(range(len(tab))); ax.set_yticklabels(labels, fontsize=9)
        # annotate log2FC at right edge
        xr = ax.get_xlim()[1]
        for i, (_, r) in enumerate(tab.iterrows()):
            ax.text(xr, i, f" log2FC {r.log2FC_Glucose_vs_Starch:+.1f}",
                    va="center", fontsize=7.5, color="#333")
        ax.set_xlabel("relative abundance (log scale)")
        ax.set_title(title, fontsize=12, loc="left")
        ax.grid(axis="x", alpha=.3)
        ax.set_ylim(-0.6, len(tab) - 0.4)
    leg = [Patch(facecolor=GC["Starch"], alpha=.5, label="Starch"),
           Patch(facecolor=GC["Glucose"], alpha=.5, label="Glucose")]
    axes[0].legend(handles=leg, loc="lower right", title="Carbon source")
    fig.suptitle("Confident annotated metabolites differing by carbon source\n"
                 f"(GNPS MQScore ≥ {MQ_MIN}, adducts collapsed; * = isomeric/degenerate annotation)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def main():
    np.random.seed(42)
    R, grp, fid2col = load_R()
    sig, short = build_shortlist()

    short.to_csv(f"{TAB}/confident_shortlist.csv", index=False)
    short[short.direction == "higher_in_Glucose"].to_csv(f"{TAB}/confident_shortlist_glucose.csv", index=False)
    short[short.direction == "higher_in_Starch"].to_csv(f"{TAB}/confident_shortlist_starch.csv", index=False)

    boxplots(short, R, grp, fid2col, f"{FIG}/Q4_boxplots_shortlist.png")
    # remove the old grid figures
    for old in ["Q4_boxplots_glucose_top10.png", "Q4_boxplots_starch_top10.png"]:
        p = os.path.join(FIG, old)
        if os.path.exists(p): os.remove(p)

    summ = {
        "MQ_MIN": MQ_MIN, "RT_TOL_min": RT_TOL,
        "annotated_significant_in": int(len(sig)),
        "after_MQ_filter": int((sig.MQScore >= MQ_MIN).sum()),
        "confident_metabolites_after_dedup": int(len(short)),
        "higher_in_glucose": int((short.direction == "higher_in_Glucose").sum()),
        "higher_in_starch": int((short.direction == "higher_in_Starch").sum()),
        "flagged_isomeric_degenerate": int(short.isomers_same_annotation.sum()),
        "flagged_direction_conflict": int(short.direction_conflict_in_annotation.sum()),
    }
    with open(f"{TAB}/confident_shortlist_summary.json", "w") as f:
        json.dump(summ, f, indent=2)
    print(json.dumps(summ, indent=2))
    print("\nTop confident GLUCOSE-enriched:")
    print(short[short.direction == "higher_in_Glucose"].head(12)[
        ["feature_id", "clean_name", "log2FC_Glucose_vs_Starch", "student_FDR",
         "MQScore", "n_collapsed_adducts", "isomers_same_annotation"]].to_string(index=False))
    print("\nTop confident STARCH-enriched:")
    print(short[short.direction == "higher_in_Starch"].head(12)[
        ["feature_id", "clean_name", "log2FC_Glucose_vs_Starch", "student_FDR",
         "MQScore", "n_collapsed_adducts", "isomers_same_annotation"]].to_string(index=False))


if __name__ == "__main__":
    main()
