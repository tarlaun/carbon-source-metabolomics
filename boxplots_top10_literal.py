#!/usr/bin/env python3
"""
Literal-request figure: boxplot of the TEN annotated features with the greatest
statistically significant differences (Student's t-test), for glucose-prevalent
and starch-prevalent features -- WITHOUT the MQScore filter or adduct collapsing
(i.e. raw top-10 per direction, exactly as the collaborator phrased it).

Same single-figure style as Q4_boxplots_shortlist.png.
"""
import os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
CSV = os.path.join(ROOT, "2026-06-06T21-12_export.csv")
FIG = os.path.join(ROOT, "results", "figures"); TAB = os.path.join(ROOT, "results", "tables")
GC = {"Starch": "#D55E00", "Glucose": "#0072B2"}
TOP_N = 10


def clean_name(s):
    s = str(s).strip().strip('"').strip()
    for pat in [r"^Spectral Match to (.+?) from ", r"^Massbank:\S+\s+(.+)$", r"^ReSpect:\S+\s+(.+)$"]:
        m = re.match(pat, s)
        if m: s = m.group(1); break
    s = re.split(r"\s+CollisionEnergy", s)[0]
    s = re.split(r"\s+-\s+[\d.]+\s*eV", s)[0]
    s = s.split("|")[0].strip().strip('"').strip()
    return s if s else "(unnamed)"


def load_R():
    df = pd.read_csv(CSV).rename(columns={pd.read_csv(CSV, nrows=0).columns[0]: "sample"})
    meta = ["sample", "ATTRIBUTE_Sample", "ATTRIBUTE_Blank", "ATTRIBUTE_Blank_Ctrl", "ATTRIBUTE_CarbonSource"]
    feat = [c for c in df.columns if c not in meta]
    grp = df["ATTRIBUTE_CarbonSource"].replace({"Startch": "Starch"}).to_numpy()
    R = df[feat].apply(pd.to_numeric, errors="coerce"); R.index = df["sample"].to_numpy()
    fid2col = {int(c.split("_")[0]): c for c in feat}
    return R, grp, fid2col


def main():
    np.random.seed(42)
    R, grp, fid2col = load_R()
    m = {g: grp == g for g in ("Starch", "Glucose")}
    sig = pd.read_csv(f"{TAB}/annotated_significant_features.csv")
    sig["clean_name"] = sig.Compound_Name.map(clean_name)

    glu = sig[sig.direction == "higher_in_Glucose"].sort_values("student_p").head(TOP_N).iloc[::-1]
    sta = sig[sig.direction == "higher_in_Starch"].sort_values("student_p").head(TOP_N).iloc[::-1]

    fig, axes = plt.subplots(2, 1, figsize=(12, 12),
                             gridspec_kw={"height_ratios": [len(glu), len(sta)]})
    for ax, tab, title in [(axes[0], glu, f"Top {len(glu)} annotated features — higher in GLUCOSE"),
                           (axes[1], sta, f"Top {len(sta)} annotated features — higher in STARCH")]:
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
        ax.set_yticks(range(len(tab)))
        ax.set_yticklabels([f"{r.clean_name[:30]}  (id {int(r.feature_id)})" for _, r in tab.iterrows()],
                           fontsize=9)
        xr = ax.get_xlim()[1]
        for i, (_, r) in enumerate(tab.iterrows()):
            ax.text(xr, i, f" log2FC {r.log2FC_Glucose_vs_Starch:+.1f}  (p={r.student_p:.1e})",
                    va="center", fontsize=7.5, color="#333")
        ax.set_xlabel("relative abundance (log scale)")
        ax.set_title(title, fontsize=12, loc="left")
        ax.grid(axis="x", alpha=.3); ax.set_ylim(-0.6, len(tab) - 0.4)
    axes[0].legend(handles=[Patch(facecolor=GC["Starch"], alpha=.5, label="Starch"),
                            Patch(facecolor=GC["Glucose"], alpha=.5, label="Glucose")],
                   loc="lower right", title="Carbon source")
    fig.suptitle("Top-10 annotated significant features by Student's t-test\n"
                 "(all annotated hits ranked by p; no MQScore filter / no adduct de-duplication)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = f"{FIG}/Q4_boxplots_top10.png"
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)
    print("\nGLUCOSE top-10 (id | name | p):")
    for _, r in glu.iloc[::-1].iterrows():
        print(f"  {int(r.feature_id):>6}  {r.clean_name[:40]:<40}  {r.student_p:.2e}")
    print("STARCH top-10 (id | name | p):")
    for _, r in sta.iloc[::-1].iterrows():
        print(f"  {int(r.feature_id):>6}  {r.clean_name[:40]:<40}  {r.student_p:.2e}")


if __name__ == "__main__":
    main()
