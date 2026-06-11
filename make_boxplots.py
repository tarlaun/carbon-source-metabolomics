#!/usr/bin/env python3
"""
(Re)generate the Phase-2 boxplots from the ALREADY-ENRICHED tables (uses
resolved_name + PubChem-resolved labels; does NOT rebuild/overwrite the tables).

Figures:
  Q4_boxplots_grid.png      CONVENTIONAL metabolomics format: faceted grid, one panel
                            per metabolite, groups on x, abundance on log-y, individual
                            replicate points overlaid, numeric BH-FDR per panel (no stars,
                            since at n=4 technical reps the p-values are descriptive only).
                            Top-10 per direction from the confident, non-flagged set.
  Q4_boxplots_shortlist.png compact horizontal overview (regenerated with resolved_name)
  Q4_boxplots_top10.png     literal raw top-10 per direction (regenerated with resolved_name)
"""
import os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "2026-06-06T21-12_export.csv")
FIG = os.path.join(ROOT, "results", "figures"); TAB = os.path.join(ROOT, "results", "tables")
GC = {"Starch": "#D55E00", "Glucose": "#0072B2"}


def load_R():
    df = pd.read_csv(CSV).rename(columns={pd.read_csv(CSV, nrows=0).columns[0]: "sample"})
    meta = ["sample", "ATTRIBUTE_Sample", "ATTRIBUTE_Blank", "ATTRIBUTE_Blank_Ctrl", "ATTRIBUTE_CarbonSource"]
    feat = [c for c in df.columns if c not in meta]
    grp = df["ATTRIBUTE_CarbonSource"].replace({"Startch": "Starch"}).to_numpy()
    R = df[feat].apply(pd.to_numeric, errors="coerce"); R.index = df["sample"].to_numpy()
    fid2col = {int(c.split("_")[0]): c for c in feat}
    return R, grp, fid2col


def lab(name, n=26):
    s = str(name).strip()
    return s if len(s) <= n else s[:n - 1] + "…"


R, grp, fid2col = load_R()
mS, mG = grp == "Starch", grp == "Glucose"
np.random.seed(42)


def panel(ax, r, leftcol):
    col = fid2col[int(r.feature_id)]
    sv, gv = R[col].values[mS], R[col].values[mG]
    bp = ax.boxplot([sv, gv], positions=[1, 2], widths=0.55, patch_artist=True, showfliers=False)
    for patch, g in zip(bp["boxes"], ["Starch", "Glucose"]):
        patch.set_facecolor(GC[g]); patch.set_alpha(.35); patch.set_edgecolor(GC[g])
    for med in bp["medians"]:
        med.set_color("k")
    for pos, vals, g in [(1, sv, "Starch"), (2, gv, "Glucose")]:
        ax.scatter(np.random.normal(pos, 0.06, len(vals)), vals, color=GC[g],
                   edgecolor="k", s=34, lw=.6, zorder=3)
    ax.set_yscale("log")
    ax.set_xlim(0.4, 2.6); ax.set_xticks([1, 2]); ax.set_xticklabels(["Starch", "Glucose"], fontsize=8.5)
    ax.tick_params(axis="y", labelsize=7)
    top = max(sv.max(), gv.max())
    ax.set_ylim(top=top * 4)   # headroom above the boxes
    ax.set_title(f"{lab(r.resolved_name)}\nid {int(r.feature_id)} · log2FC "
                 f"{r.log2FC_Glucose_vs_Starch:+.1f} · FDR {r.student_FDR:.0e}", fontsize=8.0)
    if leftcol:
        ax.set_ylabel("relative abundance", fontsize=8)


def grid_figure(short, path):
    clean = short[~short.likely_synthetic_or_contaminant]
    glu = clean[clean.direction == "higher_in_Glucose"].sort_values("student_p").head(10)
    sta = clean[clean.direction == "higher_in_Starch"].sort_values("student_p").head(10)
    fig = plt.figure(figsize=(19, 15))
    subfigs = fig.subfigures(2, 1, hspace=0.06)
    for sf, tab, title in [(subfigs[0], glu, "Higher in GLUCOSE"),
                           (subfigs[1], sta, "Higher in STARCH")]:
        sf.suptitle(f"{title} — top {len(tab)} confident annotated metabolites "
                    f"(Student's t-test; BH-FDR per panel)", fontsize=13, fontweight="bold")
        axes = sf.subplots(2, 5)
        for i, ax in enumerate(axes.flat):
            if i < len(tab):
                panel(ax, tab.iloc[i], leftcol=(i % 5 == 0))
            else:
                ax.axis("off")
    fig.suptitle("Carbon-source–differential annotated metabolites  (n=4 technical replicates/group; "
                 "FDR = BH-adjusted Student's t — descriptive, not biological inference)", fontsize=14, y=1.005)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def horizontal(tab_by_dir, path, suptitle, topn, exclude_flagged):
    fig, axes = plt.subplots(2, 1, figsize=(12, 12.5),
                             gridspec_kw={"height_ratios": [topn, topn]})
    for ax, (tab, title) in zip(axes, tab_by_dir):
        t = tab[~tab.likely_synthetic_or_contaminant] if exclude_flagged else tab
        t = t.sort_values("student_p").head(topn).iloc[::-1]
        ax.set_xscale("log")
        for i, (_, r) in enumerate(t.iterrows()):
            col = fid2col[int(r.feature_id)]
            for g, off in [("Starch", -0.20), ("Glucose", +0.20)]:
                v = R[col].values[mS if g == "Starch" else mG]
                bp = ax.boxplot([v], positions=[i + off], widths=0.34, vert=False,
                                patch_artist=True, showfliers=False, manage_ticks=False)
                bp["boxes"][0].set_facecolor(GC[g]); bp["boxes"][0].set_alpha(.35)
                bp["medians"][0].set_color("k")
                ax.scatter(v, np.full(len(v), i + off) + np.random.normal(0, .03, len(v)),
                           color=GC[g], edgecolor="k", s=26, zorder=3)
        flag = ["  ⚠" if r.likely_synthetic_or_contaminant else "" for _, r in t.iterrows()]
        ax.set_yticks(range(len(t)))
        ax.set_yticklabels([f"{lab(r.resolved_name,30)}  (id {int(r.feature_id)}){flag[i]}"
                            for i, (_, r) in enumerate(t.iterrows())], fontsize=9)
        xr = ax.get_xlim()[1]
        for i, (_, r) in enumerate(t.iterrows()):
            ax.text(xr, i, f" log2FC {r.log2FC_Glucose_vs_Starch:+.1f} · FDR {r.student_FDR:.0e}",
                    va="center", fontsize=7.5, color="#333")
        ax.set_xlabel("relative abundance (log scale)"); ax.set_title(title, fontsize=12, loc="left")
        ax.grid(axis="x", alpha=.3); ax.set_ylim(-0.6, len(t) - 0.4)
    axes[0].legend(handles=[Patch(facecolor=GC["Starch"], alpha=.5, label="Starch"),
                            Patch(facecolor=GC["Glucose"], alpha=.5, label="Glucose")],
                   loc="lower right", title="Carbon source")
    fig.suptitle(suptitle, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def main():
    short = pd.read_csv(f"{TAB}/confident_shortlist.csv")
    sig = pd.read_csv(f"{TAB}/annotated_significant_features.csv")

    grid_figure(short, f"{FIG}/Q4_boxplots_grid.png")

    horizontal([(short[short.direction == "higher_in_Glucose"], "Higher in GLUCOSE"),
                (short[short.direction == "higher_in_Starch"], "Higher in STARCH")],
               f"{FIG}/Q4_boxplots_shortlist.png",
               "Confident annotated metabolites by carbon source (MQScore≥0.85, adducts collapsed; ⚠ = likely contaminant/synthetic)",
               12, exclude_flagged=False)

    horizontal([(sig[sig.direction == "higher_in_Glucose"], "Higher in GLUCOSE"),
                (sig[sig.direction == "higher_in_Starch"], "Higher in STARCH")],
               f"{FIG}/Q4_boxplots_top10.png",
               "Top-10 annotated significant features by Student's t (literal request; no filtering; ⚠ = likely contaminant/synthetic)",
               10, exclude_flagged=False)
    print("wrote Q4_boxplots_grid.png, Q4_boxplots_shortlist.png, Q4_boxplots_top10.png")


if __name__ == "__main__":
    main()
