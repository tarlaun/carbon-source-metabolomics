#!/usr/bin/env python3
"""
Phase 1 — Unsupervised + univariate analysis of bacterial metabolite profiles
under two carbon sources (Starch vs Glucose).

Input : 2026-06-06T21-12_export.csv  (already blank-removed, imputed, TIC-normalised)
Design: 8 samples (4 Starch 'EXT_*' + 4 Glucose 'aG_*'), 7032 features.
        Rows sum to 1.0 -> relative abundances (compositional data).

Answers three questions:
  Q1  Is the metabolite profile different between groups?
        -> PCA, PCoA (Bray-Curtis), PERMANOVA (exact), hierarchical clustering
  Q2  Do the data follow a normal distribution?
        -> aggregate + per-feature normality assessment (with small-n caveats)
  Q3  Which features differ between groups?
        -> Welch t-test (log10) + BH-FDR  (parametric, primary)
           Mann-Whitney U + BH-FDR        (non-parametric, sensitivity)
           Volcano plot of the full dataset

All preprocessing is done ONCE here so every method sees a consistent matrix.
"""

import os
import re
import json
import itertools
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform, pdist
from sklearn.decomposition import PCA

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore", category=RuntimeWarning)
sns.set_style("whitegrid")
RNG = np.random.default_rng(42)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CSV = os.path.join(ROOT, "2026-06-06T21-12_export.csv")
FIG = os.path.join(ROOT, "results", "figures")
TAB = os.path.join(ROOT, "results", "tables")
os.makedirs(FIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

GROUP_COLORS = {"Starch": "#D55E00", "Glucose": "#0072B2"}
summary = {}   # collected numeric results -> JSON at the end


# ----------------------------------------------------------------------------
# 0. Load & validate
# ----------------------------------------------------------------------------
def load():
    df = pd.read_csv(CSV)
    df = df.rename(columns={df.columns[0]: "sample"})
    meta_cols = ["sample", "ATTRIBUTE_Sample", "ATTRIBUTE_Blank",
                 "ATTRIBUTE_Blank_Ctrl", "ATTRIBUTE_CarbonSource"]
    feat_cols = [c for c in df.columns if c not in meta_cols]

    # 'Startch' in the file -> normalise label to 'Starch'
    grp = df["ATTRIBUTE_CarbonSource"].replace({"Startch": "Starch"}).to_numpy()
    X = df[feat_cols].apply(pd.to_numeric, errors="coerce")
    X.index = df["sample"].to_numpy()

    assert X.isna().sum().sum() == 0, "unexpected NaN"
    assert (X.values >= 0).all(), "unexpected negative"
    rowsums = X.sum(axis=1)
    print(f"[load] {X.shape[0]} samples x {X.shape[1]} features | "
          f"groups={dict(pd.Series(grp).value_counts())} | "
          f"rowsum range [{rowsums.min():.4f}, {rowsums.max():.4f}]")

    # parse feature id / m-z / rt
    feat_info = []
    for c in feat_cols:
        m = re.match(r"^(\d+)_([\d.]+)_([\d.]+)", c)
        feat_info.append((c, m.group(1), float(m.group(2)), float(m.group(3))))
    feat_info = pd.DataFrame(feat_info, columns=["col", "feature_id", "mz", "rt"])

    summary["n_samples"] = int(X.shape[0])
    summary["n_features"] = int(X.shape[1])
    summary["groups"] = {k: int(v) for k, v in pd.Series(grp).value_counts().items()}
    return X, grp, feat_info


# ----------------------------------------------------------------------------
# helpers: Bray-Curtis, PCoA, PERMANOVA (exact), BH-FDR
# ----------------------------------------------------------------------------
def bray_curtis(M):
    return squareform(pdist(M, metric="braycurtis"))


def pcoa(D):
    """Classical MDS (Gower) on a distance matrix."""
    n = D.shape[0]
    A = -0.5 * (D ** 2)
    J = np.eye(n) - np.ones((n, n)) / n
    B = J @ A @ J
    vals, vecs = np.linalg.eigh(B)
    idx = np.argsort(vals)[::-1]
    vals, vecs = vals[idx], vecs[:, idx]
    pos = vals > 1e-9
    coords = vecs[:, pos] * np.sqrt(vals[pos])
    prop = vals[pos] / vals[pos].sum()
    return coords, prop, vals


def permanova_exact(D, grp):
    """Anderson (2001) PERMANOVA; exact test over all label permutations."""
    n = D.shape[0]
    labels = np.unique(grp)
    a = len(labels)
    iu = np.triu_indices(n, 1)
    d2 = D[iu] ** 2

    def ss(group_vec):
        sst = d2.sum() / n
        ssw = 0.0
        for lab in np.unique(group_vec):
            members = np.where(group_vec == lab)[0]
            if len(members) < 2:
                continue
            sub = D[np.ix_(members, members)]
            iu_s = np.triu_indices(len(members), 1)
            ssw += (sub[iu_s] ** 2).sum() / len(members)
        ssa = sst - ssw
        return ssa, ssw, sst

    ssa, ssw, sst = ss(grp)
    F_obs = (ssa / (a - 1)) / (ssw / (n - a))
    R2 = ssa / sst

    # exact: enumerate all C(n, n1) group assignments (2 groups)
    n1 = int((grp == labels[0]).sum())
    perms = list(itertools.combinations(range(n), n1))
    F_perm = []
    for combo in perms:
        gp = np.array([labels[1]] * n)
        gp[list(combo)] = labels[0]
        ssa_p, ssw_p, _ = ss(gp)
        F_perm.append((ssa_p / (a - 1)) / (ssw_p / (n - a)))
    F_perm = np.array(F_perm)
    # p = fraction of permuted F >= observed (observed is one of the perms)
    p = (F_perm >= F_obs - 1e-12).sum() / len(F_perm)
    # smallest attainable p = (# permutations tied at the global max F)/n_perm.
    # For a balanced 2-group design a partition and its complement give the SAME
    # pseudo-F, so the floor is 2/n_perm (= 2/70 = 0.0286), not 1/70.
    max_F = F_perm.max()
    n_at_max = int((F_perm >= max_F - 1e-9).sum())
    return dict(F=float(F_obs), R2=float(R2), p=float(p),
                n_perm=len(perms), min_attainable_p=n_at_max / len(perms),
                observed_is_max=bool(F_obs >= max_F - 1e-9))


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


# ----------------------------------------------------------------------------
# 1. preprocessing matrices
# ----------------------------------------------------------------------------
def preprocess(X):
    R = X.copy()                       # relative abundance (rows sum to 1)
    L = np.log10(R)                    # log10 (no zeros -> safe)
    # Pareto scaling of L: mean-centre, divide by sqrt(SD)  (metabolomics std)
    mu = L.mean(axis=0)
    sd = L.std(axis=0, ddof=1).replace(0, np.nan)
    P = (L - mu) / np.sqrt(sd)
    P = P.fillna(0.0)
    # autoscale (UV) for sensitivity
    A = (L - mu) / sd
    A = A.fillna(0.0)
    return R, L, P, A


# ----------------------------------------------------------------------------
# Q1 — group separation
# ----------------------------------------------------------------------------
def q1_pca(P, grp):
    pca = PCA(n_components=min(P.shape[0] - 1, 7))
    sc = pca.fit_transform(P.values)
    ev = pca.explained_variance_ratio_ * 100
    summary["pca_explained_var_pct"] = [round(float(x), 2) for x in ev]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
    for lab in np.unique(grp):
        m = grp == lab
        ax1.scatter(sc[m, 0], sc[m, 1], s=120, label=lab,
                    color=GROUP_COLORS[lab], edgecolor="k", zorder=3)
    for i, name in enumerate(P.index):
        ax1.annotate(name, (sc[i, 0], sc[i, 1]), fontsize=7,
                     xytext=(4, 4), textcoords="offset points")
    ax1.set_xlabel(f"PC1 ({ev[0]:.1f}%)")
    ax1.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
    ax1.set_title("PCA scores (log10 + Pareto scaling)")
    ax1.legend(title="Carbon source")
    ax1.axhline(0, color="grey", lw=.6); ax1.axvline(0, color="grey", lw=.6)

    ax2.bar(np.arange(1, len(ev) + 1), ev, color="#444")
    ax2.plot(np.arange(1, len(ev) + 1), np.cumsum(ev), "-o", color="#D55E00")
    ax2.set_xlabel("Principal component"); ax2.set_ylabel("Variance explained (%)")
    ax2.set_title("Scree (bars) + cumulative (line)")
    fig.tight_layout(); fig.savefig(f"{FIG}/Q1_PCA.png", dpi=150); plt.close(fig)


def q1_pcoa_permanova(R, L, grp):
    D_bc = bray_curtis(R.values)
    D_eu = squareform(pdist(L.values, metric="euclidean"))
    np.savetxt(f"{TAB}/distance_braycurtis.csv", D_bc, delimiter=",")

    res = {}
    for name, D in [("BrayCurtis(rel.abund.)", D_bc), ("Euclidean(log10)", D_eu)]:
        pm = permanova_exact(D, grp)
        res[name] = pm
        print(f"[PERMANOVA | {name}] F={pm['F']:.3f}  R2={pm['R2']:.3f}  "
              f"p={pm['p']:.4f}  (exact, {pm['n_perm']} perms; floor "
              f"min-p={pm['min_attainable_p']:.4f}; observed_is_max="
              f"{pm['observed_is_max']})")
    summary["permanova"] = res

    # PCoA on Bray-Curtis
    coords, prop, _ = pcoa(D_bc)
    summary["pcoa_explained_var_pct"] = [round(float(x * 100), 2) for x in prop[:5]]
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for lab in np.unique(grp):
        m = grp == lab
        ax.scatter(coords[m, 0], coords[m, 1], s=120, label=lab,
                   color=GROUP_COLORS[lab], edgecolor="k", zorder=3)
    for i, name in enumerate(R.index):
        ax.annotate(name, (coords[i, 0], coords[i, 1]), fontsize=7,
                    xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel(f"PCoA1 ({prop[0]*100:.1f}%)")
    ax.set_ylabel(f"PCoA2 ({prop[1]*100:.1f}%)")
    pm = res["BrayCurtis(rel.abund.)"]
    ax.set_title(f"PCoA (Bray-Curtis)\nPERMANOVA: F={pm['F']:.2f}, "
                 f"R²={pm['R2']:.2f}, p={pm['p']:.3f}")
    ax.legend(title="Carbon source")
    ax.axhline(0, color="grey", lw=.6); ax.axvline(0, color="grey", lw=.6)
    fig.tight_layout(); fig.savefig(f"{FIG}/Q1_PCoA_BrayCurtis.png", dpi=150); plt.close(fig)
    return D_bc


def q1_clustering(P, L, grp, D_bc):
    # (a) Ward dendrogram on Pareto-scaled log data (Euclidean)
    Z_eu = linkage(P.values, method="ward")
    # (b) dendrogram on Bray-Curtis distances (average linkage)
    Z_bc = linkage(squareform(D_bc, checks=False), method="average")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, Z, ttl in [(axes[0], Z_eu, "Ward / Euclidean (log10+Pareto)"),
                       (axes[1], Z_bc, "Average / Bray-Curtis")]:
        dn = dendrogram(Z, labels=list(P.index), ax=ax, leaf_rotation=90)
        ax.set_title(f"Hierarchical clustering\n{ttl}")
        for lbl in ax.get_xmajorticklabels():
            samp = lbl.get_text()
            g = grp[list(P.index).index(samp)]
            lbl.set_color(GROUP_COLORS[g])
    fig.tight_layout(); fig.savefig(f"{FIG}/Q1_dendrograms.png", dpi=150); plt.close(fig)

    # cluster agreement with true labels (2 clusters from Ward)
    cl = fcluster(Z_eu, t=2, criterion="maxclust")
    # map clusters to groups by majority and compute accuracy
    agree = max(
        (cl == c).astype(int).tolist() == (grp == g).astype(int).tolist()
        for c in np.unique(cl) for g in np.unique(grp)
    )
    summary["ward_two_clusters"] = cl.tolist()
    summary["ward_recovers_groups"] = bool(
        len(set(map(tuple, [tuple(cl[grp == g]) for g in np.unique(grp)]))) == 2
        and all(len(set(cl[grp == g])) == 1 for g in np.unique(grp))
    )

    # (c) clustered heatmap of top-variance features (z-scored)
    var = L.var(axis=0).sort_values(ascending=False)
    top = var.index[:50]
    Z = ((L[top] - L[top].mean()) / L[top].std()).T   # features x samples
    col_colors = [GROUP_COLORS[g] for g in grp]
    g = sns.clustermap(Z, cmap="RdBu_r", center=0, figsize=(9, 11),
                       col_colors=col_colors, xticklabels=L.index,
                       yticklabels=False, cbar_kws={"label": "z-score (log10)"})
    g.ax_heatmap.set_title("Top-50 variance features (z-scored log10)", pad=60)
    g.savefig(f"{FIG}/Q1_heatmap_top50.png", dpi=150); plt.close("all")


# ----------------------------------------------------------------------------
# Q2 — normality
# ----------------------------------------------------------------------------
def q2_normality(R, L, grp):
    # aggregate: histograms raw vs log
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    raw = R.values.flatten()
    logv = L.values.flatten()
    axes[0, 0].hist(np.log10(raw), bins=80, color="#999")
    axes[0, 0].set_title("All values (x = log10 intensity)\nbimodal: near-imputation-floor mode + real-signal mode")
    axes[0, 0].set_xlabel("log10(relative abundance)")

    # pooled within-group standardized residuals -> QQ plots
    def resid(M):
        out = []
        for g in np.unique(grp):
            sub = M[grp == g]
            r = sub - sub.mean(axis=0)
            s = sub.std(axis=0, ddof=1)
            s[s == 0] = np.nan
            out.append((r / s).values.flatten())
        v = np.concatenate(out)
        return v[~np.isnan(v)]

    res_raw = resid(R)
    res_log = resid(L)
    stats.probplot(res_raw, dist="norm", plot=axes[0, 1])
    axes[0, 1].set_title("QQ — RAW within-group residuals\n(tails flatten: n=4 mechanically bounds residuals)")
    stats.probplot(res_log, dist="norm", plot=axes[1, 0])
    axes[1, 0].set_title("QQ — LOG10 within-group residuals\n(closer to line in the centre than raw)")

    # per-feature skewness raw vs log
    sk_raw = stats.skew(R.values, axis=0)
    sk_log = stats.skew(L.values, axis=0)
    axes[1, 1].hist(sk_raw, bins=60, alpha=.6, label="raw", color="#D55E00")
    axes[1, 1].hist(sk_log, bins=60, alpha=.6, label="log10", color="#0072B2")
    axes[1, 1].axvline(0, color="k", lw=.8)
    axes[1, 1].set_title("Per-feature skewness")
    axes[1, 1].set_xlabel("skewness"); axes[1, 1].legend()
    fig.tight_layout(); fig.savefig(f"{FIG}/Q2_normality.png", dpi=150); plt.close(fig)

    # Shapiro on pooled residuals (aggregate normality read)
    sh_raw = stats.shapiro(res_raw[:5000] if len(res_raw) > 5000 else res_raw)
    sh_log = stats.shapiro(res_log[:5000] if len(res_log) > 5000 else res_log)

    # per-feature Shapiro within each group on LOG data (n=4 -> low power; descriptive)
    def frac_normal(M):
        cnt = 0
        tot = 0
        for g in np.unique(grp):
            sub = M[grp == g].values
            for j in range(sub.shape[1]):
                col = sub[:, j]
                if np.ptp(col) == 0:
                    continue
                try:
                    _, p = stats.shapiro(col)
                    cnt += int(p > 0.05); tot += 1
                except Exception:
                    pass
        return cnt / tot if tot else float("nan"), tot

    fn_log, tot_log = frac_normal(L)
    fn_raw, tot_raw = frac_normal(R)

    summary["normality"] = {
        "median_skew_raw": float(np.median(sk_raw)),
        "median_skew_log": float(np.median(sk_log)),
        "shapiro_pooled_raw_resid_p": float(sh_raw.pvalue),
        "shapiro_pooled_log_resid_p": float(sh_log.pvalue),
        "perfeature_shapiro_frac_pass_raw": round(float(fn_raw), 3),
        "perfeature_shapiro_frac_pass_log": round(float(fn_log), 3),
        "perfeature_shapiro_n_tests": int(tot_log),
        "note": "n=4/group -> per-feature Shapiro is severely underpowered; descriptive only.",
    }
    print(f"[normality] median skew raw={np.median(sk_raw):.2f} log={np.median(sk_log):.2f} | "
          f"Shapiro pooled resid p raw={sh_raw.pvalue:.2e} log={sh_log.pvalue:.2e} | "
          f"per-feat Shapiro pass raw={fn_raw:.2f} log={fn_log:.2f}")


# ----------------------------------------------------------------------------
# Q3 — differential features + volcano
# ----------------------------------------------------------------------------
def q3_differential(R, L, grp, feat_info):
    g0, g1 = "Starch", "Glucose"   # FC direction: log2(Glucose / Starch)
    m0 = grp == g0
    m1 = grp == g1

    mean0 = R.values[m0].mean(axis=0)   # mean relative abundance Starch
    mean1 = R.values[m1].mean(axis=0)   # mean relative abundance Glucose
    log2fc_arith = np.log2(mean1 / mean0)            # ratio of arithmetic means

    # geometric-mean log2FC = difference of mean-log10s, rescaled to log2.
    # This is the effect size the Welch-on-log test actually evaluates, so it is
    # the consistent x-axis for the volcano.
    mlog0 = L.values[m0].mean(axis=0)
    mlog1 = L.values[m1].mean(axis=0)
    log2fc = (mlog1 - mlog0) / np.log10(2.0)

    # Welch t-test on log10 values (data are log-normal)
    t_stat, t_p = stats.ttest_ind(L.values[m1], L.values[m0],
                                  axis=0, equal_var=False)
    t_fdr = bh_fdr(t_p)

    # Mann-Whitney U (non-parametric, exact for n small)
    mw_p = np.empty(L.shape[1])
    for j in range(L.shape[1]):
        try:
            mw_p[j] = stats.mannwhitneyu(L.values[m1, j], L.values[m0, j],
                                         alternative="two-sided").pvalue
        except ValueError:
            mw_p[j] = 1.0
    mw_fdr = bh_fdr(mw_p)

    res = feat_info.copy()
    res["mean_Starch"] = mean0
    res["mean_Glucose"] = mean1
    res["log2FC_Glucose_vs_Starch"] = log2fc          # geometric (test-consistent)
    res["log2FC_arith_mean_ratio"] = log2fc_arith     # arithmetic-mean ratio (reference)
    res["t_stat"] = t_stat
    res["t_p"] = t_p
    res["t_FDR"] = t_fdr
    res["mw_p"] = mw_p
    res["mw_FDR"] = mw_fdr
    res = res.sort_values("t_p").reset_index(drop=True)

    FC_THR, FDR_THR = 1.0, 0.05
    res["significant"] = (res["t_FDR"] < FDR_THR) & (res["log2FC_Glucose_vs_Starch"].abs() >= FC_THR)
    res.drop(columns=["col"]).to_csv(f"{TAB}/differential_features.csv", index=False)

    n_t_fdr = int((t_fdr < FDR_THR).sum())
    n_t_sig = int(res["significant"].sum())
    n_mw_fdr = int((mw_fdr < FDR_THR).sum())
    n_t_raw = int((t_p < 0.05).sum())
    n_mw_raw = int((mw_p < 0.05).sum())
    summary["differential"] = {
        "fc_direction": "log2 geometric-mean ratio (Glucose/Starch); + = higher in Glucose",
        "n_features": int(L.shape[1]),
        "welch_raw_p_lt_0.05": n_t_raw,
        "welch_FDR_lt_0.05": n_t_fdr,
        "welch_FDR_and_abs_log2FC_ge_1": n_t_sig,
        "mannwhitney_raw_p_lt_0.05": n_mw_raw,
        "mannwhitney_FDR_lt_0.05": n_mw_fdr,
        "mannwhitney_min_possible_p": round(2.0 / 70, 4),
        "up_in_glucose_sig": int(((res["significant"]) & (res["log2FC_Glucose_vs_Starch"] > 0)).sum()),
        "up_in_starch_sig": int(((res["significant"]) & (res["log2FC_Glucose_vs_Starch"] < 0)).sum()),
        "fc_thr": FC_THR, "fdr_thr": FDR_THR,
    }
    print(f"[differential] Welch raw p<0.05: {n_t_raw} | FDR<0.05: {n_t_fdr} | "
          f"FDR & |log2FC|>=1: {n_t_sig}  || MW raw p<0.05: {n_mw_raw} | FDR<0.05: {n_mw_fdr}")

    # volcano
    fig, ax = plt.subplots(figsize=(11, 8.5))
    x = res["log2FC_Glucose_vs_Starch"].values
    y = -np.log10(res["t_p"].values)
    sig = res["significant"].values
    ax.scatter(x[~sig], y[~sig], s=12, c="#cccccc", alpha=.5,
               label=f"not significant (n={int((~sig).sum())})", rasterized=True)
    up = sig & (x > 0)
    dn = sig & (x < 0)
    ax.scatter(x[dn], y[dn], s=20, c="#D55E00", alpha=.75, edgecolor="none",
               label=f"higher in Starch (n={int(dn.sum())})")
    ax.scatter(x[up], y[up], s=20, c="#0072B2", alpha=.75, edgecolor="none",
               label=f"higher in Glucose (n={int(up.sum())})")
    ax.axvline(FC_THR, ls="--", c="grey", lw=.8)
    ax.axvline(-FC_THR, ls="--", c="grey", lw=.8)
    ax.axvline(0, c="grey", lw=.5)
    passing = res["t_p"][res["t_FDR"] < FDR_THR]
    if len(passing):
        p_line = passing.max()
        ax.axhline(-np.log10(p_line), ls="--", c="green", lw=.9,
                   label=f"FDR=0.05 threshold (raw p≈{p_line:.1e})")
    # annotate the 12 most significant features by feature_id
    for _, r in res.head(12).iterrows():
        ax.annotate(str(r["feature_id"]),
                    (r["log2FC_Glucose_vs_Starch"], -np.log10(r["t_p"])),
                    fontsize=7, xytext=(3, 3), textcoords="offset points",
                    color="black")
    ax.set_xlabel("log2 fold change  (Glucose / Starch)  — geometric-mean ratio")
    ax.set_ylabel("-log10(raw p, Welch t-test on log10)")
    ax.set_title(f"Volcano plot — all {L.shape[1]} features\n"
                 f"colored = BH-FDR < {FDR_THR} AND |log2FC| ≥ {FC_THR}")
    ax.legend(loc="upper center", fontsize=9, ncol=2, framealpha=.9)
    fig.tight_layout()
    fig.savefig(f"{FIG}/Q3_volcano.png", dpi=150)
    plt.close(fig)

    # top-25 table preview
    top = res.head(25)[["feature_id", "mz", "rt", "log2FC_Glucose_vs_Starch",
                        "t_p", "t_FDR", "mw_p", "significant"]]
    top.to_csv(f"{TAB}/top25_features.csv", index=False)
    return res


# ----------------------------------------------------------------------------
def main():
    X, grp, feat_info = load()
    R, L, P, A = preprocess(X)
    print("\n=== Q1: group separation ===")
    q1_pca(P, grp)
    D_bc = q1_pcoa_permanova(R, L, grp)
    q1_clustering(P, L, grp, D_bc)
    print("\n=== Q2: normality ===")
    q2_normality(R, L, grp)
    print("\n=== Q3: differential features ===")
    q3_differential(R, L, grp, feat_info)

    with open(f"{TAB}/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n[done] figures ->", FIG)
    print("[done] tables  ->", TAB)
    print("\n=== SUMMARY ===")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
