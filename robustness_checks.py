#!/usr/bin/env python3
"""
Robustness & sensitivity companion to metabolite_analysis.py.

Adds the checks raised by an independent statistical review of the
Starch-vs-Glucose analysis:

  1. PERMDISP  - is the PERMANOVA separation a location shift or a dispersion
                 (spread) artifact?
  2. Leave-one-out PERMANOVA - does dropping the within-group outlier
                 (EXT_wo_ami_1) change Q1?
  3. CLR / Aitchison sensitivity - does the compositional (sum=1) constraint
                 change PCA / PERMANOVA / differential conclusions?
  4. Moderated t-test (limma-style empirical-Bayes variance shrinkage) - a
                 stable per-feature ranking for n=4, vs the unstable Welch t.
  5. Silhouette + per-sample QC + imputation-floor flags + an actionable
                 shortlist + BY-FDR sensitivity.

Same input and preprocessing as the primary script (consistent matrices).
"""
import os, re, json, itertools, warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.special import psi as digamma, polygamma
from scipy.optimize import brentq
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
CSV = os.path.join(ROOT, "2026-06-06T21-12_export.csv")
FIG = os.path.join(ROOT, "results", "figures"); TAB = os.path.join(ROOT, "results", "tables")
GC = {"Starch": "#D55E00", "Glucose": "#0072B2"}
out = {}


def trigamma(x):
    return polygamma(1, x)


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    q = np.empty(n); q[o] = np.clip(r, 0, 1); return q


def by_fdr(p):
    """Benjamini-Yekutieli (valid under arbitrary dependence)."""
    n = len(p); c = np.sum(1.0 / np.arange(1, n + 1))
    return np.clip(bh_fdr(p) * c, 0, 1)


def load():
    df = pd.read_csv(CSV).rename(columns=lambda c: c)
    df = df.rename(columns={df.columns[0]: "sample"})
    meta = ["sample", "ATTRIBUTE_Sample", "ATTRIBUTE_Blank",
            "ATTRIBUTE_Blank_Ctrl", "ATTRIBUTE_CarbonSource"]
    feat = [c for c in df.columns if c not in meta]
    grp = df["ATTRIBUTE_CarbonSource"].replace({"Startch": "Starch"}).to_numpy()
    X = df[feat].apply(pd.to_numeric, errors="coerce"); X.index = df["sample"].to_numpy()
    info = []
    for c in feat:
        m = re.match(r"^(\d+)_([\d.]+)_([\d.]+)", c)
        info.append((c, m.group(1), float(m.group(2)), float(m.group(3))))
    info = pd.DataFrame(info, columns=["col", "feature_id", "mz", "rt"])
    return X, grp, info


# ---- PERMANOVA (exact) + PERMDISP -----------------------------------------
def permanova(D, grp):
    # integer group codes avoid numpy string-dtype truncation during permutation
    n = D.shape[0]; labs = np.unique(grp); a = len(labs)
    code = np.where(grp == labs[0], 0, 1)
    iu = np.triu_indices(n, 1); d2 = D[iu] ** 2

    def ss(gv):
        sst = d2.sum() / n; ssw = 0.0
        for l in np.unique(gv):
            mem = np.where(gv == l)[0]
            if len(mem) < 2: continue
            sub = D[np.ix_(mem, mem)]; s = np.triu_indices(len(mem), 1)
            ssw += (sub[s] ** 2).sum() / len(mem)
        return sst - ssw, ssw, sst
    ssa, ssw, sst = ss(code); F = (ssa / (a - 1)) / (ssw / (n - a))
    n1 = int((code == 0).sum())
    Fp = []
    for c in itertools.combinations(range(n), n1):
        gv = np.ones(n, int); gv[list(c)] = 0
        s = ss(gv); Fp.append((s[0] / (a - 1)) / (s[1] / (n - a)))
    Fp = np.array(Fp); p = (Fp >= F - 1e-12).mean()
    return dict(F=float(F), R2=float(ssa / sst), p=float(p),
                n_perm=len(Fp), is_max=bool(F >= Fp.max() - 1e-9))


def pcoa_coords(D):
    n = D.shape[0]; A = -0.5 * D ** 2; J = np.eye(n) - np.ones((n, n)) / n
    B = J @ A @ J; w, v = np.linalg.eigh(B); idx = np.argsort(w)[::-1]
    w, v = w[idx], v[:, idx]; pos = w > 1e-9
    return v[:, pos] * np.sqrt(w[pos])


def permdisp(D, grp):
    """Anderson PERMDISP: test homogeneity of within-group dispersions.
    Distances-to-centroid are fixed; group labels are permuted (integer codes)."""
    coords = pcoa_coords(D)
    labs = np.unique(grp); a = len(labs); n = len(grp)
    dist = np.zeros(n)
    for l in labs:
        m = grp == l; dist[m] = np.sqrt(((coords[m] - coords[m].mean(0)) ** 2).sum(1))
    code = np.where(grp == labs[0], 0, 1)
    def F(d, gv):
        gm = d.mean(); cats = np.unique(gv)
        ssb = sum((gv == c).sum() * (d[gv == c].mean() - gm) ** 2 for c in cats)
        ssw = sum(((d[gv == c] - d[gv == c].mean()) ** 2).sum() for c in cats)
        return (ssb / (a - 1)) / (ssw / (n - a))
    Fobs = F(dist, code); n1 = int((code == 0).sum())
    Fp = []
    for c in itertools.combinations(range(n), n1):
        gv = np.ones(n, int); gv[list(c)] = 0; Fp.append(F(dist, gv))
    Fp = np.array(Fp); p = (Fp >= Fobs - 1e-12).mean()
    means = {l: float(dist[grp == l].mean()) for l in labs}
    return dict(F=float(Fobs), p=float(p), group_dispersion=means)


# ---- moderated t-test (limma / Smyth 2004) --------------------------------
def moderated_t(L, grp):
    g0, g1 = "Starch", "Glucose"; m0, m1 = grp == g0, grp == g1
    n0, n1 = m0.sum(), m1.sum(); d = n0 + n1 - 2
    x0, x1 = L.values[m0], L.values[m1]
    diff = x1.mean(0) - x0.mean(0)
    ss = ((x0 - x0.mean(0)) ** 2).sum(0) + ((x1 - x1.mean(0)) ** 2).sum(0)
    s2 = ss / d
    ok = s2 > 0
    z = np.log(s2[ok])
    e = z - digamma(d / 2) + np.log(d / 2)
    ebar = e.mean(); G = ok.sum()
    target = np.mean((e - ebar) ** 2) * G / (G - 1) - trigamma(d / 2)
    if target <= 0:
        d0 = np.inf; s02 = np.exp(ebar)
    else:
        d0 = float(brentq(lambda x: trigamma(x / 2) - target, 1e-3, 1e6))
        s02 = float(np.exp(ebar + digamma(d0 / 2) - np.log(d0 / 2)))
    if np.isinf(d0):
        s2_post = np.full(L.shape[1], s02); df_post = np.inf
    else:
        s2_post = (d0 * s02 + d * s2) / (d0 + d); df_post = d0 + d
    se = np.sqrt(s2_post * (1 / n0 + 1 / n1))
    t = np.where(se > 0, diff / se, 0.0)
    if np.isinf(df_post):
        p = 2 * stats.norm.sf(np.abs(t))
    else:
        p = 2 * stats.t.sf(np.abs(t), df_post)
    return diff, t, p, dict(prior_df=float(d0), s0=float(np.sqrt(s02)),
                            posterior_df=(None if np.isinf(df_post) else float(df_post)))


def main():
    X, grp, info = load()
    R = X.copy(); L = np.log10(R)
    D_bc = squareform(pdist(R.values, metric="braycurtis"))

    # 1) PERMDISP
    pd_res = permdisp(D_bc, grp); out["permdisp"] = pd_res
    print(f"[PERMDISP] F={pd_res['F']:.3f} p={pd_res['p']:.3f} "
          f"dispersion={ {k: round(v,3) for k,v in pd_res['group_dispersion'].items()} }")

    # 2) leave-one-out PERMANOVA (drop outlier EXT_wo_ami_1)
    out["loo_permanova"] = {}
    for drop in ["EXT_wo_ami_1"]:
        keep = X.index != drop
        Dk = squareform(pdist(R.values[keep], metric="braycurtis"))
        r = permanova(Dk, grp[keep]); out["loo_permanova"][f"drop_{drop}"] = r
        print(f"[LOO PERMANOVA | drop {drop}] F={r['F']:.2f} R2={r['R2']:.3f} "
              f"p={r['p']:.4f} (n_perm={r['n_perm']})")
    # silhouette on Bray-Curtis
    sil = float(silhouette_score(D_bc, grp, metric="precomputed"))
    out["silhouette_braycurtis"] = sil
    print(f"[silhouette | Bray-Curtis] {sil:.3f}")

    # 3) CLR / Aitchison sensitivity
    logR = np.log(R.values)
    clr = logR - logR.mean(axis=1, keepdims=True)        # CLR per sample
    clr = pd.DataFrame(clr, index=R.index, columns=R.columns)
    D_ait = squareform(pdist(clr.values, metric="euclidean"))
    perm_ait = permanova(D_ait, grp); out["aitchison_permanova"] = perm_ait
    # CLR PCA (Pareto) vs log10 PCA (Pareto)
    def pareto_pca(M):
        mu = M.mean(0); sd = M.std(0, ddof=1).replace(0, np.nan)
        P = ((M - mu) / np.sqrt(sd)).fillna(0.0)
        return PCA(n_components=2).fit(P.values).explained_variance_ratio_
    ev_log = pareto_pca(L); ev_clr = pareto_pca(clr)
    out["pca_pc1_log10_vs_clr"] = [round(float(ev_log[0]) * 100, 2),
                                   round(float(ev_clr[0]) * 100, 2)]
    print(f"[Aitchison PERMANOVA] F={perm_ait['F']:.2f} R2={perm_ait['R2']:.3f} "
          f"p={perm_ait['p']:.4f}  | PCA PC1 log10={ev_log[0]*100:.1f}% clr={ev_clr[0]*100:.1f}%")

    # differential concordance: Welch on log10 vs Welch on CLR
    m0, m1 = grp == "Starch", grp == "Glucose"
    _, p_log = stats.ttest_ind(L.values[m1], L.values[m0], axis=0, equal_var=False)
    _, p_clr = stats.ttest_ind(clr.values[m1], clr.values[m0], axis=0, equal_var=False)
    q_log = bh_fdr(p_log); q_clr = bh_fdr(p_clr)
    s_log = set(np.where(q_log < 0.05)[0]); s_clr = set(np.where(q_clr < 0.05)[0])
    jac = len(s_log & s_clr) / len(s_log | s_clr)
    rho = stats.spearmanr(p_log, p_clr).statistic
    out["clr_concordance"] = {"jaccard_sig_sets": round(jac, 3),
                              "spearman_pvalues": round(float(rho), 3),
                              "n_sig_log10": len(s_log), "n_sig_clr": len(s_clr)}
    print(f"[CLR concordance] Jaccard(sig sets)={jac:.3f} Spearman(p)={rho:.3f} "
          f"(log10 {len(s_log)} vs CLR {len(s_clr)} sig)")

    # 4) moderated t + annotate the differential table
    diff, mt, mp, mod_meta = moderated_t(L, grp)
    out["moderated_t"] = mod_meta
    mod_q = bh_fdr(mp)
    n_mod = int((mod_q < 0.05).sum())
    print(f"[moderated t] prior_df={mod_meta['prior_df']:.2f} "
          f"posterior_df={mod_meta['posterior_df']} | FDR<0.05={n_mod}")

    # floor flag + abundance + within-group CV, build annotated + shortlist
    floor = R.values.min()
    near = R.values <= floor * 10
    floor_any = (near[m0].any(0) | near[m1].any(0))
    log2fc = (L.values[m1].mean(0) - L.values[m0].mean(0)) / np.log10(2)
    mean_min = np.minimum(R.values[m0].mean(0), R.values[m1].mean(0))
    cv = lambda M: M.std(0, ddof=1) / M.mean(0)
    maxcv = np.maximum(cv(R.values[m0]), cv(R.values[m1]))

    ann = info.drop(columns=["col"]).copy()
    ann["log2FC_Glucose_vs_Starch"] = log2fc
    ann["mean_min_group"] = mean_min
    ann["max_within_group_CV"] = maxcv
    ann["welch_FDR"] = q_log
    ann["welch_BY_FDR"] = by_fdr(p_log)
    ann["moderated_t"] = mt
    ann["moderated_p"] = mp
    ann["moderated_FDR"] = mod_q
    ann["near_imputation_floor"] = floor_any
    ann = ann.sort_values("moderated_p").reset_index(drop=True)
    ann.to_csv(f"{TAB}/differential_features_annotated.csv", index=False)

    # actionable shortlist: stable + strong + not floor-driven + above detection
    short = ann[(ann.moderated_FDR < 0.05) & (ann.log2FC_Glucose_vs_Starch.abs() >= 1)
                & (~ann.near_imputation_floor) & (ann.mean_min_group >= 1e-5)
                & (ann.max_within_group_CV <= 0.5)].copy()
    short.to_csv(f"{TAB}/shortlist_robust_features.csv", index=False)
    out["counts"] = {
        "welch_BH_FDR_lt_0.05": int((q_log < 0.05).sum()),
        "welch_BY_FDR_lt_0.05": int((by_fdr(p_log) < 0.05).sum()),
        "moderated_BH_FDR_lt_0.05": n_mod,
        "shortlist_robust": int(len(short)),
        "n_features_touching_floor": int(floor_any.sum()),
    }
    print(f"[counts] Welch BH={out['counts']['welch_BH_FDR_lt_0.05']} "
          f"BY={out['counts']['welch_BY_FDR_lt_0.05']} "
          f"moderated={n_mod} | robust shortlist={len(short)} "
          f"| features touching floor={int(floor_any.sum())}")

    # ---- figure: robustness panel
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    # (a) log10 vs CLR -log10 p
    ax[0].scatter(-np.log10(p_log), -np.log10(p_clr), s=6, alpha=.3, c="#444")
    lim = max(ax[0].get_xlim()[1], ax[0].get_ylim()[1])
    ax[0].plot([0, lim], [0, lim], "r--", lw=.8)
    ax[0].set_xlabel("-log10 p (log10 transform)"); ax[0].set_ylabel("-log10 p (CLR)")
    ax[0].set_title(f"Differential p concordance\nSpearman={rho:.2f}, Jaccard={jac:.2f}")
    # (b) PERMDISP distances to centroid
    coords = pcoa_coords(D_bc); dist = np.zeros(len(grp))
    for l in np.unique(grp):
        mm = grp == l; dist[mm] = np.sqrt(((coords[mm] - coords[mm].mean(0)) ** 2).sum(1))
    for l in np.unique(grp):
        ax[1].scatter(np.full((grp == l).sum(), l), dist[grp == l], s=80,
                      color=GC[l], edgecolor="k")
    ax[1].set_ylabel("distance to group centroid (PCoA)")
    ax[1].set_title(f"PERMDISP (dispersion)\nF={pd_res['F']:.2f}, p={pd_res['p']:.2f} (n.s. = no spread artifact)")
    # (c) variance shrinkage: raw vs moderated SE proxy
    ax[2].scatter(diff, mt, s=6, alpha=.3, c="#0072B2")
    ax[2].axhline(0, color="grey", lw=.6); ax[2].axvline(0, color="grey", lw=.6)
    ax[2].set_xlabel("mean log10 difference (Glucose - Starch)")
    ax[2].set_ylabel("moderated t-statistic")
    ax[2].set_title(f"Moderated (empirical-Bayes) t\nprior df={mod_meta['prior_df']:.1f}")
    fig.tight_layout(); fig.savefig(f"{FIG}/Q3b_robustness.png", dpi=150); plt.close(fig)

    # QC table
    qc = pd.DataFrame({
        "sample": X.index, "group": grp,
        "row_sum": R.sum(1).values,
        "n_features_at_floor": near.sum(1),
        "max_feature_abundance": R.max(1).values,
        "median_log10": np.median(L.values, 1),
    })
    qc.to_csv(f"{TAB}/qc_per_sample.csv", index=False)

    with open(f"{TAB}/summary_robustness.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\n[done robustness] ->", json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
