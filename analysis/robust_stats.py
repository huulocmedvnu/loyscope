# Statistical robustness of the main results. Inputs are the per-group fits written by loyfit (f, bootstrap 95% interval,
# mu, mean UMI); this script only combines them. Outputs: results/tables/robust_*.tsv and robust_summary.json.
#  1. age slopes: inverse-variance weighted least squares of donor values on age (+ disease when present),
#     two-sided permutation P (age permuted within disease strata), Benjamini-Hochberg q within dataset
#  2. blood meta-analysis: REML between-cohort variance with Hartung-Knapp confidence intervals
#  3. all disease contrasts (every cell type with >= 4 donors per arm), weighted, permutation P, BH q across all
#  4. sensitivity of key results to the expected-Y threshold (2, 3, 5, 10)
#  5. overdispersion: share of donors above 10% compared with a null of one common fraction and the observed precision
import csv, os, re, json, math, random
import numpy as np

RES = os.environ.get("LOY_RESULTS", "results"); EXT = os.path.join(RES, "extract"); TAB = os.path.join(RES, "tables")
os.makedirs(TAB, exist_ok=True)
DEC = {"third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75, "ninth": 85}
def age(s):
    m = re.match(r"(\d+)-year-old", s)
    if m: return float(m.group(1))
    m = re.match(r"(\d+)-(\d+) year-old", s)
    if m: return (float(m.group(1)) + float(m.group(2))) / 2
    m = re.match(r"(\w+) decade", s); return DEC.get(m.group(1)) if m else None
def rows(fn): return list(csv.DictReader(open(fn), delimiter="\t"))
SE_FLOOR = 0.005  # minimum standard error of a group estimate (0.5 points), avoids infinite weights at the boundary

def donor_values(fn, thr, sel=lambda r: True):
    """inverse-variance combination of group estimates per (cell type, donor); value = f(Y) - median f(control)"""
    out = {}
    for r in rows(fn):
        if r["sex"] != "male" or r["k"] == "NaN" or float(r["mu"]) * float(r["mean_umi"]) < thr or not sel(r): continue
        if not (r["assay"].startswith("10x") or "GEO" in r["assay"] or "ENCODE" in r["assay"]): continue
        v = float(r["f"]) - float(r["pseudo_f_med"])
        lo, hi = float(r["f_lo"]), float(r["f_hi"])
        se = max((hi - lo) / 3.92 if math.isfinite(lo) and math.isfinite(hi) else 0.05, SE_FLOOR)
        out.setdefault((r["cell_type"], r["donor"]), []).append((v, se, r["disease"], age(r["stage"]), int(r["n"])))
    res = {}
    for k, g in out.items():
        w = [1 / s ** 2 for _, s, _, _, _ in g]; v = sum(a * x[0] for a, x in zip(w, g)) / sum(w)
        res[k] = dict(v=v, se=math.sqrt(1 / sum(w)), disease=g[0][2], age=g[0][3], n=sum(x[4] for x in g))
    return res

def wls(X, y, w):
    W = np.sqrt(w); b, *_ = np.linalg.lstsq(X * W[:, None], y * W, rcond=None); return b
def re_weights(X, y, se):
    """random-effects weights 1/(se^2 + tau^2); tau^2 by the Paule-Mandel estimator for meta-regression"""
    n, p = X.shape
    def Q(t):
        w = 1 / (se ** 2 + t); b = wls(X, y, w); return (w * (y - X @ b) ** 2).sum()
    if Q(0.0) <= n - p: return 1 / se ** 2, 0.0
    lo, hi = 0.0, max(1e-6, float(np.var(y)) * 10)
    while Q(hi) > n - p: hi *= 4
    for _ in range(80):
        mid = (lo + hi) / 2
        if Q(mid) > n - p: lo = mid
        else: hi = mid
    t = (lo + hi) / 2; return 1 / (se ** 2 + t), t
def perm_p(X, y, w, j, strata, nperm=2000, seed=1):
    """two-sided permutation P for coefficient j, column j permuted within strata"""
    b0 = wls(X, y, w)[j]; rng = random.Random(seed); ge = 0
    idx = {}
    for i, s in enumerate(strata): idx.setdefault(s, []).append(i)
    for _ in range(nperm):
        Xp = X.copy(); col = X[:, j].copy()
        for s, ii in idx.items():
            vals = [col[i] for i in ii]; rng.shuffle(vals)
            for i, v in zip(ii, vals): Xp[i, j] = v
        ge += abs(wls(Xp, y, w)[j]) >= abs(b0) - 1e-12
    return b0, (ge + 1) / (nperm + 1)
def bh(ps):
    m = len(ps); order = sorted(range(m), key=lambda i: ps[i]); q = [0.0] * m; run = 1.0
    for rank, i in reversed(list(enumerate(order, 1))): run = min(run, ps[i] * m / rank); q[i] = run
    return q

HEART = lambda r: any(k in r["tissue"].lower() for k in ("heart", "ventric", "atri", "myocard"))
def age_slopes(name, thr=5, min_donors=8, unweighted=False):
    dv = donor_values(os.path.join(EXT, f"{name}.fitamb.tsv"), thr, HEART if name == "encode_heart" else (lambda r: True))
    by = {}
    for (ct, d), x in dv.items():
        if x["age"] is not None: by.setdefault(ct, []).append(x)
    out = []
    for ct, g in by.items():
        if len(g) < min_donors or max(x["age"] for x in g) - min(x["age"] for x in g) < 20: continue
        dis = sorted({x["disease"] for x in g}); y = np.array([x["v"] for x in g]); a = np.array([x["age"] for x in g]) / 10
        cols = [np.ones(len(g)), a] + [np.array([1.0 if x["disease"] == d else 0.0 for x in g]) for d in dis[1:]]
        X = np.column_stack(cols); w = np.ones(len(g)) if unweighted else re_weights(X, y, np.array([x["se"] for x in g]))[0]
        b, p = perm_p(X, y, w, 1, [x["disease"] for x in g])
        out.append(dict(dataset=name, cell_type=ct, donors=len(g), slope_per_decade=100 * b, p=p, weighted=not unweighted, adjusted_for_disease=len(dis) > 1))
    q = bh([o["p"] for o in out])
    for o, qq in zip(out, q): o["q"] = qq
    return out

# ---------- 1. age slopes ----------
slopes = []
for ds in ["ad_mssm", "ad_radc", "retina", "aida", "iaging", "onek1k_geo", "blood", "tsapiens", "ad_microglia", "encode_heart"]:
    if os.path.exists(os.path.join(EXT, f"{ds}.fitamb.tsv")):
        thr = 1.5 if ds in ("iaging", "onek1k_geo", "blood") else 2 if ds == "encode_heart" else 5
        slopes += age_slopes(ds, thr)
unw = {}
for ds in ["ad_mssm", "ad_radc", "retina", "aida", "iaging", "onek1k_geo", "blood", "tsapiens", "ad_microglia", "encode_heart"]:
    if os.path.exists(os.path.join(EXT, f"{ds}.fitamb.tsv")):
        for o in age_slopes(ds, 1.5 if ds in ("iaging", "onek1k_geo", "blood") else 2 if ds == "encode_heart" else 5, unweighted=True): unw[(ds, o["cell_type"])] = o
with open(os.path.join(TAB, "robust_age_slopes.tsv"), "w") as f:
    f.write("dataset\tcell_type\tdonors\tslope_per_decade\tp\tq\tadjusted_for_disease\tunweighted_slope\tunweighted_p\n")
    for o in sorted(slopes, key=lambda o: (o["dataset"], o["p"])):
        u = unw.get((o["dataset"], o["cell_type"]), {})
        f.write(f"{o['dataset']}\t{o['cell_type']}\t{o['donors']}\t{o['slope_per_decade']:.3f}\t{o['p']:.4g}\t{o['q']:.4g}\t{o['adjusted_for_disease']}\t{u.get('slope_per_decade', float('nan')):.3f}\t{u.get('p', float('nan')):.4g}\n")

# ---------- 2. blood meta-analysis, REML + Hartung-Knapp ----------
CLS = [("classical monocyte", lambda c: ("CD14-positive" in c and "CD16-positive" not in c and "monocyte" in c) or c == "classical monocyte"),
       ("non-classical monocyte", lambda c: "CD16-positive" in c and "monocyte" in c), ("NK cell", lambda c: "natural killer" in c and "T cell" not in c),
       ("regulatory T cell", lambda c: "regulatory T" in c), ("naive B cell", lambda c: c == "naive B cell"), ("memory B cell", lambda c: "memory B cell" in c),
       ("CD4 T cell", lambda c: "CD4-positive" in c and "regulatory" not in c), ("CD8 T cell", lambda c: "CD8-positive" in c and "regulatory" not in c)]
def cohort_slope(name, fnc, thr=1.5):
    dv = donor_values(os.path.join(EXT, f"{name}.fitamb.tsv"), thr, lambda r: r["disease"] == "normal" and fnc(r["cell_type"]))
    dd = {}
    for (ct, d), x in dv.items():
        if x["age"] is None: continue
        e = dd.setdefault(d, [0.0, 0.0, x["age"]]); w = 1 / x["se"] ** 2; e[0] += w * x["v"]; e[1] += w
    pts = [(e[2] / 10, e[0] / e[1], e[1]) for e in dd.values()]
    if len(pts) < 10: return None
    x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts]); se0 = 1 / np.sqrt(np.array([p[2] for p in pts]))
    X = np.c_[np.ones_like(x), x]; w = re_weights(X, y, se0)[0]; b = wls(X, y, w)
    r = y - X @ b; s2 = (w * r ** 2).sum() / (len(y) - 2); cov = s2 * np.linalg.inv(X.T @ (X * w[:, None]))
    return b[1], math.sqrt(cov[1, 1]), len(pts)
def reml_hk(est):
    y = np.array([e[0] for e in est]); v = np.array([e[1] ** 2 for e in est]); k = len(y); tau2 = 0.0
    for _ in range(200):
        w = 1 / (v + tau2); mu = (w * y).sum() / w.sum()
        tau2_new = max(0.0, ((w ** 2) * ((y - mu) ** 2 - v)).sum() / (w ** 2).sum() + 1 / w.sum())
        if abs(tau2_new - tau2) < 1e-12: break
        tau2 = tau2_new
    w = 1 / (v + tau2); mu = (w * y).sum() / w.sum()
    q_hk = (w * (y - mu) ** 2).sum() / (k - 1); se = math.sqrt(q_hk / w.sum())
    from math import lgamma
    def t_ppf975(df):  # Student t 97.5% quantile by bisection on the regularised incomplete beta
        def cdf(t):
            x = df / (df + t * t); a, b_ = df / 2, 0.5
            # continued fraction for incomplete beta (Numerical Recipes)
            def betacf(a, b, x):
                qab, qap, qam = a + b, a + 1, a - 1; c, d = 1.0, 1 - qab * x / qap; d = 1 / d; h = d
                for m in range(1, 200):
                    m2 = 2 * m; aa = m * (b - m) * x / ((qam + m2) * (a + m2)); d = 1 + aa * d; c = 1 + aa / c; d = 1 / d; h *= d * c
                    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2)); d = 1 + aa * d; c = 1 + aa / c; d = 1 / d; h *= d * c
                return h
            bt = math.exp(lgamma(a + b_) - lgamma(a) - lgamma(b_) + a * math.log(x) + b_ * math.log(1 - x))
            ib = bt * betacf(a, b_, x) / a if x < (a + 1) / (a + b_ + 2) else 1 - bt * betacf(b_, a, 1 - x) / b_
            return 1 - 0.5 * ib
        lo, hi = 0.0, 100.0
        for _ in range(100):
            mid = (lo + hi) / 2
            if cdf(mid) < 0.975: lo = mid
            else: hi = mid
        return (lo + hi) / 2
    t = t_ppf975(k - 1)
    Q = ((1 / v) * (y - (y / v).sum() / (1 / v).sum()) ** 2).sum(); I2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    # two-sided P from t with k-1 df
    tstat = mu / se
    lo, hi = 0.0, 1.0
    return mu, mu - t * se, mu + t * se, tau2, I2, tstat
meta = []; per_cohort = []
for cl, fnc in CLS:
    est = []
    for n in ("onek1k_geo", "aida", "iaging", "blood"):
        e = cohort_slope(n, fnc)
        if e: est.append(e); per_cohort.append((cl, n, 100 * e[0], 100 * e[1], e[2]))
    if len(est) < 2: continue
    mu, lo, hi, tau2, I2, tstat = reml_hk(est)
    meta.append(dict(cell_class=cl, cohorts=len(est), slope=100 * mu, ci_lo=100 * lo, ci_hi=100 * hi, I2=100 * I2, excludes_zero=lo > 0 or hi < 0))
with open(os.path.join(TAB, "robust_meta.tsv"), "w") as f:
    f.write("cell_class\tcohorts\tslope_per_decade\tci_lo\tci_hi\tI2\texcludes_zero\n")
    for m in meta: f.write(f"{m['cell_class']}\t{m['cohorts']}\t{m['slope']:.2f}\t{m['ci_lo']:.2f}\t{m['ci_hi']:.2f}\t{m['I2']:.0f}\t{m['excludes_zero']}\n")

with open(os.path.join(TAB, "robust_meta_cohorts.tsv"), "w") as f:
    f.write("cell_class\tcohort\tslope_per_decade\tse\tdonors\n")
    for c in per_cohort: f.write("%s\t%s\t%.3f\t%.3f\t%d\n" % c)

# ---------- 3. all disease contrasts ----------
DIS = [("kidney", "chronic kidney disease", 2), ("kidney", "acute kidney injury", 2), ("ckd_sn", "chronic kidney disease", 5), ("ckd_sn", "acute kidney injury", 5),
       ("dcm", "dilated cardiomyopathy", 2), ("ad_microglia", "Alzheimer disease", 2), ("ad_mssm", "Alzheimer disease", 5), ("ad_radc", "Alzheimer disease", 5), ("covid", "COVID-19", 5)]
con = []
for ds, case, thr in DIS:
    fn = os.path.join(EXT, f"{ds}.fitamb.tsv")
    if not os.path.exists(fn): continue
    dv = donor_values(fn, thr, lambda r: r["disease"] in (case, "normal"))
    raw = {}
    for r in rows(fn):
        if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= thr and r["disease"] in (case, "normal"):
            raw.setdefault((r["cell_type"], r["donor"]), []).append((float(r["naive_zero"]), int(r["n"]), float(r["mean_umi"])))
    by = {}
    for (ct, d), x in dv.items():
        if x["age"] is not None: by.setdefault(ct, []).append((d, x))
    for ct, g in by.items():
        nc = sum(1 for _, x in g if x["disease"] == case); nn = len(g) - nc
        if nc < 4 or nn < 4: continue
        y = np.array([x["v"] for _, x in g]); dz = np.array([1.0 if x["disease"] == case else 0.0 for _, x in g]); a = np.array([x["age"] for _, x in g]) / 10
        X = np.column_stack([np.ones(len(g)), dz, a]); w = re_weights(X, y, np.array([x["se"] for _, x in g]))[0]
        b, p = perm_p(X, y, w, 1, ["all"] * len(g))
        yn = np.array([sum(z * n for z, n, _ in raw[(ct, d)]) / sum(n for _, n, _ in raw[(ct, d)]) for d, _ in g])
        # difference in library size between disease and control cells (cell-weighted mean UMI per donor), same model
        yu = np.array([sum(u * n for _, n, u in raw[(ct, d)]) / sum(n for _, n, _ in raw[(ct, d)]) for d, _ in g])
        bu, pu = perm_p(X, yu, np.ones(len(g)), 1, ["all"] * len(g))
        bn, pn = perm_p(X, yn, np.ones(len(g)), 1, ["all"] * len(g))
        # donor bootstrap 95% intervals (weights re-estimated in each resample)
        se_all = np.array([x["se"] for _, x in g]); rb = np.random.default_rng(7); bm, bz = [], []
        for _ in range(500):
            i = rb.integers(0, len(g), len(g))
            if not (0 < dz[i].sum() < len(i)): continue
            Xi = X[i]; bm.append(wls(Xi, y[i], re_weights(Xi, y[i], se_all[i])[0])[1]); bz.append(wls(Xi, yn[i], np.ones(len(i)))[1])
        ci = lambda v: (100 * np.percentile(v, 2.5), 100 * np.percentile(v, 97.5))
        con.append(dict(dataset=ds, disease=case, cell_type=ct, n_case=nc, n_ctrl=nn, model_diff=100 * b, p=p, zero_read_diff=100 * bn, p_zero=pn, umi_diff=bu, p_umi=pu,
                        m_lo=ci(bm)[0], m_hi=ci(bm)[1], z_lo=ci(bz)[0], z_hi=ci(bz)[1]))
q = bh([c["p"] for c in con]); qn = bh([c["p_zero"] for c in con])
for c, a_, b_ in zip(con, q, qn): c["q"] = a_; c["q_zero"] = b_
with open(os.path.join(TAB, "robust_disease_contrasts.tsv"), "w") as f:
    f.write("dataset\tdisease\tcell_type\tn_case\tn_ctrl\tmodel_diff\tm_lo\tm_hi\tp\tq\tzero_read_diff\tz_lo\tz_hi\tp_zero\tq_zero\tumi_diff\tp_umi\n")
    for c in sorted(con, key=lambda c: c["p"]): f.write(f"{c['dataset']}\t{c['disease']}\t{c['cell_type']}\t{c['n_case']}\t{c['n_ctrl']}\t{c['model_diff']:.2f}\t{c['m_lo']:.2f}\t{c['m_hi']:.2f}\t{c['p']:.4g}\t{c['q']:.4g}\t{c['zero_read_diff']:.2f}\t{c['z_lo']:.2f}\t{c['z_hi']:.2f}\t{c['p_zero']:.4g}\t{c['q_zero']:.4g}\t{c['umi_diff']:.0f}\t{c['p_umi']:.4g}\n")

# ---------- 4. sensitivity to the expected-Y threshold ----------
sens = []
KEY = [("ad_mssm", "vascular leptomeningeal cell"), ("ad_mssm", "endothelial cell"), ("ad_mssm", "microglial cell"), ("retina", "endothelial cell of arteriole"), ("retina", "pericyte")]
for thr in (2, 3, 5, 10):
    for ds in sorted({d for d, _ in KEY}):
        for o in age_slopes(ds, thr):
            if (ds, o["cell_type"]) in KEY: sens.append(dict(threshold=thr, **o))
# grouping by sequencing library in place of donor (loyregroup outputs <name>_lib)
for o in age_slopes("retina_lib", 5) if os.path.exists(os.path.join(EXT, "retina_lib.fitamb.tsv")) else []:
    if ("retina", o["cell_type"]) in KEY: sens.append(dict(threshold="library", **o))
with open(os.path.join(TAB, "robust_threshold_sensitivity.tsv"), "w") as f:
    f.write("threshold\tdataset\tcell_type\tdonors\tslope_per_decade\tp\n")
    for s in sens: f.write(f"{s['threshold']}\t{s['dataset']}\t{s['cell_type']}\t{s['donors']}\t{s['slope_per_decade']:.2f}\t{s['p']:.4g}\n")

# ---------- 5. overdispersion against a common-fraction null ----------
od = []
T = []
for ds in ["ad_mssm", "ad_radc", "retina", "aida", "iaging", "onek1k_geo", "blood", "tsapiens", "ckd_sn", "kidney"]:
    fn = os.path.join(EXT, f"{ds}.fitamb.tsv")
    if not os.path.exists(fn): continue
    for (ct, d), x in donor_values(fn, 5, lambda r: r["disease"] == "normal").items():
        if x["age"] is not None and x["age"] >= 55: T.append((ds, ct, x))
by = {}
for ds, ct, x in T: by.setdefault(ct, []).append(x)
rng = np.random.default_rng(1)
for ct, g in by.items():
    if len(g) < 15: continue
    v = np.array([max(0.0, x["v"]) for x in g]); se = np.array([x["se"] for x in g]); w = 1 / se ** 2
    common = max(0.0, (w * v).sum() / w.sum())
    obs = (v > 0.10).mean()
    sims = np.maximum(0, common + rng.standard_normal((5000, len(g))) * se) > 0.10
    null_share = sims.mean(axis=1); p = ((null_share >= obs).sum() + 1) / 5001
    od.append(dict(cell_type=ct, donors=len(g), share_above_10=100 * obs, null_mean_share=100 * null_share.mean(), p=p))
qo = bh([o["p"] for o in od])
for o, qq in zip(od, qo): o["q"] = qq
with open(os.path.join(TAB, "robust_overdispersion.tsv"), "w") as f:
    f.write("cell_type\tdonors\tshare_above_10pct\tnull_share\tp\tq\n")
    for o in sorted(od, key=lambda o: -o["share_above_10"]): f.write(f"{o['cell_type']}\t{o['donors']}\t{o['share_above_10']:.1f}\t{o['null_mean_share']:.2f}\t{o['p']:.4g}\t{o['q']:.4g}\n")
print("age slopes", len(slopes), "| meta", len(meta), "| contrasts", len(con), "| sensitivity", len(sens), "| overdispersion", len(od))
