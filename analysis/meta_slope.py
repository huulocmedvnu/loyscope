# random-effects (DerSimonian-Laird) meta-analysis of donor-level age slopes of corrected LOY, harmonized blood classes
import csv, re, math, statistics as st, numpy as np, sys, os
EXT = os.path.join(os.environ.get("LOY_RESULTS", "results"), "extract")
DEC = {"third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75, "ninth": 85}
def age(s):
    m = re.match(r"(\d+)-year-old", s)
    if m: return float(m.group(1))
    m = re.match(r"(\w+) decade", s); return DEC.get(m.group(1)) if m else None
CLS = [("classical monocyte", lambda c: ("CD14-positive" in c and "CD16-positive" not in c and "monocyte" in c) or c == "classical monocyte"),
       ("non-classical monocyte", lambda c: "CD16-positive" in c and "monocyte" in c),
       ("NK cell", lambda c: "natural killer" in c and "T cell" not in c),
       ("regulatory T cell", lambda c: "regulatory T" in c),
       ("naive B cell", lambda c: c == "naive B cell"),
       ("memory B cell", lambda c: "memory B cell" in c),
       ("CD4 T cell", lambda c: "CD4-positive" in c and "regulatory" not in c),
       ("CD8 T cell", lambda c: "CD8-positive" in c and "regulatory" not in c)]
cohorts = [("OneK1K (GEO)", os.path.join(EXT, "onek1k_geo.fit.tsv")), ("AIDA", os.path.join(EXT, "aida.fit.tsv")), ("Immune Aging", os.path.join(EXT, "iaging.fit.tsv")), ("Indonesia blood", os.path.join(EXT, "blood.fit.tsv"))]
thr = 1.5
data = {}
for name, fn in cohorts:
    try: R = [r for r in csv.DictReader(open(fn), delimiter="\t") if r["sex"] == "male" and r["disease"] == "normal" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= thr]
    except FileNotFoundError: continue
    for cl, fnc in CLS:
        dd = {}
        for r in R:
            if fnc(r["cell_type"]): dd.setdefault(r["donor"], []).append(r)
        pts = []
        for v in dd.values():
            a = age(v[0]["stage"]); n = sum(int(x["n"]) for x in v)
            if a is None: continue
            pts.append((a, sum((float(x["f"]) - float(x["pseudo_f_med"])) * int(x["n"]) for x in v) / n))
        if len(pts) < 10: continue
        x = np.array([p[0] for p in pts]) / 10; y = np.array([p[1] for p in pts])
        X = np.c_[np.ones_like(x), x]; b, res, *_ = np.linalg.lstsq(X, y, rcond=None)
        s2 = ((y - X @ b) ** 2).sum() / (len(y) - 2); se = math.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
        data.setdefault(cl, []).append((name, b[1], se, len(pts)))
for cl, rows in data.items():
    if len(rows) < 2: continue
    w = [1 / s ** 2 for _, _, s, _ in rows]; bf = sum(wi * b for wi, (_, b, _, _) in zip(w, rows)) / sum(w)
    Q = sum(wi * (b - bf) ** 2 for wi, (_, b, _, _) in zip(w, rows)); df = len(rows) - 1
    C = sum(w) - sum(wi ** 2 for wi in w) / sum(w); tau2 = max(0, (Q - df) / C)
    wr = [1 / (s ** 2 + tau2) for _, _, s, _ in rows]; br = sum(wi * b for wi, (_, b, _, _) in zip(wr, rows)) / sum(wr); ser = math.sqrt(1 / sum(wr))
    I2 = max(0, (Q - df) / Q) * 100 if Q > 0 else 0
    z = br / ser; p = math.erfc(abs(z) / math.sqrt(2))
    print(f"{cl:24s} pooled slope/10y {br*100:+.2f}% [{(br-1.96*ser)*100:+.2f}, {(br+1.96*ser)*100:+.2f}] P {p:.1e} I2 {I2:.0f}% | " + "; ".join(f"{n} {b*100:+.2f}±{s*100:.2f} (n{k})" for n, b, s, k in rows))
