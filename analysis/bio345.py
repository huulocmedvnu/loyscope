# analyses 3-5 on consolidated table (numbers from Rust fits)
import csv, statistics as st, random, math, glob
import os
RES = os.environ.get("LOY_RESULTS", "results")
EXT = os.path.join(RES, "extract")
TAB = os.path.join(RES, "tables")
os.makedirs(TAB, exist_ok=True)

R = [r for r in csv.DictReader(open(os.path.join(TAB, "loy_table.tsv")), delimiter="\t")]
# --- 4: clonal overdispersion across donors (age >= 55, normal), per cell type with >= 15 donors
print("## 4. overdispersion across donors (age>=55): share of donors with LOY>10%, CV, donor-level dispersion rho")
don = {}
for r in R:
    if r["disease"] != "normal" or r["age"] == "None" or float(r["age"]) < 55: continue
    e = don.setdefault((r["cell_type"], r["dataset"], r["donor"]), [0, 0.0]); n = int(r["n"]); e[0] += n; e[1] += n * max(0.0, float(r["loy"]))
ct = {}
for (c, ds, d), (n, s) in don.items(): ct.setdefault(c, []).append((s / n, n))
rows = []
for c, v in ct.items():
    if len(v) < 15: continue
    p = st.mean(x for x, _ in v); vo = st.pvariance([x for x, _ in v]); vb = st.mean(p * (1 - p) / n for _, n in v)
    rho = (vo - vb) / (p * (1 - p)) if p > 0 else float("nan")
    rows.append((c, len(v), p, sum(1 for x, _ in v if x > 0.10) / len(v), (vo ** 0.5) / p if p > 0 else float("nan"), rho))
for c, n, p, hi, cv, rho in sorted(rows, key=lambda x: -x[5])[:14]: print(f"  {c[:36]:36s} donors {n:4d} mean {p:.3f} >10%: {hi:.2f} CV {cv:5.2f} rho {rho:.3f}")
# --- 3: hematopoietic hierarchy, bone marrow dataset, within donor: progenitors vs mature
print("## 3. hematopoietic hierarchy (bone marrow atlas + TS bone marrow), donor-paired")
PROG = ["hematopoietic stem", "multipotent progenitor", "common myeloid progenitor", "common lymphoid progenitor", "megakaryocyte progenitor", "erythroid progenitor", "granulocyte monocyte progenitor"]
MAT = {"monocyte": ["monocyte"], "erythroid": ["erythroblast"], "B lineage": ["B cell", "pre-B", "pro-B"], "T/NK": ["T cell", "natural killer"]}
bm = [r for r in R if r["dataset"] in ("bone_marrow",) or (r["dataset"] == "tsapiens" and r["tissue"] == "bone marrow")]
dd = {}
for r in bm:
    c = r["cell_type"]; k = "progenitor" if any(w in c for w in PROG) else next((m for m, ws in MAT.items() if any(w in c for w in ws)), None)
    if not k: continue
    e = dd.setdefault((r["dataset"] + r["donor"], k), [0, 0.0]); n = int(r["n"]); e[0] += n; e[1] += n * float(r["loy"])
donors = {d for d, _ in dd}
for m in MAT:
    pairs = [(dd[(d, "progenitor")][1] / dd[(d, "progenitor")][0], dd[(d, m)][1] / dd[(d, m)][0]) for d in donors if (d, "progenitor") in dd and (d, m) in dd]
    if len(pairs) < 4: continue
    diff = [b - a for a, b in pairs]; pos = sum(1 for x in diff if x > 0)
    print(f"  {m:10s} donors {len(pairs):2d} | progenitor mean {st.mean(a for a,_ in pairs):.4f} mature mean {st.mean(b for _,b in pairs):.4f} | mature>progenitor in {pos}/{len(pairs)}")
# --- 5: selection: across cell types within dataset, LOY (age>=55) vs intact Y expression per UMI (mu)
print("## 5. LOY vs intact Y expression rate across cell types (within dataset, Spearman)")
mu = {}
for fn in glob.glob(os.path.join(EXT, "*.fitamb.tsv")):
    ds = fn.split("/")[-1].replace(".fitamb.tsv", "")
    for r in csv.DictReader(open(fn), delimiter="\t"):
        if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 5: mu.setdefault((ds, r["cell_type"]), []).append(float(r["mu"]))
def spear(x, y):
    def rk(v): s = sorted(range(len(v)), key=lambda i: v[i]); o = [0]*len(v); [o.__setitem__(i, j) for j, i in enumerate(s)]; return o
    a, b = rk(x), rk(y); n = len(x); return 1 - 6 * sum((a[i]-b[i])**2 for i in range(n)) / (n*(n*n-1))
for ds in ["ad_mssm", "aida", "iaging", "onek1k_geo", "retina", "tsapiens", "ckd_sn"]:
    pts = [(st.mean(x for x, _ in ct2), st.median(mu[(ds, c)])) for c, ct2 in
           [(c, [(s / n, n) for (cc, dd2, d), (n, s) in don.items() if cc == c and dd2 == ds]) for c in {c for (c, dd2, _) in don if dd2 == ds}] if len(ct2) >= 5 and (ds, c) in mu]
    if len(pts) >= 6: print(f"  {ds:12s} cell types {len(pts):3d} rho(LOY, Y expression) {spear([a for a,_ in pts],[b for _,b in pts]):+.2f}")
