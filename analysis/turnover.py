# analysis 1: LOY by cell class vs cell lifespan (Sender & Milo 2021); donors aged >= 55, normal tissue
import csv, statistics as st, math, random
import os
RES = os.environ.get("LOY_RESULTS", "results")
EXT = os.path.join(RES, "extract")
TAB = os.path.join(RES, "tables")
os.makedirs(TAB, exist_ok=True)

life = {r["cell_type"]: (float(r["lifespan_days"]) if r["lifespan_days"] != "None" else float(r["lifespan_rodent_days"]) if r["lifespan_rodent_days"] != "None" else None) for r in csv.DictReader(open(os.path.join(os.path.dirname(__file__), "..", "resources", "lifespan.tsv")), delimiter="\t")}
MAP = [("Neurons", ["neuron", "interneuron", "projecting", "glutamatergic", "corticothalamic", "intratelencephalic", "Purkinje", "granule cell"]),
       ("Glial cells", ["astrocyte", "oligodendrocyte"]), ("Endothelial cells", ["endothelial"]), ("Monocytes", ["monocyte"]),
       ("Neutrophils", ["neutrophil"]), ("Mature B cells", ["B cell"]), ("Mature T cells", ["T cell"]), ("Cardiomyocytes", ["cardiac muscle cell", "cardiomyocyte"]),
       ("Myocytes", ["skeletal muscle fiber", "fast muscle", "slow muscle"]), ("Adipocytes", ["adipocyte"]), ("Hepatocytes", ["hepatocyte"]),
       ("Dermal fibroblasts", ["fibroblast"]), ("Alveolar macrophages", ["alveolar macrophage"]), ("Alveolar epithelial cells", ["pneumocyte", "alveolar type"]),
       ("Bronchial epithelial cells", ["club cell", "ciliated", "respiratory basal"]), ("Epidermal cells", ["keratinocyte"]),
       ("Small intestine epithelia", ["enterocyte"]), ("Kupffer cells", ["Kupffer"]), ("Stellate cells", ["hepatic stellate"])]
def cls(c):
    for k, kws in MAP:
        if any(w.lower() in c.lower() for w in kws): return k
    return None
R = [r for r in csv.DictReader(open(os.path.join(TAB, "loy_table.tsv")), delimiter="\t") if r["disease"] == "normal" and r["age"] != "None" and float(r["age"]) >= 55]
don = {}
for r in R:
    k = cls(r["cell_type"])
    if not k: continue
    e = don.setdefault((k, r["dataset"], r["donor"]), [0, 0.0]); n = int(r["n"]); e[0] += n; e[1] += n * float(r["loy"])
by = {}
for (k, ds, d), (n, s) in don.items(): by.setdefault(k, []).append((s / n, ds))
pts = []
print(f"{'class':28s} donors datasets  median LOY   mean   lifespan_days")
for k, v in sorted(by.items(), key=lambda x: -st.mean(a for a, _ in x[1])):
    if len(v) < 5: continue
    L = life.get(k); pts.append((k, st.mean(a for a, _ in v), L))
    print(f"{k:28s} {len(v):6d} {len({ds for _, ds in v}):8d}  {st.median(a for a,_ in v):9.3f} {st.mean(a for a,_ in v):7.3f}   {L}")
p = [(m, math.log10(L)) for _, m, L in pts if L]
def rank(x): s = sorted(range(len(x)), key=lambda i: x[i]); r = [0]*len(x); [r.__setitem__(i, j) for j, i in enumerate(s)]; return r
a, b = rank([x[0] for x in p]), rank([x[1] for x in p]); n = len(p)
rho = 1 - 6 * sum((a[i] - b[i]) ** 2 for i in range(n)) / (n * (n * n - 1))
random.seed(1); ge = 0; bb = b[:]
for _ in range(10000):
    random.shuffle(bb); r2 = 1 - 6 * sum((a[i] - bb[i]) ** 2 for i in range(n)) / (n * (n * n - 1)); ge += abs(r2) >= abs(rho)
print(f"Spearman rho(LOY, log lifespan) = {rho:+.2f}, n = {n} classes, permutation P = {(ge+1)/10001:.3f}")
