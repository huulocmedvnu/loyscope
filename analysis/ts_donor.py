# analysis 2: Tabula Sapiens, same donor across tissues: does donor explain LOY variance (systemic clone) for immune vs endothelial/stromal?
import csv, statistics as st, random
import os
RES = os.environ.get("LOY_RESULTS", "results")
EXT = os.path.join(RES, "extract")
TAB = os.path.join(RES, "tables")
os.makedirs(TAB, exist_ok=True)

R = [r for r in csv.DictReader(open(os.path.join(TAB, "loy_table.tsv")), delimiter="\t") if r["dataset"] == "tsapiens"]
G = {"myeloid": ["monocyte", "macrophage", "dendritic", "neutrophil"], "lymphoid": ["T cell", "B cell", "natural killer", "plasma"],
     "endothelial": ["endothelial"], "stromal": ["fibroblast", "pericyte", "smooth muscle"]}
def grp(c):
    for k, w in G.items():
        if any(x in c for x in w): return k
for g in G:
    cell = {}
    for r in R:
        if grp(r["cell_type"]) != g: continue
        e = cell.setdefault((r["donor"], r["tissue"]), [0, 0.0]); n = int(r["n"]); e[0] += n; e[1] += n * float(r["loy"])
    pts = [(d, t, s / n) for (d, t), (n, s) in cell.items()]
    donors = {d for d, _, _ in pts}
    multi = [d for d in donors if sum(1 for x in pts if x[0] == d) >= 3]
    pts = [p for p in pts if p[0] in multi]
    if len(multi) < 3: print(g, "too few donors"); continue
    def r2(pp):
        m = st.mean(x[2] for x in pp); tot = sum((x[2] - m) ** 2 for x in pp)
        bd = sum(sum(1 for x in pp if x[0] == d) * (st.mean(x[2] for x in pp if x[0] == d) - m) ** 2 for d in {x[0] for x in pp})
        return bd / tot if tot > 0 else 0
    obs = r2(pts); random.seed(1); lab = [p[0] for p in pts]; ge = 0
    for _ in range(5000):
        random.shuffle(lab); ge += r2([(lab[i], pts[i][1], pts[i][2]) for i in range(len(pts))]) >= obs
    dm = {d: st.mean(x[2] for x in pts if x[0] == d) for d in multi}
    print(f"{g:12s} donors {len(multi)} donor-tissue points {len(pts)} | donor share of variance {obs:.2f} P {(ge+1)/5001:.3f} | donor means " + ", ".join(f"{d}:{v:.3f}" for d, v in sorted(dm.items(), key=lambda x: -x[1])))
