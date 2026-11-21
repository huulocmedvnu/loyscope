# share of male 10x groups retained by expected Y count, and control estimates by expected-count band
# writes tables/depth_retention.json and tables/Supplementary_Table_3_depth_retention.tsv
import csv, glob, os, json, statistics as st
res = os.environ.get("LOY_RESULTS", "results"); E = os.path.join(res, "extract"); T = os.path.join(res, "tables")
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from datasets import EXCLUDE as SKIP
tot = [0, 0, 0]; per = {}; band = {"1.5-2": [], "2-3": [], "3-5": [], ">=5": []}
for fn in sorted(glob.glob(os.path.join(E, "*.fitamb.tsv"))):
    ds = os.path.basename(fn).split(".")[0]
    if ds in SKIP: continue
    R = [r for r in csv.DictReader(open(fn), delimiter="\t") if r["sex"] == "male" and r["k"] != "NaN" and (r["assay"].startswith("10x") or "GEO" in r["assay"])]
    if not R: continue
    e = [float(r["mu"]) * float(r["mean_umi"]) for r in R]
    a = sum(x >= 5 for x in e); b = sum(x >= 1.5 for x in e)
    per[ds] = (len(R), round(100 * a / len(R)), round(100 * b / len(R)))
    tot[0] += len(R); tot[1] += a; tot[2] += b
    for r, x in zip(R, e):
        k = ">=5" if x >= 5 else "3-5" if x >= 3 else "2-3" if x >= 2 else "1.5-2" if x >= 1.5 else None
        if k: band[k].append(float(r["pseudo_f_med"]))
for v in band.values(): v.sort()
out = dict(total=tot, per=per, band={k: (len(v), st.median(v), v[int(0.95 * len(v))]) for k, v in band.items() if v})
os.makedirs(T, exist_ok=True)
json.dump(out, open(os.path.join(T, "depth_retention.json"), "w"), indent=1)
with open(os.path.join(T, "Supplementary_Table_3_depth_retention.tsv"), "w") as w:
    w.write("dataset\tmale_groups\tretained_expY_ge5_pct\tretained_expY_ge1.5_pct\n")
    for k, (n, a, b) in sorted(per.items()): w.write(f"{k}\t{n}\t{a}\t{b}\n")
    w.write("\nexpected_Y_band\tgroups\tcontrol_median\tcontrol_p95\n")
    for k, (n, m, p) in out["band"].items(): w.write(f"{k}\t{n}\t{m:.4f}\t{p:.4f}\n")
print("male groups", tot[0], "retained >=5: %.0f%%" % (100 * tot[1] / tot[0]), ">=1.5: %.0f%%" % (100 * tot[2] / tot[0]))
