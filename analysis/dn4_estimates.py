# kidney donor DN4 (published LOY study): fraction of cells without Y from chromatin (atacmix), the detection model,
# the zero-read rule and the density trough (tenx) in the same donor; writes tables/dn4_estimates.json
import csv, json, os, sys
res = os.environ.get("LOY_RESULTS", "results"); E = os.path.join(res, "extract")
atac = json.load(open(os.path.join(E, "gb2024_DN4.atacmix.json")))
R = [r for r in csv.DictReader(open(os.path.join(E, "gb2024.fitamb.tsv")), delimiter="\t")
     if r["donor"] == "GSE195460_DN4" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 5]
n = sum(int(r["n"]) for r in R)
model = 100 * sum(int(r["n"]) * (float(r["f"]) - float(r["pseudo_f_med"])) for r in R) / n
naive = 100 * sum(int(r["n"]) * float(r["naive_zero"]) for r in R) / n
tr = [l.rstrip("\n").split("\t") for l in open(os.path.join(E, "gb2024.trough.tsv")) if l.startswith("GSE195460_DN4\tALL")][0]
d = dict(atac=100 * atac["loy_fraction"], model=model, naive=naive, trough=100 * float(tr[3]))
os.makedirs(os.path.join(res, "tables"), exist_ok=True)
json.dump(d, open(os.path.join(res, "tables", "dn4_estimates.json"), "w"), indent=1)
print(json.dumps(d))
