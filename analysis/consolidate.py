# consolidate corrected LOY (fitamb, 10x, expY>=5, male, normal) into one donor x cell type table with age
import csv, glob, re, os
import os
RES = os.environ.get("LOY_RESULTS", "results")
EXT = os.path.join(RES, "extract")
TAB = os.path.join(RES, "tables")
os.makedirs(TAB, exist_ok=True)

DEC = {"third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75, "ninth": 85}
def age(s):
    m = re.match(r"(\d+)-year-old", s)
    if m: return float(m.group(1))
    m = re.match(r"(\d+)-(\d+) year-old", s)
    if m: return (float(m.group(1)) + float(m.group(2))) / 2
    m = re.match(r"(\w+) decade", s); return DEC.get(m.group(1)) if m else None
skip = {"onek1k"}  # onek1k_geo is the same cohort
rows = []
for fn in sorted(glob.glob(os.path.join(EXT, "*.fitamb.tsv"))):
    ds = os.path.basename(fn).replace(".fitamb.tsv", "")
    if ds in skip: continue
    for r in csv.DictReader(open(fn), delimiter="\t"):
        if r["sex"] != "male" or r["k"] == "NaN" or not r["assay"].startswith("10x") and "GEO" not in r["assay"]: continue
        if float(r["mu"]) * float(r["mean_umi"]) < 5: continue
        rows.append(dict(dataset=ds, tissue=r["tissue"], cell_type=r["cell_type"], disease=r["disease"], donor=r["donor"], age=age(r["stage"]), n=int(r["n"]),
                         loy=float(r["f"]) - float(r["pseudo_f_med"]), suspension=r["suspension"], mean_umi=float(r["mean_umi"]), expY=float(r["mu"]) * float(r["mean_umi"])))
with open(os.path.join(TAB, "loy_table.tsv"), "w") as f:
    k = list(rows[0].keys()); f.write("\t".join(k) + "\n")
    for r in rows: f.write("\t".join(str(r[x]) for x in k) + "\n")
print("rows", len(rows), "datasets", len({r['dataset'] for r in rows}), "donors", len({(r['dataset'], r['donor']) for r in rows}), "cell types", len({r['cell_type'] for r in rows}))
