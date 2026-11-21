# collect every number quoted in the manuscript from tool outputs (single source of truth)
import json, math, random, numpy as np
from common import *
N = {}
# pseudo-Y zero-read rule vs model by expected Y bins
bins = [(0, 1), (1, 2), (2, 5), (5, 1e9)]; pn = {b: [] for b in bins}; pm = {b: [] for b in bins}; ng = 0; datasets = set()
for fn in [ext(f"{k}.fitamb.tsv") for k in canonical()]:
    for r in rows(fn):
        if r["sex"] != "male" or r["k"] == "NaN" or not (r["assay"].startswith("10x") or "GEO" in r["assay"]): continue
        e = float(r["mu"]) * float(r["mean_umi"]); ng += 1; datasets.add(os.path.basename(fn))
        for b in bins:
            if b[0] <= e < b[1]: pn[b].append(float(r["pseudo_naive_med"])); pm[b].append(float(r["pseudo_f_med"]))
N["pseudo_zero_rule_by_bin"] = {f"{b[0]}-{b[1]}": round(100 * st.median(pn[b]), 1) for b in bins}
N["pseudo_model_by_bin"] = {f"{b[0]}-{b[1]}": round(100 * st.median(pm[b]), 2) for b in bins}
N["male_groups_all"] = ng
# thinning per atlas
th = {}
for a in ["aida", "retina", "covid", "hlca", "tsapiens"]:
    by = {}
    for r in rows(ext(f"{a}.bench.tsv")): by.setdefault(r["gid"], {})[r["scenario"]] = r
    th[a] = {m: round(st.median(100 * (float(g["thin0.5"][m]) - float(g["observed"][m])) for g in by.values() if "thin0.5" in g and float(g["observed"]["expY"]) >= 2), 2) for m in ("naive", "rf", "zinb")}
    th[a]["groups"] = len(by)
N["thinning"] = th
N["rf_accuracy"] = {a: open(ext(f"{a}.rf_info.txt")).read().split("test_accuracy ")[1].split()[0] for a in ["aida", "retina"] if os.path.exists(ext(f"{a}.rf_info.txt"))}
json.dump(N, open(os.path.join(TAB, "numbers.json"), "w"), indent=1); print(json.dumps(N, indent=1))
# published-claim comparison (cell-weighted, expected Y >= 1.5) and disease contrasts (same code as Figure 4)
def wmean(fn, sel, key):
    R = [r for r in rows(fn) if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 1.5 and sel(r)]
    n = sum(int(r["n"]) for r in R)
    return round(100 * sum(int(r["n"]) * (max(0.0, float(r["f"]) - float(r["pseudo_f_med"])) if key == "model" else float(r["naive_zero"])) for r in R) / n, 1)
cl = {"onek1k_cd14": (ext("onek1k_geo.fitamb.tsv"), lambda r: r["cell_type"] == "CD14-positive monocyte"),
      "kidney_control1": (ext("gb2024.fitamb.tsv"), lambda r: r["donor"] == "GSE151302_Control1"),
      "microglia_ad": (ext("ad_microglia.fitamb.tsv"), lambda r: r["cell_type"] == "microglial cell" and r["disease"] == "Alzheimer disease"),
      "microglia_ctrl": (ext("ad_microglia.fitamb.tsv"), lambda r: r["cell_type"] == "microglial cell" and r["disease"] == "normal")}
N["published_comparison"] = {k: {"zero_read": wmean(f, s, "naive"), "model": wmean(f, s, "model")} for k, (f, s) in cl.items()}
# donor-median OneK1K CD14 monocytes
R = [r for r in rows(ext("onek1k_geo.fitamb.tsv")) if r["sex"] == "male" and r["k"] != "NaN" and r["cell_type"] == "CD14-positive monocyte" and float(r["mu"]) * float(r["mean_umi"]) >= 1.5]
N["onek1k_cd14_donor_median_model"] = round(100 * st.median(max(0.0, float(r["f"]) - float(r["pseudo_f_med"])) for r in R), 1)
N["onek1k_cd14_donors"] = len({r["donor"] for r in R})
N["age_slope_meta"] = rows(os.path.join(TAB, "robust_meta.tsv"))
N["dn4"] = json.load(open(os.path.join(TAB, "dn4_estimates.json")))
json.dump(N, open(os.path.join(TAB, "numbers.json"), "w"), indent=1)
print(json.dumps({k: N[k] for k in ["published_comparison", "onek1k_cd14_donor_median_model", "onek1k_cd14_donors", "dn4"]}, indent=1))
