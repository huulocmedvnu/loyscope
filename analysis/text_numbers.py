# Every number quoted in the Results and Discussion, computed from the pipeline outputs with the definitions stated here.
# writes tables/text_numbers.json and a readable tables/text_numbers.txt
# Definitions
#   level of LOY in a cell type: mean over donors of the donor value (inverse-variance pooled groups, control estimate
#     subtracted, negative values set to zero), men aged 55 or older, groups with >= 5 expected Y reads per intact cell
#   slopes, q values, meta-analysis, contrasts and overdispersion: analysis/robust_stats.py
import csv, glob, json, math, os, re, statistics as st

RES = os.environ.get("LOY_RESULTS", "results"); E = os.path.join(RES, "extract"); T = os.path.join(RES, "tables")
rows = lambda f: list(csv.DictReader(open(f), delimiter="\t"))
DEC = {"third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75, "ninth": 85}
def age(s):
    m = re.match(r"(\d+)-year-old", s)
    if m: return float(m.group(1))
    m = re.match(r"(\d+)-(\d+) year-old", s)
    if m: return (float(m.group(1)) + float(m.group(2))) / 2
    m = re.match(r"(\w+) decade", s); return DEC.get(m.group(1)) if m else None
OUT = {}; TXT = []
def put(key, val, note=""):
    OUT[key] = val; TXT.append(f"{key}: {val}{'   # ' + note if note else ''}")
SKIP = {"onek1k", "retina_lib"}
HEART = lambda r: any(k in r["tissue"].lower() for k in ("heart", "ventric", "atri", "myocard"))

def male_groups(fn, thr, sel=lambda r: True):
    return [r for r in rows(fn) if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= thr and sel(r)]
def donor_levels(fn, ct, thr=5, sel=lambda r: True, minage=55):
    d = {}
    for r in male_groups(fn, thr, lambda r: r["cell_type"] == ct and sel(r)):
        a = age(r["stage"])
        if a is not None and a >= minage: d.setdefault(r["donor"], []).append(r)
    vals = []
    for g in d.values():
        w = [1 / max((float(x["f_hi"]) - float(x["f_lo"])) / 3.92, 0.005) ** 2 if x["f_lo"] not in ("NaN", "nan") else 1.0 for x in g]
        v = sum(wi * (float(x["f"]) - float(x["pseudo_f_med"])) for wi, x in zip(w, g)) / sum(w)
        vals.append(max(0.0, v))
    return (round(100 * st.mean(vals), 1), len(vals)) if vals else (None, 0)

# ---- scale of the study ----
ORG = {"blood": ["blood"], "bone marrow": ["bone marrow"], "lymphoid organs": ["lymph node", "spleen", "thymus"], "brain": ["cortex", "prefrontal"],
       "eye": ["eye", "retin", "choroid", "cornea", "sclera", "conjunctiva", "lacrimal", "eyelid", "chorioretinal"],
       "heart": ["heart", "cardiac", "ventric", "septum", "atrium", "myocard"], "blood vessels": ["aorta", "coronary", "vena cava", "vena_cava", "vasculature"],
       "lung and airway": ["lung", "airway", "trachea", "nose"], "kidney": ["kidney", "renal"], "bladder": ["bladder"], "prostate": ["prostate"], "testis": ["testis"],
       "stomach": ["stomach"], "small intestine": ["duodenum", "jejunum", "ileum", "small intestine"], "large intestine": ["colon", "large intestine"], "liver": ["liver"],
       "pancreas": ["pancrea", "islet"], "salivary glands": ["parotid", "sublingual", "submandibular", "saliva"], "tongue and mouth": ["tongue", "buccal"], "skin": ["skin"],
       "adipose tissue": ["adipose", "fat pad"], "skeletal muscle": ["muscle of", "skeletal muscle", "muscle tissue", "rectus"], "mammary gland": ["mammary"],
       "uterus": ["uterus", "endometrium", "myometrium"], "ovary": ["ovary"], "inner ear": ["crista", "utricle"], "adrenal gland": ["adrenal"], "placenta": ["placenta"]}
cells = donors = 0; alld = set(); mc = 0; men = set(); organs = set(); sites = set(); dis = set(); ct5 = set(); g5 = 0; c5 = 0; m5 = set()
for fn in sorted(glob.glob(os.path.join(E, "*.fitamb.tsv"))):
    n = os.path.basename(fn).split(".")[0]
    if n in ("onek1k_geo", "retina_lib"): continue  # OneK1K counted once (CELLxGENE version); library regrouping duplicates retina
    for r in rows(fn):
        k = int(r["n"]); cells += k; alld.add((n, r["donor"]))
        if r["sex"] != "male": continue
        t = r["tissue"].lower(); o = [a for a, v in ORG.items() if any(x in t for x in v)]
        organs.add(o[0] if o else t); sites.add(t); mc += k; men.add((n, r["donor"]))
        dis.update(a for a in r["disease"].split(" || ") if a != "normal")
        if r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 5 and r["assay"].startswith("10x") and n != "gb2024":
            ct5.add(r["cell_type"]); g5 += 1; c5 += k; m5.add((n, r["donor"]))
datasets = [os.path.basename(f).split(".")[0] for f in glob.glob(os.path.join(E, "*.fitamb.tsv"))]
put("datasets", len([d for d in datasets if d not in ("onek1k_geo", "retina_lib", "gb2024")]), "CELLxGENE atlases and the ENCODE collection")
put("mssm_men", len({r["donor"] for r in rows(os.path.join(E, "ad_mssm.fitamb.tsv")) if r["sex"] == "male"}))
put("total_cells_millions", round(cells / 1e6, 1)); put("total_donors", len(alld)); put("male_cells_millions", round(mc / 1e6, 1))
put("men", len(men)); put("organs_male", len(organs)); put("anatomical_sites_male", len(sites)); put("diseases_male", len(dis))
put("cell_types_expY5", len(ct5)); put("groups_expY5", g5); put("men_expY5", len(m5)); put("cells_expY5_millions", round(c5 / 1e6, 1))

# ---- measurement ----
N = json.load(open(os.path.join(T, "numbers.json")))
put("male_groups_all", N["male_groups_all"]); put("pseudo_zero_rule_by_bin", N["pseudo_zero_rule_by_bin"]); put("pseudo_model_by_bin", N["pseudo_model_by_bin"])
th = N["thinning"]; put("thinning_zero_rule_range", [min(v["naive"] for v in th.values()), max(v["naive"] for v in th.values())])
put("thinning_rf_range", [min(v["rf"] for v in th.values()), max(v["rf"] for v in th.values())]); put("thinning_model_max", max(v["zinb"] for v in th.values()))
p95 = {}
for fn in glob.glob(os.path.join(E, "*.fitamb.tsv")):
    n = os.path.basename(fn).split(".")[0]
    if n in SKIP: continue
    v = sorted(float(r["pseudo_f_med"]) for r in male_groups(fn, 5) if r["assay"].startswith("10x") or "GEO" in r["assay"])
    if len(v) >= 20: p95[n] = v[int(0.95 * len(v))]
put("atlases_control_p95_below_1pct", [sum(x < 0.01 for x in p95.values()), len(p95)])
DR = json.load(open(os.path.join(T, "depth_retention.json")))
put("retained_expY5", [round(100 * DR["total"][1] / DR["total"][0]), DR["total"][0]])
put("control_p95_by_band_pct", {k: round(100 * v[2], 1) for k, v in DR["band"].items()})

# ---- loss of the chromosome ----
par = {}
cur = None
for l in open(os.path.join(T, "par_dosage.txt")):
    m = re.match(r"== .*/([a-z_0-9]+)\.cellde\.tsv", l)
    if m: cur = m.group(1); continue
    m = re.match(r"\s+(PAR1|X gametologs).*LOY log2FC median ([-+0-9.]+)", l)
    if m and cur: par.setdefault(cur, {})[m.group(1)] = float(m.group(2))
put("par1_log2", {k: v.get("PAR1") for k, v in par.items() if "PAR1" in v}); put("gametolog_log2", {k: v.get("X gametologs") for k, v in par.items() if "X gametologs" in v})
dn4 = json.load(open(os.path.join(T, "dn4_estimates.json"))); am = json.load(open(os.path.join(E, "gb2024_DN4.atacmix.json")))
put("dn4", {k: round(v, 1) for k, v in dn4.items()})
put("dn4_msy_per_10k", {"intact": round(am["intact_rate_per_fragment"] * 1e4, 1), "female": round(am["background_per_fragment"] * 1e4, 2)})

# ---- where loss occurs ----
for ds, ct in [("ad_mssm", "vascular leptomeningeal cell"), ("ad_mssm", "endothelial cell"), ("ad_mssm", "microglial cell"),
               ("retina", "endothelial cell of arteriole"), ("retina", "endothelial cell of venule"), ("retina", "pericyte"), ("retina", "smooth muscle cell")]:
    put(f"level_55plus[{ds}|{ct}]", donor_levels(os.path.join(E, f"{ds}.fitamb.tsv"), ct))
neur = []
for ds in ("ad_mssm", "ad_radc"):
    for ct in {r["cell_type"] for r in rows(os.path.join(E, f"{ds}.fitamb.tsv"))}:
        if re.search(r"neuron|^oligodendrocyte$|^astrocyte$", ct):
            lv, n = donor_levels(os.path.join(E, f"{ds}.fitamb.tsv"), ct)
            if n >= 10: neur.append(lv)
put("neural_level_range_55plus", [min(neur), max(neur)])
S = rows(os.path.join(T, "robust_age_slopes.tsv"))
def slope(ds, ct):
    r = [x for x in S if x["dataset"] == ds and x["cell_type"] == ct]
    return {k: (round(float(r[0][k]), 4) if k != "donors" else int(r[0][k])) for k in ("slope_per_decade", "p", "q", "donors")} if r else None
for ds, ct in [("ad_mssm", "vascular leptomeningeal cell"), ("ad_mssm", "endothelial cell"), ("ad_mssm", "microglial cell"), ("ad_mssm", "L2/3-6 intratelencephalic projecting glutamatergic neuron"),
               ("retina", "endothelial cell of arteriole"), ("retina", "endothelial cell of venule"), ("retina", "pericyte"), ("encode_heart", "cardiomyocyte")]:
    put(f"slope[{ds}|{ct}]", slope(ds, ct))
radc = [x for x in S if x["dataset"] == "ad_radc" and not re.search("neuron", x["cell_type"])]
put("radc_nonneuronal_min_p_and_donor_range", [min(float(x["p"]) for x in radc), min(int(x["donors"]) for x in radc), max(int(x["donors"]) for x in radc)])
put("sensitivity_rows", [(r["threshold"], r["dataset"], r["cell_type"], round(float(r["slope_per_decade"]), 2), round(float(r["p"]), 4)) for r in rows(os.path.join(T, "robust_threshold_sensitivity.tsv"))])
# heart (ENCODE, cardiac samples, >= 2 expected reads)
H = male_groups(os.path.join(E, "encode_heart.fitamb.tsv"), 2, lambda r: HEART(r) and r["cell_type"] == "cardiomyocyte")
d = {}
for r in H: d.setdefault(r["donor"], []).append(r)
pts = [((age(g[0]["stage"]) or 0), sum((float(x["f"]) - float(x["pseudo_f_med"])) * int(x["n"]) for x in g) / sum(int(x["n"]) for x in g), sum(int(x["n"]) for x in g),
        sum(float(x["naive_zero"]) * int(x["n"]) for x in g) / sum(int(x["n"]) for x in g)) for g in d.values()]
old = [max(0.0, p[1]) for p in pts if p[0] >= 55]
put("heart_cardiomyocytes", {"donors": len(pts), "nuclei": sum(p[2] for p in pts), "age_range": [min(p[0] for p in pts if p[0] > 0), max(p[0] for p in pts)],
                             "mean_55plus_pct": round(100 * st.mean(old), 2), "max_55plus_pct": round(100 * max(old), 1), "n_55plus": len(old),
                             "zero_read_pct": round(100 * sum(p[3] * p[2] for p in pts) / sum(p[2] for p in pts), 1)})
# blood
MC = rows(os.path.join(T, "robust_meta_cohorts.tsv")); M = rows(os.path.join(T, "robust_meta.tsv"))
put("meta", {r["cell_class"]: [round(float(r["slope_per_decade"]), 2), round(float(r["ci_lo"]), 2), round(float(r["ci_hi"]), 2), int(float(r["I2"]))] for r in M})
put("cohort_slopes", {r["cell_class"] + "|" + r["cohort"]: round(float(r["slope_per_decade"]), 2) for r in MC})
put("blood_celltype_slopes_q", {f"{x['dataset']}|{x['cell_type']}": [round(float(x["slope_per_decade"]), 2), round(float(x["q"]), 4)] for x in S if x["dataset"] in ("onek1k_geo", "aida", "iaging", "blood") and re.search("monocyte|natural killer|regulatory", x["cell_type"])})
put("level_all_ages[iaging|memory regulatory T cell]", donor_levels(os.path.join(E, "iaging.fitamb.tsv"), "memory regulatory T cell", thr=1.5, minage=0))

# ---- systemic or local ----
put("within_donor", {l.split()[0]: float(re.search(r"variance ([0-9.]+) P ([0-9.]+)", l).group(1)) for l in open(os.path.join(T, "within_donor.txt")) if "variance" in l})
put("within_donor_p", {l.split()[0]: float(re.search(r"variance ([0-9.]+) P ([0-9.]+)", l).group(2)) for l in open(os.path.join(T, "within_donor.txt")) if "variance" in l})
m = re.search(r"rho\(LOY, log lifespan\) = ([-0-9.]+), n = (\d+) classes, permutation P = ([0-9.]+)", open(os.path.join(T, "turnover.txt")).read())
put("lifespan", {"rho": float(m.group(1)), "classes": int(m.group(2)), "p": float(m.group(3))})
O = rows(os.path.join(T, "robust_overdispersion.tsv"))
put("overdispersion", {r["cell_type"]: [float(r["share_above_10pct"]), float(r["null_share"]), float(r["p"])] for r in O[:6]})

# ---- cells without Y ----
dd = rows(os.path.join(E, "iaging.de_depth.tsv")); dc = rows(os.path.join(E, "iaging.de_depth_control.tsv"))
Y = {"RPS4Y1", "ZFY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "TXLNGY", "NLGN4Y", "PRKY", "TMSB4Y"}
mt = [r for r in dd if r["cell_type"] == "memory regulatory T cell"]
put("de_depth_memory_treg_q01_nonY", sum(1 for r in mt if float(r["q_bh"]) < 0.1 and r["gene"] not in Y))
put("de_depth_control_memory_treg_q01", sum(1 for r in dc if r["cell_type"] == "memory regulatory T cell" and float(r["q_bh"]) < 0.1))
put("de_depth_other_types_q01", {ct: sum(1 for r in dd if r["cell_type"] == ct and float(r["q_bh"]) < 0.1) for ct in sorted({r["cell_type"] for r in dd}) if ct != "memory regulatory T cell"})
put("de_depth_genes", {r["gene"]: [round(float(r["median_log2fc"]), 2), float(r["q_bh"])] for r in mt if r["gene"] in ("PRDM1", "RTKN2", "HELLS", "S1PR4", "SATB1", "IL7R", "ZFP36L1", "ZFP36L2", "CD99", "SLC25A6", "LTBP3")})

# ---- published estimates and disease ----
put("published_comparison", N["published_comparison"]); put("onek1k_cd14_donor_median_model", N["onek1k_cd14_donor_median_model"]); put("onek1k_cd14_donors", N["onek1k_cd14_donors"])
D = rows(os.path.join(T, "robust_disease_contrasts.tsv"))
put("disease_contrasts", {"n": len(D), "model_p05": sum(float(r["p"]) < 0.05 for r in D), "zero_p05": sum(float(r["p_zero"]) < 0.05 for r in D),
                          "model_q01": sum(float(r["q"]) < 0.1 for r in D), "zero_q01": sum(float(r["q_zero"]) < 0.1 for r in D)})
pick = lambda ds, dis, ct: next(({k: (round(float(r[k]), 2) if k in ("model_diff", "zero_read_diff") else float(r[k])) for k in ("model_diff", "p", "q", "zero_read_diff", "p_zero", "n_case", "n_ctrl")}
                                 for r in D if r["dataset"] == ds and r["disease"] == dis and r["cell_type"] == ct), None)
for key, args in {"kidney_TAL_CKD": ("kidney", "chronic kidney disease", "kidney loop of Henle thick ascending limb epithelial cell"),
                  "nuclei_PT_AKI": ("ckd_sn", "acute kidney injury", "epithelial cell of proximal tubule"),
                  "dcm_unannotated": ("dcm", "dilated cardiomyopathy", "unknown"), "dcm_cardiomyocyte": ("dcm", "dilated cardiomyopathy", "cardiac muscle cell"),
                  "microglia_fresh": ("ad_microglia", "Alzheimer disease", "microglial cell"), "microglia_mssm": ("ad_mssm", "Alzheimer disease", "microglial cell"),
                  "microglia_radc": ("ad_radc", "Alzheimer disease", "microglial cell"), "opc_mssm": ("ad_mssm", "Alzheimer disease", "oligodendrocyte precursor cell")}.items():
    put(f"contrast[{key}]", pick(*args))
best = min(D, key=lambda r: float(r["p"])); put("contrast_smallest_p", [best["dataset"], best["cell_type"], round(float(best["model_diff"]), 2), float(best["p"]), float(best["q"])])

# ---- ambient sensitivity: contamination fraction x3 and x10, and the relative parameterisation ----
def compare(suffix):
    d = []; pd = []
    for fn in glob.glob(os.path.join(E, "*.fitamb.tsv")):
        n = os.path.basename(fn).split(".")[0]; alt = os.path.join(E, f"{n}.{suffix}.tsv")
        if n in SKIP or not os.path.exists(alt): continue
        a = {r["gid"]: r for r in male_groups(fn, 5)}; b = {r["gid"]: r for r in rows(alt) if r["k"] != "NaN"}
        for g, r in a.items():
            if g in b:
                d.append(100 * ((float(b[g]["f"]) - float(b[g]["pseudo_f_med"])) - (float(r["f"]) - float(r["pseudo_f_med"]))))
                pd.append(100 * (float(b[g]["pseudo_f_med"]) - float(r["pseudo_f_med"])))
    if not d: return None
    d.sort(); pd.sort(); ad = sorted(abs(x) for x in d)
    return {"groups": len(d), "loy_change_mean": round(st.mean(d), 2), "loy_change_p5": round(d[int(0.05 * len(d))], 2), "loy_change_p95": round(d[int(0.95 * len(d))], 2),
            "loy_abs_p95": round(ad[int(0.95 * len(ad))], 2), "control_change_mean": round(st.mean(pd), 2), "control_change_p95": round(pd[int(0.95 * len(pd))], 2)}
for s in ("fitamb_x3", "fitamb_x10", "fitrel"): put(f"ambient_sensitivity[{s}]", compare(s))
json.dump(OUT, open(os.path.join(T, "text_numbers.json"), "w"), indent=1)
open(os.path.join(T, "text_numbers.txt"), "w").write("\n".join(TXT) + "\n")
print("\n".join(TXT))
