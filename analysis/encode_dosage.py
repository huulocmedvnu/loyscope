# Chromatin LOY by dosage in ENCODE heart multiome libraries, compared with the RNA estimate of the same nuclei.
# Per-nucleus mixtures are not identifiable in these data (about 3 MSY fragments per intact nucleus, female background
# 29% of the intact rate), so loss is estimated from group-level dosage, as in DNA array studies:
#   q = (MSY fragments - b * total fragments) / mean control-bin fragments,  f = 1 - q / q_ref(cell type)
# with b the female MSY rate per fragment and q_ref pooled over donors aged 35 or younger (including embryonic).
# Negative control: each of the 20 autosomal control bin sets treated as Y against the mean of the others.
# Inputs under <results>: encode/<acc>.multi.groups_cmp.tsv, .multi.nuclei.tsv.gz, .atac.tsv.gz; manifest from $MANIFEST
import csv, gzip, glob, os, sys, json, math
import numpy as np

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
NP = 20; NB = 300; MIN_N = 100
# female background: MSY fragments per fragment in the female libraries (encmulti --female outputs)
_fem = glob.glob(os.path.join(ROOT, "encode", "fem_*.tsv"))
if _fem:
    _t = _y = 0.0
    for fn in _fem:
        for r in csv.DictReader(open(fn), delimiter="\t"):
            _t += float(r["total"]); _y += float(r["chrY_msy"])
    BG = _y / _t
else:
    BG = float(os.environ["ATAC_BG"])
rng = np.random.default_rng(11)
man = {r[0]: r for r in csv.reader(open(os.environ.get("MANIFEST", "config/encode_manifest.tsv")), delimiter="\t")}
def age_of(stage):
    s = stage.split("-")[0]
    return float(s) if s.isdigit() else 0.0  # embryonic -> 0
rows = []
for gfile in sorted(glob.glob(os.path.join(ROOT, "encode/*.multi.groups_cmp.tsv"))):
    acc = os.path.basename(gfile).split(".")[0]
    G = {r["gid"]: r for r in csv.DictReader(open(gfile), delimiter="\t")}
    cg = {}
    with gzip.open(os.path.join(ROOT, f"encode/{acc}.multi.nuclei.tsv.gz"), "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"): cg[r["cell_id"]] = r["gid"]
    per = {}
    with gzip.open(os.path.join(ROOT, f"encode/{acc}.atac.tsv.gz"), "rt") as f:
        rd = csv.reader(f, delimiter="\t"); next(rd)
        for r in rd:
            g = cg.get(r[0])
            if g in G: per.setdefault(g, []).append([float(x) for x in r[1:]])
    for g, m in per.items():
        if len(m) < MIN_N: continue
        M = np.array(m)  # columns: total, Y, X, pseudo1..20
        rg = G[g]; donor = man.get(acc, [acc, acc])[1]
        rows.append(dict(acc=acc, donor=donor, cell_type=rg["cell_type"], age=age_of(rg["stage"]), n=len(m), M=M,
                         f_rna=float(rg["f_rna_model"]), naive=float(rg["naive_zero"]), expY=float(rg["expY_rna"])))
def q_of(S, col):  # S: summed counts vector; col: index of target; control = mean of pseudo bins except target
    pc = [3 + j for j in range(NP) if 3 + j != col]
    ctrl = S[pc].mean()
    sub = S[col] - (BG * S[0] if col == 1 else 0.0)
    return sub / ctrl
for r in rows:
    S = r["M"].sum(0); r["qY"] = q_of(S, 1); r["qP"] = [q_of(S, 3 + j) for j in range(NP)]
    # nucleus bootstrap of qY
    idx = rng.integers(0, r["n"], (NB, r["n"]))
    r["qY_bs"] = np.array([q_of(r["M"][i].sum(0), 1) for i in idx])
    r["qP_bs0"] = np.array([q_of(r["M"][i].sum(0), 3) for i in idx])
cts = sorted({r["cell_type"] for r in rows})
ref = {}
for ct in cts:
    R = [r for r in rows if r["cell_type"] == ct and r["age"] <= 35]
    if len(R) < 2: continue
    S = sum(r["M"].sum(0) for r in R)
    ref[ct] = (q_of(S, 1), [q_of(S, 3 + j) for j in range(NP)])
out = []
for r in rows:
    if r["cell_type"] not in ref: continue
    qr, qpr = ref[r["cell_type"]]
    f = 1 - r["qY"] / qr; se = np.std(1 - r["qY_bs"] / qr)
    fp = [1 - r["qP"][j] / qpr[j] for j in range(NP)]
    out.append(dict(acc=r["acc"], donor=r["donor"], cell_type=r["cell_type"], age=r["age"], n=r["n"], f_dosage=f, se=se,
                    f_control_med=float(np.median(fp)), f_control_sd=float(np.std(fp)), f_rna=r["f_rna"], naive=r["naive"], expY=r["expY"]))
os.makedirs(os.path.join(ROOT, "tables"), exist_ok=True)
with open(os.path.join(ROOT, "tables/encode_dosage.tsv"), "w") as f:
    k = list(out[0].keys()); f.write("\t".join(k) + "\n")
    for o in out: f.write("\t".join(f"{o[x]:.5g}" if isinstance(o[x], float) else str(o[x]) for x in k) + "\n")
# summaries
def spearman(a, b):
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b)); return float(np.corrcoef(ra, rb)[0, 1])
def perm_p(a, b, n=5000):
    r0 = abs(spearman(a, b)); a = np.array(a); c = 0
    for _ in range(n): c += abs(spearman(rng.permutation(a), b)) >= r0
    return (c + 1) / (n + 1)
A = [o for o in out if o["age"] > 35]
summ = dict(libraries=len({o["acc"] for o in out}), donors=len({o["donor"] for o in out}), groups=len(out), groups_over35=len(A),
            nuclei=int(sum(o["n"] for o in out)),
            noise_sd_control=float(np.median([o["f_control_sd"] for o in out])), median_se=float(np.median([o["se"] for o in out])))
fd = [o["f_dosage"] for o in A]; fr = [o["f_rna"] for o in A]; nz = [o["naive"] for o in A]; fc = [o["f_control_med"] for o in A]
summ["rho_dosage_rna"] = spearman(fd, fr); summ["p_dosage_rna"] = perm_p(fd, fr)
summ["rho_dosage_zero"] = spearman(fd, nz); summ["p_dosage_zero"] = perm_p(fd, nz)
summ["rho_control_rna"] = spearman(fc, fr); summ["p_control_rna"] = perm_p(fc, fr)
summ["mean_abs_diff_rna"] = float(np.mean(np.abs(np.array(fd) - np.array(fr)))); summ["mean_abs_diff_zero"] = float(np.mean(np.abs(np.array(fd) - np.array(nz))))
# age slopes per cell type (weighted by 1/se^2 + tau^2 not needed for summary: OLS with permutation)
sl = {}
for ct in cts:
    R = [o for o in out if o["cell_type"] == ct]
    if len(R) < 8: continue
    x = np.array([o["age"] for o in R]) / 10
    for key in ("f_dosage", "f_control_med", "f_rna", "naive"):
        y = np.array([o[key] for o in R]); b = np.polyfit(x, y, 1)[0]
        c = sum(abs(np.polyfit(rng.permutation(x), y, 1)[0]) >= abs(b) for _ in range(2000))
        sl.setdefault(ct, {})[key] = (100 * b, (c + 1) / 2001, len(R))
summ["age_slopes"] = {ct: {k: [round(v[0], 3), round(v[1], 4), v[2]] for k, v in d.items()} for ct, d in sl.items()}
summ["mean_over55"] = {ct: [round(100 * float(np.mean([o[k] for o in out if o["cell_type"] == ct and o["age"] >= 55])), 2) for k in ("f_dosage", "f_control_med", "f_rna", "naive")]
                       for ct in cts if sum(1 for o in out if o["cell_type"] == ct and o["age"] >= 55) >= 3}
json.dump(summ, open(os.path.join(ROOT, "tables/encode_dosage_summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
