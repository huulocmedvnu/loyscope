# disease vs normal, males, donor-level; corrected excess vs naive zero; age-adjusted OLS (numbers from Rust loyfit)
import csv, re, sys, statistics as st, random
import numpy as np
DEC = {"third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75, "ninth": 85, "tenth": 95}
def age(s):
    m = re.match(r"(\d+)-year-old", s)
    if m: return float(m.group(1))
    m = re.match(r"(\d+)-(\d+) year-old", s)
    if m: return (float(m.group(1)) + float(m.group(2))) / 2
    m = re.match(r"(\w+) decade", s)
    return DEC.get(m.group(1)) if m else None
fn, case = sys.argv[1], sys.argv[2]; thr = float(sys.argv[3]) if len(sys.argv) > 3 else 5
R = [r for r in csv.DictReader(open(fn), delimiter="\t") if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= thr and r["disease"] in (case, "normal")]
print(f"{fn} {case} vs normal, expY>={thr}: groups {len(R)}, donors case {len({r['donor'] for r in R if r['disease']==case})} normal {len({r['donor'] for r in R if r['disease']=='normal'})}")
ct = {}
for r in R: ct.setdefault(r["cell_type"], []).append(r)
for c, g in sorted(ct.items(), key=lambda x: -len(x[1])):
    dd = {}
    for r in g: dd.setdefault(r["donor"], []).append(r)
    rows = []
    for v in dd.values():
        a = age(v[0]["stage"]); n = sum(int(x["n"]) for x in v)
        if a is None: continue
        w = lambda k: sum(float(x[k]) * int(x["n"]) for x in v) / n
        rows.append((1.0 if v[0]["disease"] == case else 0.0, a, w("f") - w("pseudo_f_med"), w("naive_zero"), w("mean_umi")))
    nc = sum(r[0] for r in rows); nn = len(rows) - nc
    if nc < 4 or nn < 4: continue
    X = np.array([[1, r[0], r[1]] for r in rows])
    def coef(j):
        y = np.array([r[j] for r in rows]); b = np.linalg.lstsq(X, y, rcond=None)[0]
        # permutation of disease labels
        rng = random.Random(1); lab = [r[0] for r in rows]; ge = 0
        for _ in range(2000):
            rng.shuffle(lab); Xp = X.copy(); Xp[:, 1] = lab; bp = np.linalg.lstsq(Xp, y, rcond=None)[0]; ge += abs(bp[1]) >= abs(b[1])
        return b[1], (ge + 1) / 2001
    bc, pc = coef(2); bn, pn = coef(3); bu, pu = coef(4)
    print(f"  {c[:34]:34s} case {int(nc):2d} ctrl {int(nn):2d} | corrected diff {bc:+.3f} P {pc:.3f} | naive diff {bn:+.3f} P {pn:.3f} | UMI diff {bu:+.0f} P {pu:.3f}")
