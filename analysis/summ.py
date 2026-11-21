# summarize one or more <name>.fit.tsv files (numbers computed in Rust)
import csv, re, sys, statistics as st, random
def age(s):
    m = re.match(r"(\d+)-year-old", s); 
    if m: return float(m.group(1))
    m = re.match(r"(\d+)-(\d+) year-old", s)
    if m: return (float(m.group(1)) + float(m.group(2))) / 2
    DEC = {"third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75, "ninth": 85, "tenth": 95}
    m = re.match(r"(\w+) decade", s)
    return DEC.get(m.group(1)) if m else None
for fn in sys.argv[1:]:
    R = list(csv.DictReader(open(fn), delimiter="\t"))
    M0 = [r for r in R if r["sex"] == "male"]; F = [r for r in R if r["sex"] == "female"]
    # QC (amendment, post hoc): Y genes unmeasured (all-zero Y in a male group = zero-filled genes after integration);
    # identifiability: expected Y UMI per intact cell mu*mean_umi >= 2
    nz = [r for r in M0 if r["k"] == "NaN"]
    M = [r for r in M0 if r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 5]
    print(f"   QC: male groups {len(M0)}, all-zero Y {len(nz)} (donors {len({r['donor'] for r in nz})}), low expected Y {len(M0)-len(nz)-len(M)}")
    name = fn.split("/")[-1].replace(".fit.tsv", "")
    if not M: print(name, "no male groups"); continue
    pf = sorted(float(r["pseudo_f_med"]) for r in M)
    print(f"== {name}: male groups {len(M)}, female groups {len(F)}, male donors {len({r['donor'] for r in M})}")
    print(f"  C1 pseudo-Y f median {st.median(pf):.4f} p95 {pf[int(.95*len(pf))]:.4f} | naive pseudo zero median {st.median(float(r['pseudo_naive_med']) for r in M):.3f}")
    if F: print(f"  C2 female Y f median {st.median(float(r['f']) for r in F):.3f}")
    ex = sorted(float(r["f"]) - float(r["pseudo_f_med"]) for r in M)
    for t in (2, 5, 10):
        Mt = [r for r in M if float(r["mu"]) * float(r["mean_umi"]) >= t]
        if Mt: print(f"  expY>={t}: groups {len(Mt)}, pseudo f median {st.median(float(r['pseudo_f_med']) for r in Mt):.4f}, Y f median {st.median(float(r['f']) for r in Mt):.4f}, excess median {st.median(float(r['f'])-float(r['pseudo_f_med']) for r in Mt):+.4f}")
    print(f"  male Y: naive zero median {st.median(float(r['naive_zero']) for r in M):.3f}, corrected f median {st.median(float(r['f']) for r in M):.3f}")
    ct = {}
    for r in M: ct.setdefault(r["cell_type"], []).append(r)
    rows = []
    for c, g in ct.items():
        # donor-level: average f over groups of the same donor (n-weighted)
        dd = {}
        for r in g: dd.setdefault(r["donor"], []).append(r)
        don = [(age(v[0]["stage"]), sum((float(x["f"]) - float(x["pseudo_f_med"])) * int(x["n"]) for x in v) / sum(int(x["n"]) for x in v)) for v in dd.values()]
        don = [x for x in don if x[0] is not None]
        fm = st.median(float(r["f"]) - float(r["pseudo_f_med"]) for r in g)
        sl, p = float("nan"), float("nan")
        if len(don) >= 8 and max(x[0] for x in don) - min(x[0] for x in don) >= 20:
            def slope(d):
                xs = [a for a, _ in d]; ys = [b for _, b in d]; mx = st.mean(xs); my = st.mean(ys)
                return sum((a - mx) * (b - my) for a, b in d) / sum((a - mx) ** 2 for a in xs)
            sl = slope(don); ys = [b for _, b in don]; xs = [a for a, _ in don]; random.seed(1); ge = 0
            for _ in range(2000):
                random.shuffle(ys); ge += slope(list(zip(xs, ys))) >= sl
            p = (ge + 1) / 2001
        rows.append((c, len(dd), fm, sl, p))
    rows.sort(key=lambda x: -x[2])
    for c, nd, fm, sl, p in rows[:12]:
        print(f"  {c[:38]:38s} donors {nd:3d} excess {fm:+.3f} slope/10y {sl*10:+.3f} P {p:.3f}")
