# summarize loybench output: spike-in recovery and depth-thinning shift per method (numbers from Rust)
import csv, sys, statistics as st
for fn in sys.argv[1:]:
    R = list(csv.DictReader(open(fn), delimiter="\t")); by = {}
    for r in R: by.setdefault(r["gid"], {})[(r["scenario"], r["f_true"])] = r
    print(f"== {fn}: groups {len(by)}")
    for lo, hi in [(0, 2), (2, 5), (5, 1e9)]:
        G = [g for g in by.values() if lo <= float(g[("observed", "0")]["expY"]) < hi]
        if not G: continue
        o = lambda m: st.median(float(g[("observed", "0")][m]) for g in G)
        print(f"  expY [{lo},{hi}) groups {len(G)} | observed median: naive {o('naive'):.3f} rf {o('rf'):.3f} zinb {o('zinb'):.3f}")
        for ft in ["0.02", "0.05", "0.1", "0.2"]:
            s = []
            for m in ["naive", "rf", "zinb"]:
                rec = [(float(g[("spike", ft)][m]) - float(g[("observed", "0")][m])) / max(1e-9, 1 - float(g[("observed", "0")][m])) for g in G]
                s.append(f"{m} {st.median(rec):.3f}")
            print(f"    spike {ft}: recovered " + ", ".join(s))
        sh = [f"{m} {st.median(float(g[('thin0.5','0')][m]) - float(g[('observed','0')][m]) for g in G):+.3f}" for m in ["naive", "rf", "zinb"]]
        print("    half depth shift: " + ", ".join(sh))
