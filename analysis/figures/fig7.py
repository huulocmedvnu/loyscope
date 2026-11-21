# Figure 7: consequences of LOY in cells and constraint across cell types
import glob
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from common import *

fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.22, column_widths=[0.5, 0.5],
                    subplot_titles=("a  Memory regulatory T cells without Y", "b  Y expression and loss across blood cell types"))
# a: genome-wide contrast in memory Treg (Immune Aging), excluding Y, PAR and Y lncRNAs; top genes by sign test
PAR = {"PLCXD1", "GTPBP6", "PPP2R3B", "SHOX", "CRLF2", "CSF2RA", "IL3RA", "SLC25A6", "ASMTL", "P2RY8", "AKAP17A", "ASMT", "DHRSX", "ZBED1", "CD99"}
Y = {"RPS4Y1", "ZFY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "TXLNGY", "NLGN4Y", "PRKY", "TMSB4Y"}
R = [r for r in rows(ext("iaging.de.tsv")) if r["cell_type"] == "memory regulatory T cell" and r["gene"] not in PAR | Y and not r["gene"].startswith(("TTTY", "LINC00278"))]
R = sorted([r for r in R if float(r["q_bh"]) < 0.1], key=lambda r: float(r["sign_p"]))
up = sorted([r for r in R if float(r["median_log2fc"]) > 0], key=lambda r: float(r["sign_p"]))[:10]
dn = sorted([r for r in R if float(r["median_log2fc"]) < 0], key=lambda r: float(r["sign_p"]))[:10]
sel = sorted(up + dn, key=lambda r: float(r["median_log2fc"]))
fig.add_trace(go.Bar(x=[float(r["median_log2fc"]) for r in sel], y=[r["gene"] for r in sel], orientation="h",
                     marker_color=[C[1] if float(r["median_log2fc"]) > 0 else C[0] for r in sel], showlegend=False), 1, 1)
fig.update_xaxes(title_text="log2 ratio, cells without Y / with Y (median of donors)", row=1, col=1)
# b: across blood cell types, LOY (donors >= 55) vs intact Y expression per UMI
T = [r for r in rows(os.path.join(TAB, "loy_table.tsv")) if r["dataset"] in ("aida", "iaging") and r["disease"] == "normal" and r["age"] != "None" and float(r["age"]) >= 55]
mu = {}
for ds in ("aida", "iaging"):
    for r in rows(ext(f"{ds}.fitamb.tsv")):
        if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 5: mu.setdefault((ds, r["cell_type"]), []).append(float(r["mu"]))
for j, ds in enumerate(("aida", "iaging")):
    don = {}
    for r in T:
        if r["dataset"] != ds: continue
        e = don.setdefault((r["cell_type"], r["donor"]), [0, 0.0]); n = int(r["n"]); e[0] += n; e[1] += n * max(0.0, float(r["loy"]))
    ct = {}
    for (c, d), (n, s) in don.items(): ct.setdefault(c, []).append(100 * s / n)
    pts = [(1e4 * st.median(mu[(ds, c)]), st.mean(v), c) for c, v in ct.items() if len(v) >= 5 and (ds, c) in mu]
    fig.add_trace(go.Scatter(x=[p[0] for p in pts], y=[p[1] for p in pts], mode="markers", marker=dict(size=9, color=C[[2, 3][j]]),
                             name={"aida": "AIDA", "iaging": "Immune Aging"}[ds], text=[p[2] for p in pts], legend="legend"), 1, 2)
    xs_, ys_ = [p[0] for p in pts], [p[1] for p in pts]
    rk = lambda v: [sorted(v).index(x) for x in v]; a_, b_ = rk(xs_), rk(ys_); n_ = len(pts)
    print(ds, "cell types", n_, "rho", round(1 - 6 * sum((a_[i] - b_[i]) ** 2 for i in range(n_)) / (n_ * (n_ * n_ - 1)), 2))
fig.update_xaxes(title_text="Y-gene UMIs per 10,000 UMIs in intact cells", type="log", row=1, col=2)
fig.update_yaxes(title_text="Cells without Y, donors aged 55 or older (%)", row=1, col=2)
style(fig, width=1000, height=520, legends=dict(legend=leg(0.99, 0.99, "right", "top")))
fig.update_layout(margin=dict(l=90, b=70))
panel_titles_left(fig, [-0.05, 0.57])
save(fig, "fig7")
