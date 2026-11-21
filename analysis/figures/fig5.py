# Figure 5: where LOY occurs (donors aged >= 55, healthy tissue) and how it rises with age in blood (random-effects meta-analysis)
import math, numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from common import *

COMP = [("neural", ["neuron", "interneuron", "projecting", "glutamatergic", "astrocyte", "oligodendrocyte"], C[7]),
        ("vascular", ["endothelial", "pericyte", "smooth muscle", "leptomeningeal"], C[1]),
        ("myeloid", ["monocyte", "macrophage", "microglial", "dendritic", "neutrophil"], C[0]),
        ("lymphoid", ["T cell", "B cell", "natural killer", "plasma"], C[3]),
        ("epithelial and stromal", ["epithelial", "tubule", "fibroblast", "hepatocyte", "keratinocyte", "duct", "podocyte", "cardiac muscle"], C[4])]
def comp(c):
    for n, kw, col in COMP:
        if any(k.lower() in c.lower() for k in kw): return n, col
    return None, None
T = [r for r in rows(os.path.join(TAB, "loy_table.tsv")) if r["disease"] == "normal" and r["age"] != "None" and float(r["age"]) >= 55]
don = {}
for r in T:
    e = don.setdefault((r["cell_type"], r["dataset"], r["donor"]), [0, 0.0]); n = int(r["n"]); e[0] += n; e[1] += n * max(0.0, float(r["loy"]))
ct = {}
for (c, ds, d), (n, s) in don.items(): ct.setdefault(c, []).append(100 * s / n)
pick = [(c, v) for c, v in ct.items() if len(v) >= 10 and comp(c)[0]]
pick.sort(key=lambda x: -st.mean(x[1])); pick = pick[:14] + [p for p in pick[14:] if comp(p[0])[0] == "neural"][:4]
fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.32, vertical_spacing=0.17, column_widths=[0.5, 0.5],
                    specs=[[{"rowspan": 2}, {}], [None, {}]],
                    subplot_titles=("a  Men aged 55 or older", "b  Change with age, four blood cohorts", "c  Brain, MSSM cohort"))
shown = set()
for c, v in pick:
    n, col = comp(c); v = sorted(v); m = st.mean(v); se = st.pstdev(v) / math.sqrt(len(v))
    fig.add_trace(go.Scatter(x=[m], y=[f"{c} ({len(v)})"], mode="markers", marker=dict(symbol="diamond", size=10, color=col),
                             error_x=dict(type="data", array=[1.96 * se], thickness=1, width=0), name=n, legend="legend", showlegend=n not in shown), 1, 1)
    shown.add(n)
fig.update_xaxes(title_text="Cells without Y, mean of donors (%)", rangemode="tozero", row=1, col=1); fig.update_yaxes(autorange="reversed", row=1, col=1)
# b: DerSimonian-Laird pooled age slopes
CLS = [("classical monocyte", lambda c: ("CD14-positive" in c and "CD16-positive" not in c and "monocyte" in c) or c == "classical monocyte"),
       ("non-classical monocyte", lambda c: "CD16-positive" in c and "monocyte" in c), ("NK cell", lambda c: "natural killer" in c and "T cell" not in c),
       ("regulatory T cell", lambda c: "regulatory T" in c), ("naive B cell", lambda c: c == "naive B cell"), ("memory B cell", lambda c: "memory B cell" in c),
       ("CD4 T cell", lambda c: "CD4-positive" in c and "regulatory" not in c), ("CD8 T cell", lambda c: "CD8-positive" in c and "regulatory" not in c)]
cohorts = [("OneK1K", ext("onek1k_geo.fit.tsv")), ("AIDA", ext("aida.fit.tsv")), ("Immune Aging", ext("iaging.fit.tsv")), ("Indonesia blood", ext("blood.fit.tsv"))]
pooled = []
for cl, fnc in CLS:
    est = []
    for name, fn in cohorts:
        if not os.path.exists(fn): continue
        dd = {}
        for r in rows(fn):
            if r["sex"] == "male" and r["disease"] == "normal" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 1.5 and fnc(r["cell_type"]): dd.setdefault(r["donor"], []).append(r)
        pts = [(age(v[0]["stage"]), sum((float(x["f"]) - float(x["pseudo_f_med"])) * int(x["n"]) for x in v) / sum(int(x["n"]) for x in v)) for v in dd.values() if age(v[0]["stage"])]
        if len(pts) < 10: continue
        x = np.array([p[0] for p in pts]) / 10; y = np.array([p[1] for p in pts]); X = np.c_[np.ones_like(x), x]
        b = np.linalg.lstsq(X, y, rcond=None)[0]; s2 = ((y - X @ b) ** 2).sum() / (len(y) - 2); est.append((b[1], math.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])))
    if len(est) < 2: continue
    w = [1 / s ** 2 for _, s in est]; bf = sum(a * b for a, (b, _) in zip(w, est)) / sum(w)
    Q = sum(a * (b - bf) ** 2 for a, (b, _) in zip(w, est)); dfree = len(est) - 1; Cc = sum(w) - sum(a ** 2 for a in w) / sum(w); tau2 = max(0, (Q - dfree) / Cc)
    wr = [1 / (s ** 2 + tau2) for _, s in est]; br = sum(a * b for a, (b, _) in zip(wr, est)) / sum(wr); se = math.sqrt(1 / sum(wr))
    pooled.append((cl, 100 * br, 100 * se, len(est), max(0, (Q - dfree) / Q) * 100 if Q > 0 else 0))
fig.add_trace(go.Scatter(x=[p[1] for p in pooled], y=[f"{p[0]} (I2 {p[4]:.0f}%)" for p in pooled], mode="markers", marker=dict(symbol="diamond", size=10, color=C[0]),
                         error_x=dict(type="data", array=[1.96 * p[2] for p in pooled], thickness=1.2, width=0), showlegend=False), 1, 2)
fig.add_vline(x=0, line=dict(color=GREY), row=1, col=2)
fig.update_xaxes(title_text="Change per decade of age (points)", row=1, col=2); fig.update_yaxes(autorange="reversed", row=1, col=2)
# c: donor-level LOY against age in the MSSM brain cohort, with least-squares lines
import numpy as np
BR = [("vascular leptomeningeal cell", "VLMC", C[1]), ("endothelial cell", "endothelium", C[4]), ("microglial cell", "microglia", C[0]), ("L2/3-6 intratelencephalic projecting glutamatergic neuron", "neurons", C[7])]
mssm = [r for r in rows(ext("ad_mssm.fitamb.tsv")) if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 5]
for ct, lab, col in BR:
    dd = {}
    for r in mssm:
        if r["cell_type"] == ct: dd.setdefault(r["donor"], []).append(r)
    pts = [(age(v[0]["stage"]), 100 * sum((float(x["f"]) - float(x["pseudo_f_med"])) * int(x["n"]) for x in v) / sum(int(x["n"]) for x in v)) for v in dd.values() if age(v[0]["stage"])]
    if len(pts) < 5: continue
    x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts]); b = np.polyfit(x, y, 1)
    fig.add_trace(go.Scatter(x=x, y=y, mode="markers", marker=dict(size=5, color=col, opacity=0.55), name=f"{lab} ({len(pts)})", legend="legend2"), 2, 2)
    xx = np.array([x.min(), x.max()]); fig.add_trace(go.Scatter(x=xx, y=b[0] * xx + b[1], mode="lines", line=dict(color=col, width=2.5), showlegend=False), 2, 2)
fig.update_xaxes(title_text="Age (years)", row=2, col=2); fig.update_yaxes(title_text="Cells without Y (%)", range=[-2, 36], row=2, col=2)
style(fig, width=1000, height=860, legends=dict(legend=leg(0.0, -0.08, "left", "top"), legend_orientation="h", legend2=leg(0.64, 0.40, "left", "top")))
fig.update_layout(margin=dict(l=280, b=90))
panel_titles_left(fig, [-0.33, 0.62, 0.62])
save(fig, "fig5")
import json; json.dump(pooled, open(os.path.join(TAB, "age_slope_meta.json"), "w"))
