# Figure 6: origin of LOY: shared across tissues of a donor (haematopoietic) or tissue-specific; clonal overdispersion; lifespan
import random, math
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from common import *

T = rows(os.path.join(TAB, "loy_table.tsv"))
fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.18, vertical_spacing=0.2, specs=[[{}, {}], [{"colspan": 2}, None]],
                    subplot_titles=("a  Same donor, several tissues", "b  Cell lifespan", "c  Donors above 10%"))
# a: Tabula Sapiens donor share of variance across tissues
G = {"myeloid": ["monocyte", "macrophage", "dendritic", "neutrophil"], "lymphoid": ["T cell", "B cell", "natural killer", "plasma"],
     "stromal": ["fibroblast", "pericyte", "smooth muscle"], "endothelial": ["endothelial"]}
def r2(pp):
    m = st.mean(x[2] for x in pp); tot = sum((x[2] - m) ** 2 for x in pp)
    return sum(sum(1 for x in pp if x[0] == d) * (st.mean(x[2] for x in pp if x[0] == d) - m) ** 2 for d in {x[0] for x in pp}) / tot if tot > 0 else 0
xs, ys, ps = [], [], []
for g, w in G.items():
    cell = {}
    for r in T:
        if r["dataset"] != "tsapiens" or not any(k in r["cell_type"] for k in w): continue
        e = cell.setdefault((r["donor"], r["tissue"]), [0, 0.0]); n = int(r["n"]); e[0] += n; e[1] += n * float(r["loy"])
    pts = [(d, t, s / n) for (d, t), (n, s) in cell.items()]
    multi = {d for d, _, _ in pts if sum(1 for x in pts if x[0] == d) >= 3}; pts = [p for p in pts if p[0] in multi]
    if len(multi) < 3: continue
    o = r2(pts); rng = random.Random(1); lab = [p[0] for p in pts]; ge = 0
    for _ in range(5000): rng.shuffle(lab); ge += r2([(lab[i], pts[i][1], pts[i][2]) for i in range(len(pts))]) >= o
    xs.append(g); ys.append(100 * o); ps.append((ge + 1) / 5001)
fig.add_trace(go.Bar(x=xs, y=ys, marker_color=[C[0], C[3], C[4], C[1]][:len(xs)], text=[("P < 0.001" if p < 0.001 else f"P {p:.3f}") for p in ps], textposition="outside", showlegend=False), 1, 1)
fig.update_yaxes(title_text="Variance explained by donor (%)", range=[0, 100], row=1, col=1)
# b: share of donors aged >= 55 with LOY above 10% (cell types with >= 15 donors)
don = {}
for r in T:
    if r["disease"] != "normal" or r["age"] == "None" or float(r["age"]) < 55: continue
    e = don.setdefault((r["cell_type"], r["dataset"], r["donor"]), [0, 0.0]); n = int(r["n"]); e[0] += n; e[1] += n * max(0.0, float(r["loy"]))
ct = {}
for (c, ds, d), (n, s) in don.items(): ct.setdefault(c, []).append(s / n)
b = sorted([(c, 100 * sum(1 for x in v if x > 0.10) / len(v), len(v)) for c, v in ct.items() if len(v) >= 15], key=lambda x: -x[1])[:12]
fig.add_trace(go.Bar(x=[x[1] for x in b], y=[f"{x[0]} ({x[2]})" for x in b], orientation="h", marker_color=C[0], showlegend=False), 2, 1)
fig.update_xaxes(title_text="Donors with more than 10% of cells without Y (%)", row=2, col=1); fig.update_yaxes(autorange="reversed", row=2, col=1)
# c: class mean LOY (donors >= 55) vs lifespan (Sender and Milo 2021)
life = {}
for r in rows(os.path.join(os.path.dirname(__file__), "..", "..", "resources", "lifespan.tsv")):
    v = r["lifespan_days"] if r["lifespan_days"] != "None" else r["lifespan_rodent_days"]
    if v != "None": life[r["cell_type"]] = float(v)
MAP = [("Neurons", ["neuron", "interneuron", "projecting", "glutamatergic"]), ("Glial cells", ["astrocyte", "oligodendrocyte"]), ("Endothelial cells", ["endothelial"]),
       ("Monocytes", ["monocyte"]), ("Mature B cells", ["B cell"]), ("Mature T cells", ["T cell"]), ("Dermal fibroblasts", ["fibroblast"]),
       ("Alveolar macrophages", ["alveolar macrophage"]), ("Alveolar epithelial cells", ["pneumocyte", "alveolar type"]), ("Bronchial epithelial cells", ["club cell", "ciliated", "respiratory basal"])]
cls = {}
for (c, ds, d), (n, s) in don.items():
    k = next((k for k, w in MAP if any(x.lower() in c.lower() for x in w)), None)
    if k: cls.setdefault(k, {}).setdefault((ds, d), []).append(s / n)
px = [(k, life[k], 100 * st.mean(st.mean(v) for v in dv.values()), len(dv)) for k, dv in cls.items() if k in life and len(dv) >= 5]
fig.add_trace(go.Scatter(x=[p[1] for p in px], y=[p[2] for p in px], mode="markers+text", text=[p[0] for p in px], textposition=["top center" if p[0] not in ("Monocytes", "Neurons", "Mature B cells") else ("middle right" if p[0] == "Monocytes" else "bottom center" if p[0] == "Mature B cells" else "middle left") for p in px],
                         textfont=dict(size=11, color=BLACK), marker=dict(size=10, color=C[2]), showlegend=False), 1, 2)
fig.update_xaxes(title_text="Lifespan of the cell (days)", type="log", range=[0, 6], row=1, col=2); fig.update_yaxes(title_text="Cells without Y (%)", row=1, col=2)
style(fig, width=1000, height=900)
fig.update_layout(margin=dict(l=330, b=70))
panel_titles_left(fig, [-0.33, 0.59, -0.33])
save(fig, "fig6")
