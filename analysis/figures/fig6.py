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
# c: share of donors aged >= 55 with more than 10% of cells without Y, observed and under a common-fraction null
#    (analysis/robust_stats.py, Supplementary Table 4c); twelve cell types with the highest observed share
O = sorted(rows(os.path.join(TAB, "robust_overdispersion.tsv")), key=lambda r: -float(r["share_above_10pct"]))[:12]
lab = [f"{r['cell_type']} ({r['donors']})" for r in O]
fig.add_trace(go.Bar(x=[float(r["share_above_10pct"]) for r in O], y=lab, orientation="h", marker_color=C[0], name="observed", legend="legend"), 2, 1)
fig.add_trace(go.Scatter(x=[float(r["null_share"]) for r in O], y=lab, mode="markers", marker=dict(symbol="line-ns-open", size=14, color=BLACK, line=dict(width=2)),
                         name="expected from measurement error", legend="legend"), 2, 1)
fig.update_xaxes(title_text="Donors with more than 10% of cells without Y (%)", row=2, col=1); fig.update_yaxes(autorange="reversed", row=2, col=1)
# b: class mean LOY (donors >= 55) vs lifespan (Sender and Milo 2021), classes and values from analysis/turnover.py
L = rows(os.path.join(TAB, "lifespan_classes.tsv"))
# label placement: (x anchor, x shift, y shift) in pixels, chosen so that no two labels overlap
OFF = {"Endothelial cells": ("center", 0, 14), "Monocytes": ("left", 9, 9), "Dermal fibroblasts": ("center", 0, -15), "Mature B cells": ("center", 0, -15),
       "Mature T cells": ("left", 8, 0), "Alveolar macrophages": ("left", 8, 0), "Alveolar epithelial cells": ("center", 0, 14),
       "Bronchial epithelial cells": ("center", 0, -15), "Glial cells": ("right", -7, 9), "Cardiomyocytes": ("center", 0, 14), "Neurons": ("center", 0, -14)}
fig.add_trace(go.Scatter(x=[float(r["lifespan_days"]) for r in L], y=[float(r["mean_loy_pct"]) for r in L], mode="markers", marker=dict(size=10, color=C[2]), showlegend=False), 1, 2)
for r in L:
    xa, dx, dy = OFF.get(r["cell_class"], ("center", 0, 14))
    fig.add_annotation(x=math.log10(float(r["lifespan_days"])), y=float(r["mean_loy_pct"]), text=r["cell_class"], showarrow=False, xanchor=xa, xshift=dx, yshift=dy,
                       font=dict(size=11, color=BLACK), row=1, col=2)
fig.update_xaxes(title_text="Lifespan of the cell (days)", type="log", range=[0, 6.3], row=1, col=2); fig.update_yaxes(title_text="Cells without Y (%)", row=1, col=2)
style(fig, width=1000, height=900, legends=dict(legend=leg(0.99, 0.30, "right", "top")))
fig.update_layout(margin=dict(l=330, b=70))
panel_titles_left(fig, [-0.33, 0.59, -0.33])
save(fig, "fig6")
