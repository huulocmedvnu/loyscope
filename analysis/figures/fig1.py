# Figure 1: current single-cell LOY calls track sequencing depth (numbers from loyfit / loybench outputs)
import csv, glob, os, statistics as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

RES = os.environ.get("LOY_RESULTS", "results"); EXT = os.path.join(RES, "extract")
OUT = os.path.join(RES, "figures"); os.makedirs(OUT, exist_ok=True)
C = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3"]
BLACK = "#111111"
ATLASES = {"aida": "AIDA blood", "retina": "Retina/choroid", "covid": "COVID-19 blood", "hlca": "Lung (HLCA)", "tsapiens": "Tabula Sapiens"}

def rows(fn): return list(csv.DictReader(open(fn), delimiter="\t"))

fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.2, vertical_spacing=0.2,
                    subplot_titles=("a  Autosomal pseudo-Y genes", "b  Depth halved", "c  Random Forest and zero-read calls", "d  Spike-in loss"))
# a: pseudo-Y apparent loss under the zero-read rule vs model, by expected Y UMI per intact cell
bins = [(0, 1), (1, 2), (2, 5), (5, 10), (10, 20), (20, 1e9)]; lab = ["<1", "1-2", "2-5", "5-10", "10-20", ">20"]
pn = {b: [] for b in bins}; pm = {b: [] for b in bins}
for fn in glob.glob(os.path.join(EXT, "*.fitamb.tsv")):
    for r in rows(fn):
        if r["sex"] != "male" or r["k"] == "NaN" or not (r["assay"].startswith("10x") or "GEO" in r["assay"]): continue
        e = float(r["mu"]) * float(r["mean_umi"])
        for b in bins:
            if b[0] <= e < b[1]: pn[b].append(float(r["pseudo_naive_med"])); pm[b].append(float(r["pseudo_f_med"]))
for i, (d, name) in enumerate([(pn, "zero-read rule"), (pm, "detection model")]):
    y = [100 * st.median(d[b]) if d[b] else None for b in bins]
    fig.add_trace(go.Scatter(x=lab, y=y, mode="lines+markers", name=name, marker=dict(size=8, color=C[i]), line=dict(color=C[i]), legend="legend"), 1, 1)
fig.update_xaxes(title_text="Expected Y UMIs per intact cell", row=1, col=1)
fig.update_yaxes(title_text="Cells called lost (%)", row=1, col=1)
# b: shift in estimated LOY when depth is halved (median over groups), per atlas and method
meths = [("naive", "zero-read rule"), ("rf", "Random Forest"), ("zinb", "detection model")]
for j, (m, name) in enumerate(meths):
    xs, ys = [], []
    for a, aname in ATLASES.items():
        fn = os.path.join(EXT, f"{a}.bench.tsv")
        if not os.path.exists(fn): continue
        by = {}
        for r in rows(fn): by.setdefault(r["gid"], {})[r["scenario"] + r["f_true"]] = r
        sh = [100 * (float(g["thin0.50"][m]) - float(g["observed0"][m])) for g in by.values() if "thin0.50" in g and float(g["observed0"]["expY"]) >= 2]
        if sh: xs.append(aname); ys.append(st.median(sh))
    fig.add_trace(go.Scatter(y=xs, x=ys, mode="markers", name=name, marker=dict(symbol="diamond", size=11, color=C[[0, 3, 1][j]]), legend="legend2"), 1, 2)
fig.update_xaxes(title_text="Change in LOY when depth is halved (points)", row=1, col=2)
fig.update_yaxes(autorange="reversed", row=1, col=2)
fig.update_xaxes(zeroline=True, zerolinecolor="#999999", row=1, col=2)
# c: RF vs zero-read estimate per group (AIDA)
fn = os.path.join(EXT, "aida.bench.tsv")
obs = [r for r in rows(fn) if r["scenario"] == "observed"]
fig.add_trace(go.Scatter(x=[100 * float(r["naive"]) for r in obs], y=[100 * float(r["rf"]) for r in obs], mode="markers",
                         marker=dict(size=4, color=C[3], opacity=0.5), name="AIDA groups", showlegend=False), 2, 1)
fig.add_trace(go.Scatter(x=[0, 100], y=[0, 100], mode="lines", line=dict(color="#999999", dash="dot"), showlegend=False), 2, 1)
fig.update_xaxes(title_text="Zero-read rule (%)", range=[0, 60], row=2, col=1)
fig.update_yaxes(title_text="Random Forest (%)", range=[0, 60], row=2, col=1)
# d: recovery of spiked-in loss (groups with expected Y >= 5), pooled over atlases
for j, (m, name) in enumerate(meths):
    xs, ys, lo, hi = [], [], [], []
    for ft in ["0.02", "0.05", "0.1", "0.2"]:
        rec = []
        for a in ATLASES:
            fn = os.path.join(EXT, f"{a}.bench.tsv")
            if not os.path.exists(fn): continue
            by = {}
            for r in rows(fn): by.setdefault(r["gid"], {})[r["scenario"] + r["f_true"]] = r
            for g in by.values():
                if "spike" + ft in g and float(g["observed0"]["expY"]) >= 5:
                    o = float(g["observed0"][m]); rec.append(100 * (float(g["spike" + ft][m]) - o) / max(1e-9, 1 - o))
        rec.sort(); xs.append(100 * float(ft)); ys.append(st.median(rec)); lo.append(rec[len(rec) // 4]); hi.append(rec[3 * len(rec) // 4])
    fig.add_trace(go.Scatter(x=[x + (j - 1) * 0.5 for x in xs], y=ys, mode="markers", marker=dict(symbol="diamond", size=9, color=C[[0, 3, 1][j]]),
                             error_y=dict(type="data", symmetric=False, array=[h - y for h, y in zip(hi, ys)], arrayminus=[y - l for y, l in zip(ys, lo)], thickness=1, width=0),
                             name=name, legend="legend3"), 2, 2)
fig.add_trace(go.Scatter(x=[0, 21], y=[0, 21], mode="lines", line=dict(color="#999999", dash="dot"), showlegend=False), 2, 2)
fig.update_xaxes(title_text="Spiked-in loss (%)", row=2, col=2)
fig.update_yaxes(title_text="Recovered loss (%)", row=2, col=2)
fig.update_layout(template="plotly_white", width=900, height=780, font=dict(family="Arial", size=13, color=BLACK), barmode="group",
                  margin=dict(l=70, r=20, t=50, b=60),
                  legend=dict(x=0.25, y=0.98, xanchor="left", yanchor="top", bgcolor="rgba(255,255,255,0)"),
                  legend2=dict(x=0.99, y=0.99, xanchor="right", yanchor="top", bgcolor="rgba(255,255,255,0)"),
                  legend3=dict(x=0.58, y=0.40, xanchor="left", yanchor="top", bgcolor="rgba(255,255,255,0)"))
fig.update_xaxes(showgrid=False, showline=True, linecolor=BLACK, ticks="outside", tickfont=dict(color=BLACK), title_font=dict(color=BLACK))
fig.update_yaxes(showgrid=False, showline=True, linecolor=BLACK, ticks="outside", tickfont=dict(color=BLACK), title_font=dict(color=BLACK), zeroline=True, zerolinecolor="#999999")
for ann in fig.layout.annotations: ann.update(x=0.0 if ann.x < 0.5 else 0.565, xanchor="left", font=dict(size=14, color=BLACK))
fig.write_image(os.path.join(OUT, "fig1.png"), scale=3); fig.write_html(os.path.join(OUT, "fig1.html"))
print("written", OUT)
