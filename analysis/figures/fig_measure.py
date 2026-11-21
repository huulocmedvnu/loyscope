# Figure 2 (main): measuring LOY in single cells. a, false loss for autosomal control genes by depth; b, change when depth is halved;
# c, model estimates for control genes by atlas; d, CELLxGENE versus author matrices. Full benchmarks are Supplementary Figs. 1 and 2.
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from common import *

NAMES = {"aida": "AIDA", "onek1k_geo": "OneK1K", "iaging": "Immune Aging", "blood": "PBMC India", "retina": "Retina/choroid", "tsapiens": "Tabula Sapiens",
         "ad_mssm": "Brain MSSM", "ad_radc": "Brain RADC", "ad_microglia": "Microglia", "ckd_sn": "Kidney KPMP", "covid": "COVID-19", "dcm": "Heart DCM", "hlca": "Lung HLCA"}
fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.2, vertical_spacing=0.2,
                    subplot_titles=("a  Genes that cannot be lost", "b  Sequencing depth halved", "c  Detection model, control genes", "d  Two count matrices of one cohort"))
bins = [(0, 1), (1, 2), (2, 5), (5, 10), (10, 1e9)]; lab = ["<1", "1-2", "2-5", "5-10", ">10"]
pn = {b: [] for b in bins}; pm = {b: [] for b in bins}
for fn in glob.glob(ext("*.fitamb.tsv")):
    for r in rows(fn):
        if r["sex"] != "male" or r["k"] == "NaN" or not (r["assay"].startswith("10x") or "GEO" in r["assay"]): continue
        e = float(r["mu"]) * float(r["mean_umi"])
        for b in bins:
            if b[0] <= e < b[1]: pn[b].append(float(r["pseudo_naive_med"])); pm[b].append(float(r["pseudo_f_med"]))
for i, (d, name) in enumerate([(pn, "zero-read rule"), (pm, "detection model")]):
    fig.add_trace(go.Scatter(x=lab, y=[100 * st.median(d[b]) for b in bins], mode="lines+markers", name=name, marker=dict(size=8, color=C[[0, 1][i]]), line=dict(color=C[[0, 1][i]]), legend="legend"), 1, 1)
fig.update_xaxes(title_text="Expected Y UMIs per intact cell", row=1, col=1); fig.update_yaxes(title_text="Cells called lost (%)", row=1, col=1)
for j, (m, name, c) in enumerate([("naive", "zero-read rule", C[0]), ("rf", "Random Forest", C[3]), ("zinb", "detection model", C[1])]):
    xs, ys = [], []
    for a, aname in [("aida", "AIDA blood"), ("retina", "Retina/choroid"), ("covid", "COVID-19 blood"), ("hlca", "Lung HLCA"), ("tsapiens", "Tabula Sapiens")]:
        by = {}
        for r in rows(ext(f"{a}.bench.tsv")): by.setdefault(r["gid"], {})[r["scenario"]] = r
        xs.append(st.median(100 * (float(g["thin0.5"][m]) - float(g["observed"][m])) for g in by.values() if "thin0.5" in g and float(g["observed"]["expY"]) >= 2)); ys.append(aname)
    fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers", marker=dict(symbol="diamond", size=11, color=c), name=name, legend="legend2"), 1, 2)
fig.update_xaxes(title_text="Change in cells called lost (points)", zeroline=True, zerolinecolor=GREY, row=1, col=2); fig.update_yaxes(autorange="reversed", row=1, col=2)
ys, med, p95 = [], [], []
for k, n in NAMES.items():
    fn = ext(f"{k}.fitamb.tsv")
    if not os.path.exists(fn): continue
    v = sorted(100 * float(r["pseudo_f_med"]) for r in rows(fn) if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 5 and (r["assay"].startswith("10x") or "GEO" in r["assay"]))
    if len(v) >= 10: ys.append(n); med.append(st.median(v)); p95.append(v[int(0.95 * len(v))])
fig.add_trace(go.Scatter(x=med, y=ys, mode="markers", marker=dict(symbol="diamond", size=10, color=C[1]), error_x=dict(type="data", symmetric=False, array=[b - a for a, b in zip(med, p95)], arrayminus=[0] * len(med), thickness=1, width=0), showlegend=False), 2, 1)
fig.update_xaxes(title_text="Estimated loss, median to 95th percentile (%)", range=[-0.3, 5], row=2, col=1); fig.update_yaxes(autorange="reversed", row=2, col=1)
def load(fn): return {(r["cell_type"], r["donor"]): 100 * max(0.0, float(r["f"]) - float(r["pseudo_f_med"])) for r in rows(fn) if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 1.5}
a, b = load(ext("onek1k.fitamb.tsv")), load(ext("onek1k_geo.fitamb.tsv")); k = [x for x in a if x in b]
fig.add_trace(go.Scatter(x=[a[x] for x in k], y=[b[x] for x in k], mode="markers", marker=dict(size=4, color=C[2], opacity=0.5), showlegend=False), 2, 2)
fig.add_trace(go.Scatter(x=[0, 40], y=[0, 40], mode="lines", line=dict(color=GREY, dash="dot"), showlegend=False), 2, 2)
fig.update_xaxes(title_text="OneK1K, CELLxGENE matrix (%)", range=[0, 40], row=2, col=2); fig.update_yaxes(title_text="OneK1K, GEO author matrix (%)", range=[0, 40], row=2, col=2)
style(fig, height=800, legends=dict(legend=leg(0.42, 0.99, "right", "top"), legend2=leg(1.0, 0.99, "right", "top")))
panel_titles_left(fig, [-0.14, 0.6, -0.14, 0.6])
fig.update_layout(margin=dict(l=140))
save(fig, "fig_measure")
