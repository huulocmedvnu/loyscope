# Figure 2: the detection model is calibrated, depth-invariant and reproduces on author matrices
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from common import *

NAMES = {k: ATLAS_NAMES[k] for k in canonical()}
fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.17, vertical_spacing=0.2,
                    subplot_titles=("a  Pseudo-Y genes", "b  Observed groups", "c  Depth halved", "d  CELLxGENE and author matrices"))
# a: model estimate on pseudo-Y sets (true value 0): median and 95th percentile over male groups with expected Y >= 5
ys, med, p95 = [], [], []
for k, n in NAMES.items():
    fn = ext(f"{k}.fitamb.tsv")
    if not os.path.exists(fn): continue
    v = sorted(100 * float(r["pseudo_f_med"]) for r in rows(fn) if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 5 and (r["assay"].startswith("10x") or "GEO" in r["assay"]))
    if len(v) < 10: continue
    ys.append(n); med.append(st.median(v)); p95.append(v[int(0.95 * len(v))])
fig.add_trace(go.Scatter(x=med, y=ys, mode="markers", marker=dict(symbol="diamond", size=10, color=C[1]), error_x=dict(type="data", symmetric=False, array=[b - a for a, b in zip(med, p95)], arrayminus=[0] * len(med), thickness=1, width=0), name="median, 95th percentile", legend="legend"), 1, 1)
fig.update_xaxes(title_text="Estimated loss for pseudo-Y genes (%)", range=[-0.3, 5], row=1, col=1)
fig.update_yaxes(autorange="reversed", tickmode="array", tickvals=list(NAMES.values()), ticktext=list(NAMES.values()), tickfont=dict(size=10), row=1, col=1)
# b: zero-read rule vs model, male groups, expected Y >= 5, all atlases
xn, ym = [], []
for k in NAMES:
    fn = ext(f"{k}.fitamb.tsv")
    if not os.path.exists(fn): continue
    for r in rows(fn):
        if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 5:
            xn.append(100 * float(r["naive_zero"])); ym.append(100 * max(0.0, float(r["f"]) - float(r["pseudo_f_med"])))
fig.add_trace(go.Scatter(x=xn, y=ym, mode="markers", marker=dict(size=3, color=C[0], opacity=0.35), showlegend=False), 1, 2)
fig.add_trace(go.Scatter(x=[0, 60], y=[0, 60], mode="lines", line=dict(color=GREY, dash="dot"), showlegend=False), 1, 2)
fig.update_xaxes(title_text="Zero-read rule (%)", range=[0, 60], row=1, col=2); fig.update_yaxes(title_text="Detection model (%)", range=[0, 60], row=1, col=2)
# c: model estimate at full vs half depth (bench), groups with expected Y >= 2
xs, yh = [], []
for k in ["aida", "retina", "covid", "hlca", "tsapiens", "ad_mssm", "ckd_sn"]:
    fn = ext(f"{k}.bench.tsv")
    if not os.path.exists(fn): continue
    by = {}
    for r in rows(fn): by.setdefault(r["gid"], {})[r["scenario"]] = r
    for g in by.values():
        if "thin0.5" in g and float(g["observed"]["expY"]) >= 2: xs.append(100 * float(g["observed"]["zinb"])); yh.append(100 * float(g["thin0.5"]["zinb"]))
fig.add_trace(go.Scatter(x=xs, y=yh, mode="markers", marker=dict(size=3, color=C[1], opacity=0.35), showlegend=False), 2, 1)
fig.add_trace(go.Scatter(x=[0, 40], y=[0, 40], mode="lines", line=dict(color=GREY, dash="dot"), showlegend=False), 2, 1)
fig.update_xaxes(title_text="Full depth (%)", range=[0, 40], row=2, col=1); fig.update_yaxes(title_text="Half depth (%)", range=[0, 40], row=2, col=1)
# d: OneK1K CELLxGENE vs GEO author matrices, donor x cell type, expected Y >= 1.5
def load(fn):
    return {(r["cell_type"], r["donor"]): 100 * max(0.0, float(r["f"]) - float(r["pseudo_f_med"])) for r in rows(fn)
            if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 1.5}
a, b = load(ext("onek1k.fitamb.tsv")), load(ext("onek1k_geo.fitamb.tsv"))
k = [x for x in a if x in b]
fig.add_trace(go.Scatter(x=[a[x] for x in k], y=[b[x] for x in k], mode="markers", marker=dict(size=4, color=C[2], opacity=0.5), showlegend=False), 2, 2)
fig.add_trace(go.Scatter(x=[0, 40], y=[0, 40], mode="lines", line=dict(color=GREY, dash="dot"), showlegend=False), 2, 2)
fig.update_xaxes(title_text="CELLxGENE matrix (%)", range=[0, 40], row=2, col=2); fig.update_yaxes(title_text="GEO author matrix (%)", range=[0, 40], row=2, col=2)
style(fig, height=960, legends=dict(legend=leg(0.41, 0.99, "right", "top")))
panel_titles_left(fig, [0.0, 0.585, 0.0, 0.585])
save(fig, "fig2")
