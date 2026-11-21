# Figure 3: model-called LOY cells have lost the chromosome (PAR dosage, gametologs, chromatin)
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from common import *

PAR1 = ["CD99", "SLC25A6", "AKAP17A", "ASMTL", "P2RY8", "ZBED1", "DHRSX", "CSF2RA", "GTPBP6", "PLCXD1", "PPP2R3B"]
GAM = ["DDX3X", "KDM5C", "KDM6A", "EIF1AX", "ZFX", "RPS4X", "USP9X"]
SETS = {"aida": "AIDA blood", "retina": "Retina/choroid", "ad_microglia": "Microglia"}
fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.1, column_widths=[0.36, 0.30, 0.34],
                    subplot_titles=("a  PAR1 genes (X and Y copies)", "b  X gametologs", "c  Kidney DN4, one donor"))
def gene_stats(genes, col):
    out = {}
    for k in SETS:
        fn = ext(f"{k}.cellde.tsv")
        if not os.path.exists(fn): continue
        R = [r for r in rows(fn) if float(r["mean_ypos_cp10k"]) >= 0.5 and r[col] != "nan"]
        for g in genes:
            v = [float(r[col]) for r in R if r["gene"] == g]
            if len(v) >= 3: out.setdefault(g, []).extend(v)
    return out
for col, (genes, cidx) in enumerate([(PAR1, 1), (GAM, 2)]):
    loy = gene_stats(genes, "log2fc_loy"); ctl = gene_stats(genes, "log2fc_pseudo")
    gs = [g for g in genes if g in loy]
    for j, (d, name, c) in enumerate([(loy, "model-called LOY cells", C[1]), (ctl, "pseudo-zero cells (control)", C[0])]):
        y = [st.median(d.get(g, [float("nan")])) for g in gs]
        q1 = [sorted(d[g])[len(d[g]) // 4] if g in d else None for g in gs]; q3 = [sorted(d[g])[3 * len(d[g]) // 4] if g in d else None for g in gs]
        fig.add_trace(go.Scatter(x=y, y=[g for g in gs], mode="markers", marker=dict(symbol="diamond", size=9, color=c),
                                 error_x=dict(type="data", symmetric=False, array=[b - a for a, b in zip(y, q3)], arrayminus=[a - b for a, b in zip(y, q1)], thickness=1, width=0),
                                 name=name, legend="legend", showlegend=(col == 0)), 1, col + 1)
    fig.add_vline(x=-1, line=dict(color=GREY, dash="dot"), row=1, col=col + 1); fig.add_vline(x=0, line=dict(color=GREY), row=1, col=col + 1)
    fig.update_xaxes(title_text="log2 ratio, LOY / Y+ cells", range=[-2, 1], row=1, col=col + 1)
    fig.update_yaxes(autorange="reversed", row=1, col=col + 1)
# c: sample-level estimates for the same donor (ATAC fragments vs RNA calls); values from fragy + loyfit + tenx outputs
import json
v = json.load(open(os.path.join(TAB, "dn4_estimates.json")))
labs = ["ATAC fragments", "RNA, model", "RNA, zero-read", "RNA, trough"]
vals = [v["atac"], v["model"], v["naive"], v["trough"]]
fig.add_trace(go.Bar(x=vals, y=labs, orientation="h", marker_color=[C[2], C[1], C[0], C[3]], text=[f"{x:.1f}" for x in vals], textposition="outside", showlegend=False), 1, 3)
fig.update_xaxes(title_text="Cells without Y (%)", range=[0, 36], row=1, col=3)
fig.update_yaxes(autorange="reversed", row=1, col=3)
style(fig, width=1000, height=460, legends=dict(legend=leg(0.0, -0.22, "left", "top"), legend_orientation="h"))
fig.update_layout(margin=dict(l=90, b=110))
panel_titles_left(fig, [0.0, 0.41, 0.72])
save(fig, "fig3")
