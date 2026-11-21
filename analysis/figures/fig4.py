# Figure 4: published LOY estimates and disease contrasts re-analysed with the detection model
import random, numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from common import *

fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.3, vertical_spacing=0.2, column_widths=[0.55, 0.45], row_heights=[0.42, 0.58],
                    specs=[[{}, {}], [{"colspan": 2}, None]],
                    subplot_titles=("a  Published estimates", "b  Calling rules in 25 studies", "c  Disease minus control (donors: disease vs control)"))
def wmean(fn, sel, key):
    R = [r for r in rows(fn) if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= 1.5 and sel(r)]
    n = sum(int(r["n"]) for r in R)
    if key == "corrected": return 100 * sum(int(r["n"]) * max(0.0, float(r["f"]) - float(r["pseudo_f_med"])) for r in R) / n
    return 100 * sum(int(r["n"]) * float(r["naive_zero"]) for r in R) / n
# a: three published claims: value reported in the paper, zero-read rule and model on the same public data
items = [("OneK1K classical monocytes", 18.6, ext("onek1k_geo.fitamb.tsv"), lambda r: r["cell_type"] == "CD14-positive monocyte"),
         ("Kidney, published study (one control)", 15.7, ext("gb2024.fitamb.tsv"), lambda r: r["donor"] == "GSE151302_Control1"),
         ("Microglia, Alzheimer disease", 28.1, ext("ad_microglia.fitamb.tsv"), lambda r: r["cell_type"] == "microglial cell" and r["disease"] == "Alzheimer disease"),
         ("Microglia, controls", 2.51, ext("ad_microglia.fitamb.tsv"), lambda r: r["cell_type"] == "microglial cell" and r["disease"] == "normal")]
labs = [i[0] for i in items]
for j, (name, key, c) in enumerate([("published", None, C[3]), ("zero-read rule", "naive", C[0]), ("detection model", "corrected", C[1])]):
    x = [it[1] if key is None else wmean(it[2], it[3], key) for it in items]
    fig.add_trace(go.Scatter(x=x, y=labs, mode="markers", marker=dict(symbol="diamond", size=11, color=c), name=name, legend="legend"), 1, 1)
fig.update_xaxes(title_text="Cells without Y (%)", range=[0, 32], row=1, col=1); fig.update_yaxes(autorange="reversed", row=1, col=1)
# b: donor-level age-adjusted differences with donor bootstrap CI, zero-read rule vs model
def contrast(fn, ct, case, thr=5, nb=1000):
    R = [r for r in rows(fn) if r["sex"] == "male" and r["k"] != "NaN" and float(r["mu"]) * float(r["mean_umi"]) >= thr and r["disease"] in (case, "normal") and r["cell_type"] == ct]
    dd = {}
    for r in R: dd.setdefault(r["donor"], []).append(r)
    pts = []
    for v in dd.values():
        a = age(v[0]["stage"]); n = sum(int(x["n"]) for x in v)
        if a is None: continue
        w = lambda k: sum(float(x[k]) * int(x["n"]) for x in v) / n
        pts.append((float(v[0]["disease"] == case), a, w("f") - w("pseudo_f_med"), w("naive_zero")))
    def est(p, j):
        X = np.array([[1, q[0], q[1]] for q in p]); y = np.array([q[j] for q in p]); return np.linalg.lstsq(X, y, rcond=None)[0][1] * 100
    rng = random.Random(7); out = []
    for j in (3, 2):  # zero-read rule first, then model
        b = est(pts, j); bs = []
        for _ in range(nb):
            s = [pts[rng.randrange(len(pts))] for _ in pts]
            if 0 < sum(q[0] for q in s) < len(s): bs.append(est(s, j))
        bs.sort(); out.append((b, bs[int(0.025 * len(bs))], bs[int(0.975 * len(bs))]))
    return out, int(sum(q[0] for q in pts)), int(len(pts) - sum(q[0] for q in pts))
cons = [("Kidney TAL, CKD (KPMP cells)", ext("kidney.fitamb.tsv"), "kidney loop of Henle thick ascending limb epithelial cell", "chronic kidney disease", 2),
        ("Kidney principal cells, CKD (KPMP nuclei)", ext("ckd_sn.fitamb.tsv"), "kidney collecting duct principal cell", "chronic kidney disease", 5),
        ("Kidney proximal tubule, AKI (KPMP nuclei)", ext("ckd_sn.fitamb.tsv"), "epithelial cell of proximal tubule", "acute kidney injury", 5),
        ("Heart, unannotated nuclei, DCM", ext("dcm.fitamb.tsv"), "unknown", "dilated cardiomyopathy", 2),
        ("Heart cardiomyocytes, DCM", ext("dcm.fitamb.tsv"), "cardiac muscle cell", "dilated cardiomyopathy", 2),
        ("Microglia, Alzheimer disease (fresh)", ext("ad_microglia.fitamb.tsv"), "microglial cell", "Alzheimer disease", 2),
        ("Microglia, Alzheimer disease (MSSM)", ext("ad_mssm.fitamb.tsv"), "microglial cell", "Alzheimer disease", 5)]
RC = {(r["dataset"], r["disease"], r["cell_type"]): r for r in rows(os.path.join(TAB, "robust_disease_contrasts.tsv"))}
ylab = []; res = []
for name, fn, ct, case, thr in cons:
    r = RC[(os.path.basename(fn).split(".")[0], case, ct)]
    ylab.append(f"{name} ({r['n_case']} vs {r['n_ctrl']})")
    res.append(((float(r["zero_read_diff"]), float(r["z_lo"]), float(r["z_hi"])), (float(r["model_diff"]), float(r["m_lo"]), float(r["m_hi"]))))
for j, (name, c, off) in enumerate([("zero-read rule", C[0], -0.15), ("detection model", C[1], 0.15)]):
    xs = [r[j][0] for r in res]
    fig.add_trace(go.Scatter(x=xs, y=[i + off for i in range(len(res))], mode="markers", marker=dict(symbol="diamond", size=10, color=c),
                             error_x=dict(type="data", symmetric=False, array=[r[j][2] - r[j][0] for r in res], arrayminus=[r[j][0] - r[j][1] for r in res], thickness=1.2, width=0),
                             name=name, legend="legend", showlegend=False), 2, 1)
fig.add_vline(x=0, line=dict(color=GREY), row=2, col=1)
fig.update_yaxes(tickvals=list(range(len(res))), ticktext=ylab, autorange="reversed", row=2, col=1)
fig.update_xaxes(title_text="Difference in cells without Y, age-adjusted (points)", row=2, col=1)
# c: calling rules of the 25 studies found by the literature search
cat = {}
for r in rows(os.path.join(os.path.dirname(__file__), "..", "..", "resources", "literature_search.tsv")):
    d = r["depth_handling"].lower()
    k = "ML classifier" if "ml classifier" in d else "depth-scaled threshold" if "scaled" in d else "UMI filter only" if "umi filter" in d else "none stated" if d.startswith("none") else "not reported"
    cat[k] = cat.get(k, 0) + 1
order = ["none stated", "UMI filter only", "depth-scaled threshold", "ML classifier", "not reported", "detection model"]
cat["detection model"] = 0
fig.add_trace(go.Bar(x=[cat.get(k, 0) for k in order], y=order, orientation="h", marker_color=[C[0], C[0], C[3], C[3], GREY, C[1]], text=[str(cat.get(k, 0)) for k in order], textposition="outside", showlegend=False), 1, 2)
fig.update_xaxes(title_text="Studies", range=[0, 16], row=1, col=2); fig.update_yaxes(autorange="reversed", row=1, col=2)
style(fig, width=1000, height=820, legends=dict(legend=leg(0.0, -0.1, "left", "top"), legend_orientation="h"))
fig.update_layout(margin=dict(l=300, b=90))
panel_titles_left(fig, [-0.3, 0.6, -0.3])
save(fig, "fig4")
