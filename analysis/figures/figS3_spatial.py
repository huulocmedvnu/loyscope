# Supplementary Figure 3: Y dosage in Visium HD sections (spaty tiles of 16 x 16 bins of 8 um)
# a, b: tumour section with loss of Y, ratio of MSY reads to X-gametolog reads and epithelial share per tile
# c: MSY/X ratio across quintiles of epithelial share in each male section
# d: MSY reads per UMI in every section (female sections give the background of the probe set)
import os, glob, csv
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from common import C, BLACK, GREY, style, save, panel_titles_left, RES

SP = os.environ.get("LOY_SPATIAL", os.path.join(RES, "spatial"))
LAB = {"Colon_Cancer_P2": "Colorectal cancer, P2", "Colon_Normal_P3": "Normal colon, P3", "Colon_Normal_P5": "Normal colon, P5",
       "Prostate_Cancer_FFPE": "Prostate cancer", "Lung_Cancer": "Lung cancer", "Lung_Cancer_post_Xenium_Prime_5K_Experiment2": "Lung cancer, serial",
       "Tonsil_Fresh_Frozen": "Tonsil", "Breast_Cancer_Fresh_Frozen": "Breast cancer", "Kidney_FFPE": "Kidney", "Lymph_Node_FFPE": "Lymph node",
       "Colon_Cancer": "Colorectal cancer", "Pancreas": "Pancreas", "Colon_Cancer_P1": "Colorectal cancer, P1", "Colon_Cancer_P5": "Colorectal cancer, P5"}
def tiles(n):
    return [r for r in csv.DictReader(open(os.path.join(SP, f"Visium_HD_Human_{n}.t16.tiles.tsv")), delimiter="\t")]
DUP = {"Colon_Cancer"}  # Space Ranger 4 reanalysis of the P2 section
samples = sorted(os.path.basename(f).replace("Visium_HD_Human_", "").replace(".t16.tiles.tsv", "") for f in glob.glob(os.path.join(SP, "*.t16.tiles.tsv")))
samples = [n for n in samples if n not in DUP]
per_umi = {}
for n in samples:
    T = tiles(n); u = sum(float(r["umi"]) for r in T); y = sum(float(r["y"]) for r in T)
    per_umi[n] = y / u
male = [n for n in samples if per_umi[n] > 1e-5]
fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.12, vertical_spacing=0.16,
                    subplot_titles=("a  Colorectal cancer, MSY / X gametolog reads", "b  Same section, epithelial share",
                                    "c  MSY / X ratio by epithelial share", "d  MSY reads per UMI"))
# a, b: maps
T = [r for r in tiles("Colon_Cancer_P2") if float(r["y"]) + float(r["x"]) >= 20]
rr = np.array([float(r["tile_row"]) for r in T]); cc = np.array([float(r["tile_col"]) for r in T])
ratio = np.array([float(r["y"]) / max(float(r["x"]), 1) for r in T])
epi = np.array([float(r["epi"]) / max(float(r["epi"]) + float(r["stroma"]), 1) for r in T])
fig.add_trace(go.Scatter(x=cc, y=rr, mode="markers", marker=dict(symbol="square", size=4, color=np.clip(ratio, 0, 1.2), colorscale="RdBu", cmin=0, cmax=1.2,
                         colorbar=dict(title="MSY / X", len=0.38, y=0.80, x=0.445, thickness=10)), showlegend=False), 1, 1)
fig.add_trace(go.Scatter(x=cc, y=rr, mode="markers", marker=dict(symbol="square", size=4, color=epi, colorscale="Viridis", cmin=0, cmax=1,
                         colorbar=dict(title="epithelial", len=0.38, y=0.80, x=1.0, thickness=10)), showlegend=False), 1, 2)
for col in (1, 2):
    fig.update_xaxes(showticklabels=False, showline=False, zeroline=False, showgrid=False, ticks="", row=1, col=col)
    fig.update_yaxes(showticklabels=False, showline=False, zeroline=False, showgrid=False, ticks="", autorange="reversed", scaleanchor=f"x{col if col > 1 else ''}", row=1, col=col)
# c: quintiles of epithelial share
for i, n in enumerate(m for m in male if m != "Tonsil_Fresh_Frozen"):
    R = [r for r in tiles(n) if float(r["epi"]) + float(r["stroma"]) >= 50]
    if len(R) < 100: continue
    ef = np.array([float(r["epi"]) / (float(r["epi"]) + float(r["stroma"])) for r in R]); o = np.argsort(ef)
    xs, ys = [], []
    for q in np.array_split(o, 5):
        xs.append(ef[q].mean()); ys.append(sum(float(R[j]["y"]) for j in q) / sum(float(R[j]["x"]) for j in q))
    hot = n == "Colon_Cancer_P2"
    fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines+markers", line=dict(width=3 if hot else 1.5, color=C[1] if hot else C[(i + 2) % len(C)]),
                             marker=dict(size=7), name=LAB.get(n, n), legend="legend"), 2, 1)
fig.update_xaxes(title_text="Epithelial share of marker reads", range=[0, 1], row=2, col=1)
fig.update_yaxes(title_text="MSY / X gametolog reads", rangemode="tozero", row=2, col=1)
# d: reads per UMI
order = sorted(samples, key=lambda n: per_umi[n])
fig.add_trace(go.Bar(x=[max(per_umi[n], 1e-10) for n in order], y=[LAB.get(n, n) for n in order], orientation="h",
                     marker_color=[C[0] if per_umi[n] > 1e-5 else GREY for n in order], showlegend=False), 2, 2)
fig.update_xaxes(title_text="MSY reads per UMI", type="log", exponentformat="power", dtick=1, row=2, col=2)
style(fig, width=1000, height=900, legends=dict(legend=dict(x=0.0, y=-0.12, xanchor="left", yanchor="top", orientation="h", font=dict(size=11))))
fig.update_layout(margin=dict(l=60, r=40, t=50, b=120))
for col in (1, 2):
    fig.update_xaxes(showline=False, ticks="", showticklabels=False, zeroline=False, row=1, col=col)
    fig.update_yaxes(showline=False, ticks="", showticklabels=False, zeroline=False, row=1, col=col)
panel_titles_left(fig, [0.0, 0.56, 0.0, 0.56])
save(fig, "figS3_spatial")
print("written figS3_spatial", len(samples), "sections,", len(male), "male")
