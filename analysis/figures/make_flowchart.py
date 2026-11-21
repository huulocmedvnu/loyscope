"""Figure 1: study design as an SVG schematic, exported to PNG (600 dpi) and PDF.

Rows follow the order of the Results. Each row: data (left), a pictogram of the comparison (middle) and a small
plot of the readout (right). Readout marks use the estimates reported in the Results (results/tables/numbers.json).
Canvas units: 10 per mm, 180 mm wide, so 8 pt is about 28 units.
"""
import json, math, os, subprocess
from common import TAB, OUT

N = json.load(open(os.path.join(TAB, "numbers.json")))
W = 1800
INK, GREY, LINE = "#111111", "#8C96A8", "#DCE1EA"
BLUE, RED, GREEN, PURPLE, ORANGE, CYAN = "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3"
FONT = "Arial, Helvetica, sans-serif"
svg = []; add = svg.append

def text(x, y, s, size=26, weight="normal", anchor="start", fill=INK):
    add(f"<text x='{x:.1f}' y='{y:.1f}' font-family='{FONT}' font-size='{size}' font-weight='{weight}' text-anchor='{anchor}' fill='{fill}'>{s}</text>")
def line(x1, y1, x2, y2, col=INK, w=3, dash=None, head=False):
    d = f" stroke-dasharray='{dash}'" if dash else ""; h = " marker-end='url(#ah)'" if head else ""
    add(f"<line x1='{x1:.1f}' y1='{y1:.1f}' x2='{x2:.1f}' y2='{y2:.1f}' stroke='{col}' stroke-width='{w}'{d}{h}/>")
def cell(x, y, col, r=14, y_mark=True, dead=False):
    """a cell; the small bar inside marks chromosome Y, absent when the cell has lost Y"""
    if dead:
        add(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r}' fill='none' stroke='{col}' stroke-width='2.5' stroke-dasharray='4 4'/>"); return
    add(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r}' fill='{col}' fill-opacity='0.22' stroke='{col}' stroke-width='3'/>")
    if y_mark: add(f"<rect x='{x - 2.5:.1f}' y='{y - r * 0.5:.1f}' width='5' height='{r:.1f}' fill='{col}'/>")
def axes(x0, y0, w, h, xlab=None, ylab=None):
    line(x0, y0, x0 + w, y0, INK, 2.5); line(x0, y0, x0, y0 - h, INK, 2.5)
    if xlab: text(x0 + w / 2, y0 + 30, xlab, 23, anchor="middle")
    if ylab: add(f"<text transform='translate({x0 - 14},{y0 - h / 2}) rotate(-90)' font-family='{FONT}' font-size='23' text-anchor='middle' fill='{INK}'>{ylab}</text>")
def bars(x0, y0, vals, cols, labs, vmax, h=120, bw=46, gap=26, fmt="{:.1f}"):
    for i, (v, c, l) in enumerate(zip(vals, cols, labs)):
        x = x0 + 14 + i * (bw + gap); hh = h * max(v, 0) / vmax
        add(f"<rect x='{x}' y='{y0 - hh:.1f}' width='{bw}' height='{hh:.1f}' fill='{c}'/>")
        text(x + bw / 2, y0 - hh - 8, fmt.format(v), 22, anchor="middle")
        text(x + bw / 2, y0 + 28, l, 21, anchor="middle")

add("<svg xmlns='http://www.w3.org/2000/svg' width='180mm' height='@HMM@mm' viewBox='0 0 1800 @H@'>")
add(f"<defs><marker id='ah' viewBox='0 0 10 10' refX='8' refY='5' markerWidth='4' markerHeight='4' orient='auto'><path d='M0,0 L10,5 L0,10 Z' fill='{GREY}'/></marker></defs>")
add("<rect width='1800' height='@H@' fill='white'/>")
C1, C2, C3 = 40, 560, 1250
text(C1, 58, "Data", 34, "bold"); text(C2, 58, "Comparison", 34, "bold"); text(C3, 58, "Readout", 34, "bold")
RH, Y0, FIGX = 260, 90, 1760
def row_frame(k, col, title, lines_, fig):
    y = Y0 + k * RH
    add(f"<rect x='{C1 - 20}' y='{y + 12}' width='8' height='{RH - 34}' fill='{col}'/>")
    text(C1, y + 44, title, 29, "bold", fill=col)
    for i, s in enumerate(lines_): text(C1, y + 84 + i * 33, s, 25)
    text(FIGX, y + 44, fig, 27, "bold", "end", col)
    if k: line(C1, y, W - 40, y, LINE, 2)
    return y

# ---- row 1: depth and detection
y = row_frame(0, BLUE, "Measuring loss in single cells", ["23 public atlases", "9.6 million cells, 2,235 men", "23 organs, 33 diseases", f"{N['male_groups_all']:,} donor-cell-type groups"], "Fig. 2")
cy = y + 120
for i, (r, nY) in enumerate([(22, 3), (14, 1), (9, 0)]):
    x = C2 + 60 + i * 120
    add(f"<circle cx='{x}' cy='{cy}' r='{r}' fill='{BLUE}' fill-opacity='0.22' stroke='{BLUE}' stroke-width='3'/>")
    add(f"<rect x='{x - 2.5}' y='{cy - r * 0.5:.1f}' width='5' height='{r:.1f}' fill='{BLUE}'/>")
    for j in range(nY): add(f"<circle cx='{x - 10 + j * 10}' cy='{cy + r + 18}' r='4' fill='{RED}'/>")
    text(x, cy + r + 50, ["deep", "medium", "shallow"][i], 22, anchor="middle")
text(C2 + 60, y + 64, "intact cells, Y reads (red)", 22, fill=INK)
line(C2 + 390, cy, C2 + 450, cy, GREY, 4, head=True)
text(C2 + 470, cy - 20, "zero reads", 22); text(C2 + 470, cy + 10, "or loss?", 22)
bx, by, bh = C3 + 40, y + 190, 115
axes(bx, by, 330, bh, ylab="called lost (%)")
pz = list(N["pseudo_zero_rule_by_bin"].values()); pm = list(N["pseudo_model_by_bin"].values())
px = [bx + 30 + i * 90 for i in range(4)]
for vals, col in ((pz, BLUE), (pm, RED)):
    add("<polyline points='" + " ".join(f"{x:.1f},{by - v / 60 * bh:.1f}" for x, v in zip(px, vals)) + f"' fill='none' stroke='{col}' stroke-width='4'/>")
    for x, v in zip(px, vals): add(f"<circle cx='{x:.1f}' cy='{by - v / 60 * bh:.1f}' r='6' fill='{col}'/>")
for x, l in zip(px, ["&lt;1", "1-2", "2-5", "&gt;5"]): text(x, by + 28, l, 21, anchor="middle")
text(bx + 165, by + 54, "expected Y reads per cell", 22, anchor="middle")
text(bx + 345, by - bh + 20, "zero-read rule", 22, fill=BLUE); text(bx + 345, by - 10, "detection model", 22, fill=RED)

# ---- row 2: physical loss
y = row_frame(1, RED, "Loss or silencing", ["PAR1 genes, X and Y copies", "AIDA, retina, microglia", "kidney nuclei with", "chromatin fragments"], "Fig. 3")
cx0, cy = C2 + 40, y + 110
def chrom(x, yy, lab, col, par=True):
    add(f"<rect x='{x}' y='{yy}' width='130' height='22' rx='11' fill='{col}' fill-opacity='0.3' stroke='{col}' stroke-width='2.5'/>")
    if par: add(f"<rect x='{x}' y='{yy}' width='26' height='22' rx='8' fill='{ORANGE}'/>")
    text(x + 145, yy + 18, lab, 22)
chrom(cx0, cy - 40, "X", GREY); chrom(cx0, cy - 8, "Y", BLUE)
text(cx0, cy + 42, "intact: two PAR1 copies", 22)
chrom(cx0 + 260, cy - 40, "X", GREY)
add(f"<rect x='{cx0 + 260}' y='{cy - 8}' width='130' height='22' rx='11' fill='none' stroke='{BLUE}' stroke-width='2.5' stroke-dasharray='5 5'/>")
text(cx0 + 260, cy + 42, "lost: one copy", 22)
text(cx0, cy + 76, "PAR1", 22, fill=ORANGE)
d = N["dn4"]
bars(C3 + 30, y + 200, [d["atac"], d["model"], d["naive"], d["trough"]], [GREEN, RED, BLUE, PURPLE], ["ATAC", "model", "zero", "trough"], 32)
text(C3 + 330, y + 120, "one donor,", 22); text(C3 + 330, y + 150, "cells without Y (%)", 22)

# ---- row 3: published claims and disease contrasts
y = row_frame(5, GREEN, "Published estimates", ["Europe PMC: 551 records,", "25 single-cell LOY studies", "OneK1K, kidney, microglia", "kidney and heart disease"], "Fig. 7")
cy = y + 120
for i in range(4): cell(C2 + 40 + i * 34, cy - 15, GREEN, 13)
for i in range(4): cell(C2 + 40 + i * 34, cy + 22, GREEN, 9)
text(C2 + 30, y + 64, "control", 22); text(C2 + 250, y + 64, "disease", 22)
for i in range(4): cell(C2 + 260 + i * 30, cy - 15, RED, 9)
for i in range(4): cell(C2 + 260 + i * 30, cy + 22, RED, 9)
text(C2 + 30, cy + 70, "deeper cells", 22); text(C2 + 250, cy + 70, "shallower cells", 22)
text(C2 + 30, cy + 102, "zero reads rise with lower depth", 22, fill=GREEN)
bx, by, bw = C3 + 170, y + 185, 260
xs_ = lambda v: bx + bw * (v + 5) / 15
line(xs_(0), by - 140, xs_(0), by, GREY, 2, "5 4"); line(bx, by, bx + bw, by, INK, 2.5)
for v in (-5, 0, 5, 10): text(xs_(v), by + 28, f"{v}", 21, anchor="middle")
for i, (lab, a, b) in enumerate([("kidney CKD", 6.6, 1.2), ("kidney AKI", -3.7, -2.1), ("heart DCM", 5.4, 0.6)]):
    yy = by - 115 + i * 38
    text(bx - 14, yy + 8, lab, 22, anchor="end")
    add(f"<circle cx='{xs_(a):.1f}' cy='{yy}' r='8' fill='{BLUE}'/>"); add(f"<circle cx='{xs_(b):.1f}' cy='{yy}' r='8' fill='{RED}'/>")
text(bx + bw + 14, by - 50, "zero-read", 22, fill=BLUE); text(bx + bw + 14, by - 20, "model", 22, fill=RED)
text(bx + bw / 2, by + 56, "disease minus control (points)", 22, anchor="middle")

# ---- row 4: where LOY occurs
y = row_frame(2, PURPLE, "Where loss occurs", ["1,526 male donors,", "10,219 groups", "four blood cohorts"], "Fig. 4")
cy = y + 125
# neuron, vessel, blood cell pictograms
nx = C2 + 60
add(f"<circle cx='{nx}' cy='{cy}' r='16' fill='{PURPLE}' fill-opacity='0.25' stroke='{PURPLE}' stroke-width='3'/>")
for a in (0.5, 2.0, 3.6, 5.0): line(nx + 16 * math.cos(a), cy + 16 * math.sin(a), nx + 44 * math.cos(a), cy + 44 * math.sin(a), PURPLE, 3)
text(nx, cy + 76, "neurons", 22, anchor="middle")
vx = C2 + 230
add(f"<path d='M{vx - 60},{cy - 22} Q{vx},{cy - 40} {vx + 60},{cy - 22}' fill='none' stroke='{RED}' stroke-width='5'/>")
add(f"<path d='M{vx - 60},{cy + 22} Q{vx},{cy + 4} {vx + 60},{cy + 22}' fill='none' stroke='{RED}' stroke-width='5'/>")
for i in range(3): cell(vx - 36 + i * 36, cy, RED, 9, y_mark=i != 1)
text(vx, cy + 76, "vessels", 22, anchor="middle")
bxx = C2 + 420
for i in range(5): cell(bxx - 50 + (i % 3) * 34, cy - 14 + (i // 3) * 32, BLUE, 12, y_mark=i not in (1, 3))
text(bxx - 15, cy + 76, "blood, microglia", 22, anchor="middle")
bars(C3 + 30, y + 205, [17.2, 9.6, 3.2, 0.1], [RED, RED, BLUE, PURPLE], ["VLMC", "endo", "micro", "neuron"], 20)
text(C3 + 300, y + 90, "brain, 468 donors,", 22); text(C3 + 300, y + 120, "cells without Y (%)", 22)

# ---- row 5: systemic or local
y = row_frame(3, ORANGE, "Systemic or local", ["Tabula Sapiens,", "same donor,", "several tissues"], "Fig. 5")
cy = y + 120
add(f"<circle cx='{C2 + 50}' cy='{cy - 30}' r='16' fill='none' stroke='{INK}' stroke-width='3'/>")
add(f"<path d='M{C2 + 50},{cy - 14} L{C2 + 50},{cy + 30} M{C2 + 26},{cy} L{C2 + 74},{cy} M{C2 + 50},{cy + 30} L{C2 + 34},{cy + 60} M{C2 + 50},{cy + 30} L{C2 + 66},{cy + 60}' stroke='{INK}' stroke-width='3' fill='none'/>")
for i, t in enumerate(["lung", "spleen", "skin"]):
    yy = cy - 50 + i * 48
    line(C2 + 90, cy, C2 + 170, yy, GREY, 2.5, head=True)
    text(C2 + 180, yy + 8, t, 22)
    cell(C2 + 280, yy, BLUE, 11, y_mark=False); cell(C2 + 320, yy, RED, 11, y_mark=(i != 1))
text(C2 + 260, cy + 100, "myeloid", 22, fill=BLUE); text(C2 + 360, cy + 100, "endothelial", 22, fill=RED)
bars(C3 + 30, y + 205, [79, 24, 16, 8], [BLUE, ORANGE, PURPLE, RED], ["myeloid", "stroma", "lymph", "endo"], 90, fmt="{:.0f}", gap=40)
text(C3 + 330, y + 90, "variance explained", 22); text(C3 + 330, y + 120, "by donor (%)", 22)

# ---- row 6: cells without Y
y = row_frame(4, CYAN, "Cells without Y", ["Immune Aging,", "memory regulatory T cells", "Tabula Sapiens B cells"], "Fig. 6")
cy = y + 120
cell(C2 + 70, cy, CYAN, 26); text(C2 + 70, cy + 64, "with Y", 22, anchor="middle")
cell(C2 + 250, cy, CYAN, 26, y_mark=False); text(C2 + 250, cy + 64, "without Y", 22, anchor="middle")
line(C2 + 110, cy, C2 + 200, cy, GREY, 4, head=True)
text(C2 + 320, cy - 10, "same donor and", 22); text(C2 + 320, cy + 20, "cell type, all genes", 22)
bx, by, bw = C3 + 140, y + 200, 260
xs_ = lambda v: bx + bw * (v + 1) / 2
line(xs_(0), by - 150, xs_(0), by, GREY, 2); line(bx, by, bx + bw, by, INK, 2.5)
for v in (-1, 0, 1): text(xs_(v), by + 28, f"{v}", 21, anchor="middle")
for i, (g, v) in enumerate([("TCL1A", 0.99), ("PRDM1", 0.22), ("IL7R", -0.23), ("SATB1", -0.65)]):
    yy = by - 125 + i * 32; c = RED if v > 0 else BLUE
    add(f"<rect x='{min(xs_(0), xs_(v)):.1f}' y='{yy - 10}' width='{abs(xs_(v) - xs_(0)):.1f}' height='20' fill='{c}'/>")
    text(bx - 12, yy + 8, g, 22, anchor="end")
text(bx + bw / 2, by + 56, "log2 ratio, without / with Y", 22, anchor="middle")

y_end = Y0 + 6 * RH
for x in (C2 - 30, C3 - 30): line(x, 30, x, y_end - 10, LINE, 3)
add("</svg>")
Hf = y_end + 10
out = "\n".join(svg).replace("@HMM@", f"{Hf / 10:.1f}").replace("@H@", str(Hf))
p = os.path.join(OUT, "flowchart.svg"); open(p, "w").write(out + "\n")
subprocess.run(["rsvg-convert", "-w", str(round(180 / 25.4 * 600)), "-b", "white", "-o", os.path.join(OUT, "flowchart.png"), p], check=True)
subprocess.run(["rsvg-convert", "-f", "pdf", "-o", os.path.join(OUT, "flowchart.pdf"), p], check=True)
print("wrote", p)
