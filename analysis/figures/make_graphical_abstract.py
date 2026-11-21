"""Graphical abstract as an editable SVG, exported to PNG and TIFF (600 dpi) and PDF.

Canvas units: 10 per mm, 180 x 90 mm, so 8 pt is about 28 units. Values are those reported in the Results
(results/tables/numbers.json and the brain, blood and Tabula Sapiens analyses).
"""
import json, math, os, subprocess
from common import TAB, OUT

N = json.load(open(os.path.join(TAB, "numbers.json")))
import csv
meta = {r["cell_class"]: float(r["slope_per_decade"]) for r in csv.DictReader(open(os.path.join(os.environ.get("LOY_RESULTS", "results"), "tables", "robust_meta.tsv")), delimiter="\t")}
zero_rule = list(N["pseudo_zero_rule_by_bin"].values())[0]
W, H = 1800, 900
INK, MUTE, RULE = "#1B1F2A", "#5B6475", "#E3E7EE"
VAS, VAS_L = "#E4573D", "#FBE3DD"        # vasculature
MYE, MYE_L = "#4F5BD5", "#E2E5FB"        # myeloid, blood
NEU, NEU_L = "#14A37F", "#D7F2EA"        # neurons and glia
TRE = "#9B51E0"                           # regulatory T cells
XCH, YCH, PAR = "#9AA3B2", "#4F5BD5", "#F5A623"
FONT = "Arial, Helvetica, sans-serif"
svg = []; add = svg.append

def text(x, y, s, size=27, weight="normal", anchor="start", fill=INK, extra=""):
    add(f"<text x='{x:.1f}' y='{y:.1f}' font-family='{FONT}' font-size='{size}' font-weight='{weight}' text-anchor='{anchor}' fill='{fill}' {extra}>{s}</text>")
def path(d, fill="none", stroke=INK, w=3, extra=""):
    add(f"<path d='{d}' fill='{fill}' stroke='{stroke}' stroke-width='{w}' stroke-linecap='round' stroke-linejoin='round' {extra}/>")
def circ(x, y, r, fill, stroke, w=3, extra=""):
    add(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r:.1f}' fill='{fill}' stroke='{stroke}' stroke-width='{w}' {extra}/>")
def chrom_x(x, y, s=1.0, col=XCH):
    for a in (-24, 24):
        add(f"<g transform='translate({x:.1f},{y:.1f}) rotate({a}) scale({s:.3f})'><rect x='-6' y='-30' width='12' height='60' rx='6' fill='{col}'/><rect x='-6' y='-30' width='12' height='9' rx='4' fill='{PAR}'/></g>")
def chrom_y(x, y, s=1.0, col=YCH, lost=False):
    if lost:
        add(f"<g transform='translate({x:.1f},{y:.1f}) scale({s:.3f})'><rect x='-6' y='-20' width='12' height='40' rx='6' fill='none' stroke='{col}' stroke-width='2.5' stroke-dasharray='4 4'/></g>")
    else:
        add(f"<g transform='translate({x:.1f},{y:.1f}) scale({s:.3f})'><rect x='-6' y='-20' width='12' height='40' rx='6' fill='{col}'/><rect x='-6' y='-20' width='12' height='8' rx='4' fill='{PAR}'/></g>")
def nucleus(x, y, r, col_fill, col_line, y_lost=False):
    circ(x, y, r, col_fill, col_line, 3)
    chrom_x(x - r * 0.28, y, r / 60); chrom_y(x + r * 0.32, y + r * 0.08, r / 60, lost=y_lost)

add(f"<svg xmlns='http://www.w3.org/2000/svg' width='180mm' height='90mm' viewBox='0 0 {W} {H}'>")
add(f"<defs><marker id='ah' viewBox='0 0 10 10' refX='8' refY='5' markerWidth='5' markerHeight='5' orient='auto'><path d='M0,0 L10,5 L0,10 Z' fill='{MYE}'/></marker></defs>")
add(f"<rect width='{W}' height='{H}' fill='white'/>")
for x in (585, 1215): add(f"<line x1='{x}' y1='50' x2='{x}' y2='780' stroke='{RULE}' stroke-width='3'/>")

# ================= panel 1: measuring =================
x0 = 45
text(x0, 88, "A reliable count of", 36, "bold"); text(x0, 130, "cells without Y", 36, "bold")
nucleus(x0 + 90, 255, 58, "#F4F6FA", "#C9CFDA")
text(x0 + 90, 352, "Y present,", 25, anchor="middle", fill=MUTE); text(x0 + 90, 381, "no read captured", 25, anchor="middle", fill=MUTE)
nucleus(x0 + 330, 255, 58, "#F4F6FA", "#C9CFDA", y_lost=True)
text(x0 + 330, 352, "Y lost", 25, anchor="middle", fill=MUTE); text(x0 + 330, 381, "from the cell", 25, anchor="middle", fill=MUTE)
text(x0 + 210, 264, "or", 28, "bold", "middle", MUTE)
text(x0, 440, "Genes that cannot be lost", 26, "bold")
bx, by = x0 + 20, 640
for i, (v, col, lab) in enumerate([(zero_rule, "#B8BED0", "zero-read rule"), (0.0, VAS, "this study")]):
    h = max(135 * v / 60, 3)
    add(f"<rect x='{bx + i * 230}' y='{by - h:.1f}' width='130' height='{h:.1f}' rx='6' fill='{col}'/>")
    text(bx + i * 230 + 65, by - h - 14, f"{v:.0f}%", 34, "bold", "middle", INK)
    text(bx + i * 230 + 65, by + 34, lab, 24, anchor="middle", fill=MUTE)
text(x0, 728, "Cells called lost carry half", 25); text(x0, 760, "the PAR1 dosage of intact cells", 25)

# ================= panel 2: where =================
x0 = 625
text(x0, 88, "Loss concentrates in", 36, "bold"); text(x0, 130, "vessels, not neurons", 36, "bold")
bxc, byc = x0 + 165, 305
path(f"M{bxc - 125},{byc + 20} C{bxc - 140},{byc - 60} {bxc - 70},{byc - 125} {bxc + 5},{byc - 115} C{bxc + 85},{byc - 135} {bxc + 150},{byc - 70} {bxc + 140},{byc} "
     f"C{bxc + 150},{byc + 55} {bxc + 95},{byc + 80} {bxc + 45},{byc + 70} C{bxc + 20},{byc + 105} {bxc - 30},{byc + 100} {bxc - 45},{byc + 72} C{bxc - 95},{byc + 85} {bxc - 120},{byc + 60} {bxc - 125},{byc + 20} Z",
     fill="#F7F4F2", stroke="#C9B9B2", w=3)
for d in (f"M{bxc - 60},{byc - 60} C{bxc - 40},{byc - 20} {bxc - 70},{byc + 10} {bxc - 50},{byc + 50}", f"M{bxc + 20},{byc - 80} C{bxc + 10},{byc - 40} {bxc + 40},{byc - 10} {bxc + 25},{byc + 40}",
          f"M{bxc + 90},{byc - 40} C{bxc + 70},{byc - 10} {bxc + 95},{byc + 15} {bxc + 80},{byc + 45}"):
    path(d, stroke="#D9CCC6", w=2.5)
path(f"M{bxc - 150},{byc + 95} C{bxc - 60},{byc + 40} {bxc - 10},{byc + 10} {bxc + 60},{byc - 30} C{bxc + 95},{byc - 50} {bxc + 120},{byc - 85} {bxc + 150},{byc - 110}", stroke=VAS, w=10)
path(f"M{bxc - 10},{byc + 12} C{bxc - 5},{byc + 50} {bxc + 20},{byc + 70} {bxc + 30},{byc + 95}", stroke=VAS, w=6)
for (cx, cy, lost) in ((bxc - 80, byc + 52, False), (bxc + 20, byc - 8, True), (bxc + 105, byc - 66, False)):
    circ(cx, cy, 10, "white", VAS, 3)
    if lost: path(f"M{cx - 4},{cy - 4} L{cx + 4},{cy + 4} M{cx + 4},{cy - 4} L{cx - 4},{cy + 4}", stroke=VAS, w=2.5)
nx, ny = bxc - 30, byc - 58
circ(nx, ny, 10, NEU_L, NEU, 3)
for a in (0.3, 1.8, 3.3, 4.6): path(f"M{nx + 10 * math.cos(a):.1f},{ny + 10 * math.sin(a):.1f} L{nx + 28 * math.cos(a):.1f},{ny + 28 * math.sin(a):.1f}", stroke=NEU, w=2.5)
vx = x0 + 360
for i, (lab, v, col) in enumerate([("leptomeningeal", "17%", VAS), ("endothelial", "10%", VAS), ("microglia", "3%", MYE), ("neurons, glia", "0%", NEU)]):
    yy = 222 + i * 60
    text(vx, yy, v, 34, "bold", fill=col); text(vx + 76, yy - 4, lab, 22)
text(vx, 470, "brains of 468 men", 22, fill=MUTE)
ex, ey = x0 + 40, 560
path(f"M{ex - 34},{ey} C{ex - 14},{ey - 28} {ex + 14},{ey - 28} {ex + 34},{ey} C{ex + 14},{ey + 28} {ex - 14},{ey + 28} {ex - 34},{ey} Z", fill="white", stroke=INK, w=3)
circ(ex, ey, 11, VAS_L, VAS, 3)
text(x0 + 95, 552, "Retinal arterioles and pericytes", 25); text(x0 + 95, 584, "lose Y with age", 25)
dx, dy = x0 + 40, 690
path(f"M{dx},{dy - 34} C{dx + 20},{dy - 8} {dx + 26},{dy + 6} {dx + 26},{dy + 14} C{dx + 26},{dy + 30} {dx + 14},{dy + 40} {dx},{dy + 40} C{dx - 14},{dy + 40} {dx - 26},{dy + 30} {dx - 26},{dy + 14} C{dx - 26},{dy + 6} {dx - 20},{dy - 8} {dx},{dy - 34} Z", fill="#F6E4FF", stroke=TRE, w=3)
text(x0 + 95, 686, "Blood: regulatory T cells lose Y", 25)
text(x0 + 95, 718, f"in four cohorts, +{meta['regulatory T cell']:.1f} points per decade", 25)

# ================= panel 3: how =================
x0 = 1255
text(x0, 88, "Two routes of spread", 36, "bold")
hx, hy = x0 + 110, 200
circ(hx, hy, 30, "#F4F6FA", "#C9CFDA", 3)
path(f"M{hx - 55},{hy + 52} C{hx - 55},{hy + 38} {hx - 40},{hy + 34} {hx},{hy + 34} C{hx + 40},{hy + 34} {hx + 55},{hy + 38} {hx + 55},{hy + 52} L{hx + 50},{hy + 220} L{hx + 22},{hy + 220} L{hx + 18},{hy + 360} L{hx - 18},{hy + 360} L{hx - 22},{hy + 220} L{hx - 50},{hy + 220} Z", fill="#F4F6FA", stroke="#C9CFDA", w=3)
mx, my = hx, hy + 205
circ(mx, my, 13, MYE_L, MYE, 3)
for ox, oy in ((hx - 22, hy + 85), (hx + 20, hy + 112), (hx + 28, hy + 162)):
    path(f"M{mx},{my - 13} Q{(mx + ox) / 2 - 34},{(my + oy) / 2} {ox},{oy + 11}", stroke=MYE, w=2.5, extra="stroke-dasharray='5 5' marker-end='url(#ah)'")
    circ(ox, oy, 9, MYE_L, MYE, 2.5)
text(hx - 22, my + 8, "marrow", 22, anchor="end", fill=MYE)
rx = x0 + 250
circ(rx + 12, 255, 12, MYE_L, MYE, 3)
text(rx + 40, 250, "Myeloid loss", 27, "bold", fill=MYE); text(rx + 40, 282, "shared by all organs", 25)
text(rx + 40, 314, "of a man (79%)", 25, fill=MUTE)
path(f"M{rx - 6},430 L{rx + 30},430", stroke=VAS, w=8)
text(rx + 40, 425, "Endothelial loss", 27, "bold", fill=VAS); text(rx + 40, 457, "arises in each", 25)
text(rx + 40, 489, "vascular bed (8%)", 25, fill=MUTE)
text(rx - 10, 640, "Large clones in a", 25); text(rx - 10, 672, "minority of older men", 25)

# ================= footer =================
add(f"<rect x='0' y='808' width='{W}' height='92' fill='#F5F7FB'/>")
text(W / 2, 850, "Loss of chromosome Y is confined to dividing cells and concentrated in the vasculature", 30, "bold", "middle", INK)
text(W / 2, 886, "10.5 million cells   •   2,287 men   •   24 organs   •   242 cell types   •   33 diseases   •   24 datasets", 22, "normal", "middle", MUTE)
add("</svg>")
p = os.path.join(OUT, "graphical_abstract.svg"); open(p, "w").write("\n".join(svg) + "\n")
w600 = round(180 / 25.4 * 600); png = os.path.join(OUT, "graphical_abstract.png")
subprocess.run(["rsvg-convert", "-w", str(w600), "-b", "white", "-o", png, p], check=True)
subprocess.run(["rsvg-convert", "-f", "pdf", "-o", os.path.join(OUT, "graphical_abstract.pdf"), p], check=True)
try:
    from PIL import Image
    Image.open(png).convert("RGB").save(os.path.join(OUT, "graphical_abstract.tif"), dpi=(600, 600), compression="tiff_lzw")
except ImportError:
    pass
print("wrote", p)
