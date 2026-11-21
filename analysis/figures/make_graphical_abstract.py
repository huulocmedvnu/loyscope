"""Graphical abstract as an editable SVG, exported to PNG and TIFF (600 dpi) and PDF.

Canvas units: 10 per mm, 180 x 90 mm, so 8 pt is about 28 units. Values are those reported in the Results
(results/tables/numbers.json and the brain, blood and Tabula Sapiens analyses).
"""
import json, math, os, subprocess
from common import TAB, OUT

N = json.load(open(os.path.join(TAB, "numbers.json")))
meta = {p[0]: p[1] for p in N["age_slope_meta"]}
W, H = 1800, 900
BLUE, RED, GREEN, PURPLE, ORANGE, INK, GREY = "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#111111", "#8C96A8"
FONT = "Arial, Helvetica, sans-serif"
svg = []; add = svg.append
def text(x, y, s, size=29, weight="normal", anchor="start", fill=INK, style=""):
    add(f"<text x='{x:.1f}' y='{y:.1f}' font-family='{FONT}' font-size='{size}' font-weight='{weight}' text-anchor='{anchor}' fill='{fill}' {style}>{s}</text>")
def line(x1, y1, x2, y2, col=INK, w=3, dash=None, head=False):
    d = f" stroke-dasharray='{dash}'" if dash else ""; h = " marker-end='url(#ah)'" if head else ""
    add(f"<line x1='{x1:.1f}' y1='{y1:.1f}' x2='{x2:.1f}' y2='{y2:.1f}' stroke='{col}' stroke-width='{w}'{d}{h}/>")
def cell(x, y, col, r=22, y_bar=True, reads=0):
    """cell with chromosome Y drawn as a bar; dashed outline of the bar when Y is lost; red dots are Y reads"""
    add(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r}' fill='{col}' fill-opacity='0.2' stroke='{col}' stroke-width='3'/>")
    if y_bar: add(f"<rect x='{x - 3.5:.1f}' y='{y - r * 0.55:.1f}' width='7' height='{r * 1.1:.1f}' rx='3' fill='{col}'/>")
    else: add(f"<rect x='{x - 3.5:.1f}' y='{y - r * 0.55:.1f}' width='7' height='{r * 1.1:.1f}' rx='3' fill='none' stroke='{col}' stroke-width='2' stroke-dasharray='3 3'/>")
    for j in range(reads): add(f"<circle cx='{x - (reads - 1) * 6 + j * 12:.1f}' cy='{y + r + 14:.1f}' r='4.5' fill='{RED}'/>")
def bar(x, y0, h, w, col, val, lab, hmax=230, vmax=20):
    hh = hmax * val / vmax
    add(f"<rect x='{x}' y='{y0 - hh:.1f}' width='{w}' height='{hh:.1f}' fill='{col}'/>")
    text(x + w / 2, y0 - hh - 10, f"{val:g}%" if val >= 1 else "0%", 27, anchor="middle")
    text(x + w / 2, y0 + 34, lab, 26, anchor="middle")

add(f"<svg xmlns='http://www.w3.org/2000/svg' width='180mm' height='90mm' viewBox='0 0 {W} {H}'>")
add(f"<defs><marker id='ah' viewBox='0 0 10 10' refX='8' refY='5' markerWidth='4' markerHeight='4' orient='auto'><path d='M0,0 L10,5 L0,10 Z' fill='{GREY}'/></marker></defs>")
add(f"<rect width='{W}' height='{H}' fill='white'/>")
for x in (600, 1200): line(x, 40, x, 790, "#DCE1EA", 3)

# ---------------- panel 1: measuring loss ----------------
x0 = 40
text(x0, 80, "Measuring loss of Y", 38, "bold")
text(x0, 122, "in single cells", 38, "bold")
# intact cells at three depths, and a cell without Y
for i, (r, k) in enumerate([(30, 3), (22, 1), (16, 0)]):
    cell(x0 + 60 + i * 110, 230, BLUE, r, True, k)
text(x0, 312, "intact cells, fewer Y reads (red)", 25)
text(x0, 342, "at lower sequencing depth", 25)
cell(x0 + 450, 230, RED, 26, False, 0)
text(x0 + 450, 312, "Y lost", 25, anchor="middle", fill=RED)
# control genes: zero-read rule vs model
pz = N["pseudo_zero_rule_by_bin"]; pm = N["pseudo_model_by_bin"]
z = list(pz.values())[0]; mz = list(pm.values())[0]
text(x0, 410, "Genes that cannot be lost,", 28, "bold")
text(x0, 446, "cells called lost", 28)
bar(x0 + 70, 655, 0, 110, BLUE, round(z), "zero reads", hmax=160, vmax=60)
bar(x0 + 300, 655, 0, 110, RED, round(mz), "model", hmax=160, vmax=60)
text(x0, 738, "PAR1 dosage halved and", 28, fill=RED)
text(x0, 772, "chromatin agrees in cells called lost", 28, fill=RED)

# ---------------- panel 2: where loss occurs ----------------
x0 = 640
text(x0, 80, "Where men lose Y", 38, "bold")
text(x0, 122, "brains of 468 men", 30)
# small icons above bars: vessel, microglia, neuron
def vessel(cx, cy):
    add(f"<path d='M{cx - 45},{cy - 16} Q{cx},{cy - 30} {cx + 45},{cy - 16}' fill='none' stroke='{RED}' stroke-width='5'/>")
    add(f"<path d='M{cx - 45},{cy + 16} Q{cx},{cy + 2} {cx + 45},{cy + 16}' fill='none' stroke='{RED}' stroke-width='5'/>")
def neuron(cx, cy):
    add(f"<circle cx='{cx}' cy='{cy}' r='13' fill='{GREEN}' fill-opacity='0.3' stroke='{GREEN}' stroke-width='3'/>")
    for a in (0.4, 1.9, 3.4, 4.9): line(cx + 13 * math.cos(a), cy + 13 * math.sin(a), cx + 36 * math.cos(a), cy + 36 * math.sin(a), GREEN, 3)
def microglia(cx, cy):
    add(f"<circle cx='{cx}' cy='{cy}' r='12' fill='{BLUE}' fill-opacity='0.3' stroke='{BLUE}' stroke-width='3'/>")
    for a in (0.2, 1.3, 2.4, 3.5, 4.6, 5.7): line(cx + 12 * math.cos(a), cy + 12 * math.sin(a), cx + 26 * math.cos(a) + 5, cy + 26 * math.sin(a), BLUE, 2)
base = 520
vessel(x0 + 70, 190); vessel(x0 + 200, 190); microglia(x0 + 330, 190); neuron(x0 + 460, 190)
bar(x0 + 25, base, 0, 90, RED, 17.2, "VLMC")
bar(x0 + 155, base, 0, 90, RED, 9.6, "endothelium")
bar(x0 + 285, base, 0, 90, BLUE, 3.2, "microglia")
bar(x0 + 415, base, 0, 90, GREEN, 0.1, "neurons")
text(x0, 620, "Rising about 3 points per decade", 28, fill=RED)
text(x0, 656, "in vascular cells, also in the retina", 28, fill=RED)
text(x0, 730, f"Blood: regulatory T cells fastest,", 28, fill=PURPLE)
text(x0, 766, f"+{meta['regulatory T cell']:.1f} points per decade, 4 cohorts", 28, fill=PURPLE)

# ---------------- panel 3: how loss spreads ----------------
x0 = 1240
text(x0, 80, "How loss spreads", 38, "bold")
text(x0, 122, "same man, several tissues", 30)
px, py = x0 + 70, 290
add(f"<circle cx='{px}' cy='{py - 70}' r='24' fill='none' stroke='{INK}' stroke-width='4'/>")
add(f"<path d='M{px},{py - 46} L{px},{py + 40} M{px - 36},{py - 10} L{px + 36},{py - 10} M{px},{py + 40} L{px - 26},{py + 96} M{px},{py + 40} L{px + 26},{py + 96}' stroke='{INK}' stroke-width='4' fill='none'/>")
tissues = ["lung", "spleen", "skin"]
for i, t in enumerate(tissues):
    yy = 190 + i * 95
    line(px + 50, py - 10, x0 + 210, yy, GREY, 3, head=True)
    text(x0 + 222, yy + 9, t, 27)
    cell(x0 + 345, yy, BLUE, 18, False)            # myeloid: shared loss
    cell(x0 + 425, yy, RED, 18, i != 1)            # endothelial: loss in one tissue only
text(x0 + 360, 440, "myeloid", 27, "bold", "end", BLUE); text(x0 + 410, 440, "endothelial", 27, "bold", "start", RED)
text(x0, 560, "Myeloid loss shared across tissues", 28, fill=BLUE)
text(x0, 596, "(79% of variance from the donor)", 28, fill=BLUE)
text(x0, 666, "Endothelial loss arises locally", 28, fill=RED)
text(x0, 702, "(8% from the donor)", 28, fill=RED)
text(x0, 766, "Large clones in a minority of men", 28)

# ---------------- footer ----------------
text(W / 2, 862, "Loss of chromosome Y is confined to dividing cells and concentrated in the vasculature", 31, "normal", "middle", INK, "font-style='italic'")
add("</svg>")
p = os.path.join(OUT, "graphical_abstract.svg"); open(p, "w").write("\n".join(svg) + "\n")
w600 = round(180 / 25.4 * 600)
png = os.path.join(OUT, "graphical_abstract.png")
subprocess.run(["rsvg-convert", "-w", str(w600), "-b", "white", "-o", png, p], check=True)
subprocess.run(["rsvg-convert", "-f", "pdf", "-o", os.path.join(OUT, "graphical_abstract.pdf"), p], check=True)
try:
    from PIL import Image
    Image.open(png).convert("RGB").save(os.path.join(OUT, "graphical_abstract.tif"), dpi=(600, 600), compression="tiff_lzw")
except ImportError:
    pass
print("wrote", p)
