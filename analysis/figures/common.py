# shared style and readers for manuscript figures (inputs are files written by the Rust tools)
import csv, os, re, glob, statistics as st
import plotly.io as pio
pio.kaleido.scope.mathjax = None  # otherwise kaleido prints a MathJax loading message into PDF exports
RES = os.environ.get("LOY_RESULTS", "results"); EXT = os.path.join(RES, "extract"); TAB = os.path.join(RES, "tables")
OUT = os.path.join(RES, "figures"); os.makedirs(OUT, exist_ok=True)
C = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692", "#B6E880", "#FECB52"]
BLACK = "#111111"; GREY = "#999999"
DEC = {"third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75, "ninth": 85}
def age(s):
    m = re.match(r"(\d+)-year-old", s)
    if m: return float(m.group(1))
    m = re.match(r"(\d+)-(\d+) year-old", s)
    if m: return (float(m.group(1)) + float(m.group(2))) / 2
    m = re.match(r"(\w+) decade", s); return DEC.get(m.group(1)) if m else None
def rows(fn): return list(csv.DictReader(open(fn), delimiter="\t"))
def ext(name): return os.path.join(EXT, name)
def style(fig, width=900, height=760, legends=None):
    fig.update_layout(template="plotly_white", width=width, height=height, font=dict(family="Arial", size=13, color=BLACK),
                      margin=dict(l=75, r=20, t=50, b=60), **(legends or {}))
    fig.update_xaxes(showgrid=False, showline=True, linecolor=BLACK, ticks="outside", tickfont=dict(color=BLACK), title_font=dict(color=BLACK))
    fig.update_yaxes(showgrid=False, showline=True, linecolor=BLACK, ticks="outside", tickfont=dict(color=BLACK), title_font=dict(color=BLACK))
    for ann in fig.layout.annotations:
        if ann.text and len(ann.text) > 2 and ann.text[1:3] == "  ": ann.update(xanchor="left", font=dict(size=14, color=BLACK))
    return fig
def panel_titles_left(fig, xs):
    i = 0
    for ann in fig.layout.annotations:
        if ann.text and len(ann.text) > 2 and ann.text[1:3] == "  ": ann.update(x=xs[i], xanchor="left"); i += 1
def save(fig, name):
    """PNG at 600 dpi for a 180 mm print width, vector PDF, and HTML"""
    scale = (180 / 25.4 * 600) / fig.layout.width
    fig.write_image(os.path.join(OUT, name + ".png"), scale=scale)
    fig.write_image(os.path.join(OUT, name + ".pdf"))
    fig.write_html(os.path.join(OUT, name + ".html")); print("written", name, f"{round(fig.layout.width * scale)} px wide")
def leg(x, y, xa="left", ya="top"): return dict(x=x, y=y, xanchor=xa, yanchor=ya, bgcolor="rgba(255,255,255,0)")
