# summarize loycell output (numbers from Rust): LOY vs Y+ log2FC by gene class, with pseudo-zero control
import csv, sys, statistics as st
PAR1 = ["PLCXD1","GTPBP6","PPP2R3B","SHOX","CRLF2","CSF2RA","IL3RA","SLC25A6","ASMTL","P2RY8","AKAP17A","ASMT","DHRSX","ZBED1","CD99"]
PAR2 = ["SPRY3","VAMP7","IL9R"]
GAM = {"DDX3X":"DDX3Y","KDM5C":"KDM5D","KDM6A":"UTY","EIF1AX":"EIF1AY","ZFX":"ZFY","RPS4X":"RPS4Y1","USP9X":"USP9Y","NLGN4X":"NLGN4Y","TBL1X":"TBL1Y","PRKX":"PRKY","TMSB4X":"TMSB4Y","TXLNG":"TXLNGY"}
Y = set(GAM.values())
for fn in sys.argv[1:]:
    R = [r for r in csv.DictReader(open(fn), delimiter="\t") if float(r["mean_ypos_cp10k"]) >= 0.5]
    if not R: continue
    print(f"== {fn}: groups {len({r['gid'] for r in R})}")
    def show(label, genes):
        g = [r for r in R if r["gene"] in genes]
        if not g: return
        l = [float(r["log2fc_loy"]) for r in g]; z = [float(r["log2fc_pseudo"]) for r in g if r["log2fc_pseudo"] != "nan"]
        print(f"  {label:26s} n {len(g):5d} | LOY log2FC median {st.median(l):+.2f} | pseudo-zero control {st.median(z):+.2f}")
    show("Y genes (MSY)", Y); show("PAR1 (X+Y copies)", PAR1); show("PAR2", PAR2)
    show("X gametologs", set(GAM)); show("XIST", {"XIST"})
    other = {r["gene"] for r in R} - Y - set(PAR1) - set(PAR2) - set(GAM) - {"XIST","TSIX","CD99","TCL1A"}
    show("other chrX genes", other)
    for gname in ["CD99","SLC25A6","AKAP17A","DDX3X","KDM5C","EIF1AX","ZFX","RPS4X","KDM6A","TCL1A"]:
        g = [r for r in R if r["gene"] == gname]
        if len(g) >= 5: print(f"    {gname:8s} groups {len(g):4d} LOY {st.median(float(r['log2fc_loy']) for r in g):+.2f}  pseudo {st.median(float(r['log2fc_pseudo']) for r in g if r['log2fc_pseudo']!='nan'):+.2f}")
