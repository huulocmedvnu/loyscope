# pool per-library ENCODE fits into one table with donor identifiers from the manifest
import csv, glob, os, sys
src, manifest, out = sys.argv[1:4]
man = {r[0]: r for r in csv.reader(open(manifest), delimiter="\t")}
hdr = None
with open(out, "w") as w:
    for f in sorted(glob.glob(os.path.join(src, "ENCSR*.fitamb.tsv"))):
        acc = os.path.basename(f).split(".")[0]
        rd = csv.DictReader(open(f), delimiter="\t")
        if hdr is None:
            hdr = rd.fieldnames; w.write("\t".join(hdr) + "\n")
        for r in rd:
            r["donor"] = man[acc][1]; r["assay"] = "10x multiome (ENCODE)"
            w.write("\t".join(r[k] for k in hdr) + "\n")
