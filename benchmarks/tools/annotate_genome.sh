#!/usr/bin/env bash
# Quick gene annotation for TSS classification when no curated annotation is
# available: Prodigal (protein-coding), barrnap (rRNA) and ARAGORN (tRNA).
# Writes gene + child (CDS/rRNA/tRNA) features, the layout ANNOgesic expects.
# Prefer a curated (NCBI) or Bakta annotation when you have one.
#
# usage: annotate_genome.sh genome.fa SEQID out.gff
# needs: prodigal, barrnap, aragorn, python3
set -euo pipefail
fa=$1; seqid=$2; out=$3
tmp=$(mktemp -d)
prodigal -q -i "$fa" -f gff -o "$tmp/prodigal.gff" -g 11
barrnap --quiet --kingdom bac "$fa" > "$tmp/barrnap.gff" 2> /dev/null
aragorn -t -gc11 -l -w "$fa" > "$tmp/aragorn.txt"
python3 - "$tmp" "$seqid" "$out" <<'EOF'
import re, sys
tmp, seqid, out_path = sys.argv[1:4]
features = []
for line in open(tmp + "/prodigal.gff"):
    if not line.startswith("#"):
        f = line.split("\t")
        features.append(("CDS", int(f[3]), int(f[4]), f[6], "cds"))
for line in open(tmp + "/barrnap.gff"):
    if not line.startswith("#"):
        f = line.split("\t")
        features.append(("rRNA", int(f[3]), int(f[4]), f[6], re.search(r"Name=([^;]+)", f[8]).group(1)))
for line in open(tmp + "/aragorn.txt"):
    m = re.match(r"\s*\d+\s+(tRNA-\S+)\s+(c?)\[(\d+),(\d+)\]", line)
    if m:
        features.append(("tRNA", int(m.group(3)), int(m.group(4)), "-" if m.group(2) else "+", m.group(1)))
features.sort(key=lambda x: x[1])
with open(out_path, "w") as out:
    out.write("##gff-version 3\n")
    for i, (kind, start, end, strand, name) in enumerate(features, 1):
        tag = "%s_%05d" % (seqid.split(".")[0], i)
        out.write("\t".join([seqid, "predicted", "gene", str(start), str(end), ".", strand, ".",
                             "ID=gene-%s;Name=%s;locus_tag=%s" % (tag, tag, tag)]) + "\n")
        out.write("\t".join([seqid, "predicted", kind, str(start), str(end), ".", strand,
                             "0" if kind == "CDS" else ".",
                             "ID=%s-%s;Parent=gene-%s;locus_tag=%s;product=%s" % (kind.lower(), tag, tag, tag, name)]) + "\n")
counts = {k: sum(1 for f in features if f[0] == k) for k in ("CDS", "rRNA", "tRNA")}
print("annotated", counts)
EOF
rm -rf "$tmp"
