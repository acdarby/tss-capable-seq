#!/usr/bin/env python3
"""Compare TSS called in mixed samples with each organism called alone (run.sh).

For each mixture (E. coli at a given share of its reads, mapping-quality
filter q20 or q0) and organism: enriched reads counted on the genome, TSS and
high-confidence TSS in the mixture and alone, the share of the TSS called
alone that the mixture also calls (+-2 nt, "recovered") and the share of the
mixture's TSS that were also called alone ("shared"; the rest are new in
the mixture, e.g. from reads of the other organism).

usage: compare.py CALLS_DIR > mixture.tsv
"""
import bisect
import collections
import csv
import glob
import os
import sys


def read_calls(path, organism=None):
    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if organism is not None and rows and "organism" in rows[0]:
        rows = [r for r in rows if r["organism"] == organism]
    return [(r["seqid"], r["strand"], int(r["position"]), r["confidence"]) for r in rows]


def reads_on(qc_path):
    with open(qc_path) as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["role"] == "enriched"]
    return sum(float(r["reads"]) for r in rows)


def share_found(query, reference, window=2):
    index = collections.defaultdict(list)
    for seqid, strand, pos, _ in reference:
        index[(seqid, strand)].append(pos)
    for v in index.values():
        v.sort()
    hit = 0
    for seqid, strand, pos, _ in query:
        v = index.get((seqid, strand), [])
        i = bisect.bisect_left(v, pos - window)
        hit += i < len(v) and v[i] <= pos + window
    return hit / len(query) if query else float("nan")


def main(calls_dir):
    fields = ["mixture", "mapq", "organism", "reads_mixed", "reads_alone", "tss_mixed",
              "tss_alone", "recovered", "shared", "high_mixed", "high_alone",
              "high_recovered", "high_shared"]
    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(fields)
    for path in sorted(glob.glob(os.path.join(calls_dir, "mix_*_q*.tsv"))):
        name = os.path.basename(path)[:-4]                 # mix_ecoli0.1_q20
        if name.endswith(".qc"):
            continue
        tag, mapq = name[4:].rsplit("_q", 1)
        for organism, alone in (("Ecoli", "alone_%s" % tag), ("Salmonella", "alone_salmonella")):
            mixed = read_calls(path, organism)
            ref = read_calls(os.path.join(calls_dir, alone + ".tsv"))
            high_m = [c for c in mixed if c[3] == "high"]
            high_r = [c for c in ref if c[3] == "high"]
            out.writerow([tag, mapq, organism,
                          "%.0f" % reads_on(os.path.join(calls_dir, "%s.%s.qc.tsv" % (name, organism))),
                          "%.0f" % reads_on(os.path.join(calls_dir, alone + ".qc.tsv")),
                          len(mixed), len(ref),
                          "%.3f" % share_found(ref, mixed), "%.3f" % share_found(mixed, ref),
                          len(high_m), len(high_r),
                          "%.3f" % share_found(high_r, mixed), "%.3f" % share_found(high_m, ref)])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
