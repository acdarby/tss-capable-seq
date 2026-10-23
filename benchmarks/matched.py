#!/usr/bin/env python3
"""Matched-size comparison of TSS rankings (docs/benchmark_plan.md, section 5).

Fixed cut-offs and FDR calls give different numbers of TSS, so their promoter
rates cannot be compared directly. Here every candidate position (at least
`--min-reads` enriched reads) is ranked by each method's own score, and the
top N sites (after suppressing lower-ranked sites within 6 nt on the same
strand) are scored at several N on evidence that no method uses:

* box10 / box35: -10 (TATAAT) and -35 (TTGACA) boxes with <= 1 mismatch,
  as in harness.py; box10_background is the same 100 nt downstream;
* purine_1: share with A or G at +1;
* in_rna: share within 5 nt of an rRNA or tRNA gene (processed ends);
* ref_recall_5: recall of the reference TSS within 5 nt (E. coli Table S1;
  that table was made with the Ettwiller pipeline, so it favours the RPM
  ranking).

Rankings: `rpm` (enriched RPM among positions with log2 enrichment >= 0, the
Ettwiller filter), `fold` (log2 enrichment with one pseudo-read, ties by RPM)
and `pvalue_d<D>` (beta-binomial p-value at dispersion D, ties by RPM).

usage: matched.py --data-dir DIR --dataset ecoli --out results/matched_ecoli
"""
import argparse
import math
import os
import sys
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
import cappable_tss as ct  # noqa: E402
import harness  # noqa: E402


def candidates(enr, ctl, exclude, min_reads):
    """Per-position counts for one enriched and one control library."""
    ne = ct.library_size(enr, exclude) or 1.0
    nc = ct.library_size(ctl, exclude) or 1.0
    rows = []
    for (seqid, strand), positions in enr.items():
        cpos = ctl.get((seqid, strand), {})
        for pos, e in positions.items():
            if e >= min_reads:
                rows.append((seqid, strand, pos, e, cpos.get(pos, 0.0)))
    return rows, ne, nc


def top_sites(ranked, n, distance=6):
    """First n sites in rank order, skipping any within `distance` nt of one kept."""
    kept, taken = [], set()
    for seqid, strand, pos in ranked:
        if any((seqid, strand, pos + d) in taken for d in range(-distance, distance + 1)):
            continue
        taken.add((seqid, strand, pos))
        kept.append((seqid, strand, pos))
        if len(kept) == n:
            break
    return kept


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--dataset", default="ecoli")
    parser.add_argument("--manifest", default=os.path.join(HERE, "datasets.toml"))
    parser.add_argument("--cache-dir")
    parser.add_argument("--min-reads", type=int, default=3)
    parser.add_argument("--dispersions", default="0,0.01,0.02,0.05")
    parser.add_argument("--sizes", default="1000,2000,4000,7000,11000,16000")
    parser.add_argument("--out", required=True, help="Output prefix (.tsv and .md)")
    args = parser.parse_args(argv)

    with open(args.manifest, "rb") as handle:
        ds = [d for d in tomllib.load(handle)["dataset"] if d["name"] == args.dataset][0]
    if not ds.get("control"):
        parser.error("dataset %s has no control library" % args.dataset)
    path = lambda p: os.path.join(args.data_dir, p)  # noqa: E731
    cache = args.cache_dir or args.out + "_cache"
    mate, mapq, chrom = ds.get("mate", 1), ds.get("min_mapq", 0), ds["chromosome"]
    enr = harness.load_signal(path(ds["enriched"][0]), mate, mapq, cache)
    ctl = harness.load_signal(path(ds["control"][0]), mate, mapq, cache)
    exclude = ct.GeneIndex(ct.read_genes(path(ds["annotation"]), ("rRNA", "tRNA")))
    seq = harness.read_fasta(path(ds["genome"]))[chrom]
    reference = None
    if ds.get("reference"):
        reference = [(st, p) for s, st, p in ct.read_sites(path(ds["reference"]))
                     if s in (chrom, "gi|48994873|gb|U00096.2|")]

    rows, ne, nc = candidates(enr, ctl, exclude, args.min_reads)
    rows = [r for r in rows if r[0] == chrom]
    print("%d candidate positions" % len(rows), file=sys.stderr)
    rpm = {r[:3]: r[3] / ne * 1e6 for r in rows}
    fold = {r[:3]: math.log2((r[3] / ne + 1 / ne) / (r[4] / nc + 1 / nc)) for r in rows}

    rankings = {}
    # Ettwiller filter: log2(E_rpm / C_rpm) >= 0 with a missing control counted as 1 read.
    passed = [r[:3] for r in rows if r[3] / ne >= max(r[4], 1.0) / nc]
    rankings["rpm"] = sorted(passed, key=lambda k: -rpm[k])
    rankings["fold"] = sorted(rpm, key=lambda k: (-fold[k], -rpm[k]))
    for d in [float(x) for x in args.dispersions.split(",")]:
        p0 = ne / (ne + nc)
        pval = {r[:3]: ct.betabinom_sf(int(r[3]), int(r[3] + r[4]), p0, d) for r in rows}
        rankings["pvalue_d%g" % d] = sorted(pval, key=lambda k: (pval[k], -rpm[k]))
        print("ranked by p-value at dispersion %g" % d, file=sys.stderr)

    out = []
    for n in [int(x) for x in args.sizes.split(",")]:
        for name, ranked in rankings.items():
            kept = top_sites(ranked, n)
            sites = [(st, p) for _, st, p in kept]
            m10, m35 = harness.box_rates(sites, seq)
            row = dict(dataset="%s_top%d" % (args.dataset, n), method=name, n_tss=len(sites),
                       box10=m10, box35=m35,
                       box10_background=harness.box_rates(harness.shifted(sites), seq)[0],
                       purine_1=sum(1 for st, p in sites
                                    if harness.upstream(seq, st, p, 0, 0) in ("A", "G"))
                       / len(sites) if sites else float("nan"),
                       in_rna=sum(1 for s, st, p in kept
                                  if exclude.overlapping(s, p - 5, p + 5)) / len(kept)
                       if kept else float("nan"))
            if reference:
                row["ref_recall_5"] = harness.match(reference, sites, 5) / len(reference)
            out.append(row)
            print("  top %-6d %-14s box10 %.3f" % (n, name, m10), file=sys.stderr)
    harness.write_rows(out, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
