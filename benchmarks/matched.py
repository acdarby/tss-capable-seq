#!/usr/bin/env python3
"""Matched-size comparison of TSS rankings (docs/benchmark_plan.md, section 5).

Methods call different numbers of TSS, so their promoter rates cannot be
compared directly. Here all methods start from the same candidates: pooled
5'-end count >= `--min-reads`, reads in every enriched replicate and, with a
control, not depleted after enrichment (RPM at least that of the control, a
missing control counting as one read, as in Ettwiller et al. 2016). The
candidates are ranked by each method's score, sites within 6 nt of a
higher-ranked site on the same strand are skipped, and the top N are scored
on evidence that no method uses:

* core10: sigma70 -10 core (T, A, T at -12, -11, -7, hexamer starting at
  -14..-10; harness.promoter_elements); core10_bg: composition-matched random
  positions (harness.MatchedBackground); ext10 and box35_spaced: extended -10
  and a -35 at 16-18 nt from the -10;
* box10 / box35: the earlier measure (TATAAT, TTGACA with <= 1 mismatch
  anywhere in -15..-5 and -40..-28); box10_background is 100 nt downstream;
* purine_1: share with A or G at +1;
* in_rna: share within 5 nt of an rRNA or tRNA gene;
* ref_recall_5: recall of reference TSS within 5 nt (E. coli Table S1, made
  with the Ettwiller pipeline, so it favours the RPM ranking).

The band table gives the core10 rate of ranks in consecutive bands, which
shows where calls stop looking like promoters, and marks where each method's
own threshold stops (the Ettwiller cut-off, and q <= 0.05 for the tests).

Rankings: `rpm` (pooled count), `local` (local-background p-value,
cappable_tss.call_local, dispersion estimated; `local_d<D>` at fixed
dispersions), `ratio` (count / local background rate) and, with a control and one replicate,
`control_d<D>` (beta-binomial p-value against the control at dispersion D).

usage: matched.py --data-dir DIR --dataset walbb --out results/matched_walbb
"""
import argparse
import os
import sys
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
import cappable_tss as ct  # noqa: E402
import harness  # noqa: E402


def top_sites(ranked, n=None, distance=6):
    """Sites in rank order, skipping any within `distance` nt of one kept."""
    kept, taken = [], set()
    for seqid, strand, pos in ranked:
        if any((seqid, strand, pos + d) in taken for d in range(-distance, distance + 1)):
            continue
        taken.add((seqid, strand, pos))
        kept.append((seqid, strand, pos))
        if n is not None and len(kept) == n:
            break
    return kept


def score(kept, seq, exclude, reference, background=None):
    sites = [(st, p) for _, st, p in kept]
    m10, m35 = harness.box_rates(sites, seq)
    prom = harness.promoter_rates(sites, seq, background)
    row = dict(n_tss=len(sites), core10=prom["core10"], core10_bg=prom.get("core10_bg"),
               box35_spaced=prom["box35"], ext10=prom["ext10"], box10=m10, box35=m35,
               box10_background=harness.box_rates(harness.shifted(sites), seq)[0],
               purine_1=sum(1 for st, p in sites if harness.upstream(seq, st, p, 0, 0) in "AG")
               / len(sites) if sites else float("nan"),
               in_rna=sum(1 for s, st, p in kept if exclude.overlapping(s, p - 5, p + 5))
               / len(kept) if kept else float("nan"))
    if reference:
        row["ref_recall_5"] = harness.match(reference, sites, 5) / len(reference)
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--dataset", default="ecoli")
    parser.add_argument("--manifest", default=os.path.join(HERE, "datasets.toml"))
    parser.add_argument("--cache-dir")
    parser.add_argument("--min-reads", type=int, default=3)
    parser.add_argument("--with-technical-replicate", action="store_true",
                        help="Add the dataset's technical_replicate as a second library")
    parser.add_argument("--dispersions", default="0,0.05",
                        help="Control-test dispersions (one replicate with a control only)")
    parser.add_argument("--local-dispersions", default="",
                        help="Extra local rankings at fixed dispersions, comma-separated")
    parser.add_argument("--sizes", default="1000,2000,4000,7000,11000,16000")
    parser.add_argument("--band", type=int, default=2000, help="Rank band width")
    parser.add_argument("--out", required=True, help="Output prefix (.tsv and .md)")
    args = parser.parse_args(argv)

    with open(args.manifest, "rb") as handle:
        ds = [d for d in tomllib.load(handle)["dataset"] if d["name"] == args.dataset][0]
    path = lambda p: os.path.join(args.data_dir, p)  # noqa: E731
    cache = args.cache_dir or args.out + "_cache"
    mate, mapq, chrom = ds.get("mate", 1), ds.get("min_mapq", 0), ds["chromosome"]

    def load(p):
        return {k: v for k, v in harness.load_signal(path(p), mate, mapq, cache).items()
                if k[0] == chrom}

    files = list(ds["enriched"])
    if args.with_technical_replicate and ds.get("technical_replicate"):
        files.append(ds["technical_replicate"])
    enriched = [load(p) for p in files]
    controls = [load(p) for p in ds.get("control", [])]
    exclude = ct.GeneIndex(ct.read_genes(path(ds["annotation"]), ("rRNA", "tRNA")))
    seq = harness.read_fasta(path(ds["genome"]))[chrom]
    background = harness.MatchedBackground(seq)
    reference = None
    if ds.get("reference"):
        reference = [(st, p) for s, st, p in ct.read_sites(path(ds["reference"]))
                     if s in (chrom, "gi|48994873|gb|U00096.2|")]

    # Every candidate is tested (fdr 1) so each method ranks the same set.
    local = ct.call_local(enriched, controls or None, fdr=1.0, min_reads=args.min_reads,
                          exclude=exclude)
    pooled = ct.pool_signals(enriched)
    count = {t.key: pooled[(t.seqid, t.strand)][t.pos] for t in local}
    pval = {t.key: t.pvalue for t in local}
    keys = list(count)
    print("%d candidate positions" % len(keys), file=sys.stderr)
    rankings = {"rpm": sorted(keys, key=lambda k: -count[k]),
                "local": sorted(keys, key=lambda k: (pval[k], -count[k]))}
    stops = {"local": len(top_sites([t.key for t in local if t.qvalue <= 0.05]))}
    ratio = {t.key: count[t.key] / t.background for t in local}
    rankings["ratio"] = sorted(keys, key=lambda k: (-ratio[k], -count[k]))
    for d in [float(x) for x in args.local_dispersions.split(",") if x]:
        fixed = ct.call_local(enriched, controls or None, fdr=1.0, min_reads=args.min_reads,
                              exclude=exclude, dispersion=d)
        fp = {t.key: (t.pvalue, t.qvalue) for t in fixed}
        name = "local_d%g" % d
        rankings[name] = sorted(fp, key=lambda k: (fp[k][0], -count[k]))
        stops[name] = len(top_sites([k for k in fp if fp[k][1] <= 0.05]))
    ett = harness.method_ettwiller2016(enriched, controls)
    stops["rpm"] = len([t for t in ett if t.seqid == chrom])
    if controls and len(enriched) == 1:
        cand = set(keys)
        for d in [float(x) for x in args.dispersions.split(",")]:
            calls = ct.call_binomial(enriched, controls, fdr=1.0, min_reads=args.min_reads,
                                     exclude=exclude, dispersion=d)
            cp = {t.key: (t.pvalue, t.qvalue) for t in calls if t.key in cand}
            name = "control_d%g" % d
            rankings[name] = sorted(cp, key=lambda k: (cp[k][0], -count[k]))
            stops[name] = len(top_sites([k for k in cp if cp[k][1] <= 0.05]))

    out = []
    for n in [int(x) for x in args.sizes.split(",")]:
        for name, ranked in rankings.items():
            row = dict(dataset="%s_top%d" % (args.dataset, n), method=name)
            row.update(score(top_sites(ranked, n), seq, exclude, reference, background))
            out.append(row)
            print("  top %-6d %-12s core10 %.3f" % (n, name, row["core10"]), file=sys.stderr)
    harness.write_rows(out, args.out)

    # Box rate by rank band, and where each method's own threshold stops.
    full = {name: top_sites(ranked) for name, ranked in rankings.items()}
    longest = max(len(v) for v in full.values())
    with open(args.out + ".bands.tsv", "w") as handle:
        names = list(full)
        handle.write("ranks\t" + "\t".join(names) + "\n")
        for start in range(0, longest, args.band):
            cells = []
            for name in names:
                band = full[name][start:start + args.band]
                rate = (harness.promoter_rates([(st, p) for _, st, p in band], seq)["core10"]
                        if band else None)
                mark = "*" if start < stops.get(name, 0) <= start + args.band else ""
                cells.append("" if rate is None else "%.3f%s" % (rate, mark))
            handle.write("%d-%d\t%s\n" % (start + 1, start + args.band, "\t".join(cells)))
        handle.write("stop\t" + "\t".join(str(stops.get(n, "")) for n in names) + "\n")
    print("wrote %s.bands.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
