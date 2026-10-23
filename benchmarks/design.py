#!/usr/bin/env python3
"""Experimental-design comparisons for the recommendations (docs/design_recommendations.md).

Runs the standard mode (cappable_tss.standard_calls, with the genome) on
variants of one dataset and reports TSS numbers, confidence tiers and
evidence for each:

* --variants replicates: the enriched replicates as given, and the same reads
  pooled into one library (same total depth, no replicates);
* --variants depth: all reads thinned to 100, 50, 25 and 10%;
* --variants control: with and without the control library.

Evidence: -10 core rate (harness.promoter_elements) with a composition-matched
background, share within 5 nt of rRNA/tRNA genes, the wAlbB coverage rise
where available, and, against the full-data standard calls (>= 10 reads,
+-2 nt), the share of calls confirmed ("confirmed_by_full"; for the full data
itself this is trivially high and not reported).

usage: design.py --data-dir DIR --dataset walbb --variants replicates,depth --out PREFIX
"""
import argparse
import os
import sys
import tomllib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
import cappable_tss as ct  # noqa: E402
import harness  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--dataset", default="walbb")
    parser.add_argument("--manifest", default=os.path.join(HERE, "datasets.toml"))
    parser.add_argument("--cache-dir")
    parser.add_argument("--variants", default="replicates,depth,control")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    with open(args.manifest, "rb") as handle:
        ds = [d for d in tomllib.load(handle)["dataset"] if d["name"] == args.dataset][0]
    path = lambda p: os.path.join(args.data_dir, p)  # noqa: E731
    cache = args.cache_dir or args.out + "_cache"
    mate, mapq, chrom = ds.get("mate", 1), ds.get("min_mapq", 0), ds["chromosome"]

    def load(p):
        return {k: v for k, v in harness.load_signal(path(p), mate, mapq, cache).items()
                if k[0] == chrom}

    enriched = [load(p) for p in ds["enriched"]]
    controls = [load(p) for p in ds.get("control", [])]
    exclude = ct.GeneIndex(ct.read_genes(path(ds["annotation"]), ("rRNA", "tRNA")))
    seq = harness.read_fasta(path(ds["genome"]))[chrom]
    genome = {chrom: seq}
    background = harness.MatchedBackground(seq)
    cov = None
    if ds.get("orthogonal"):
        cov = harness.coverage_track(path(ds["orthogonal"]), chrom, ds.get("orthogonal_mate", 1),
                                     mapq, cache)
    reference, _, _, _ = ct.standard_calls(enriched, controls, exclude, genome)
    full = {(t.strand, t.pos + d) for t in reference if t.reads >= 10 for d in range(-2, 3)}
    rng = np.random.default_rng(args.seed)

    def evaluate(label, libs, ctls, compare=True):
        calls, design, _, _ = ct.standard_calls(libs, ctls, exclude, genome)
        sites = harness.as_sites(calls)
        prom = harness.promoter_rates(sites, seq, background)
        tiers = {k: sum(t.confidence == k for t in calls) for k in ("high", "medium", "low")}
        high = [t for t in calls if t.confidence == "high"]
        row = dict(dataset=args.dataset, method=label, design=design["method"],
                   reads=int(sum(ct.total_signal(lib) for lib in libs)),
                   n_tss=len(calls), high=tiers["high"], medium=tiers["medium"],
                   low=tiers["low"], core10=prom["core10"], core10_bg=prom.get("core10_bg"),
                   core10_high=harness.promoter_rates(harness.as_sites(high), seq)["core10"],
                   in_rna=sum(1 for t in calls if exclude.overlapping(t.seqid, t.pos - 5, t.pos + 5))
                   / len(calls) if calls else float("nan"))
        if compare:
            row["confirmed_by_full"] = (sum((t.strand, t.pos) in full for t in calls) / len(calls)
                                        if calls else float("nan"))
        if cov is not None:
            row["ortho_rise"] = harness.coverage_rise(sites, cov)
        print(row, file=sys.stderr)
        return row

    rows = []
    variants = args.variants.split(",")
    if "replicates" in variants and len(enriched) > 1:
        rows.append(evaluate("%d replicates" % len(enriched), enriched, controls, compare=False))
        rows.append(evaluate("same reads as 1 library", [ct.pool_signals(enriched)],
                             [ct.pool_signals(controls)] if controls else [], compare=False))
    if "depth" in variants:
        for frac in (0.5, 0.25, 0.1):
            libs = [harness.thin(lib, frac, rng) for lib in enriched]
            ctls = [harness.thin(lib, frac, rng) for lib in controls]
            rows.append(evaluate("%g%% of reads" % (100 * frac), libs, ctls))
    if "control" in variants and controls:
        rows.append(evaluate("with control", enriched, controls, compare=False))
        rows.append(evaluate("without control", enriched, [], compare=False))
    harness.write_rows(rows, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
