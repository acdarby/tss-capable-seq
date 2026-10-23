#!/usr/bin/env python3
"""Validation of the standard mode (`call --mode standard`).

Runs cappable_tss.standard_calls on one dataset, as the tool does, and judges
each confidence tier on evidence that the tiers do not use directly:

* core10 / core10_bg: sigma70 -10 core (harness.promoter_elements) and
  composition-matched background (only partly independent: with a genome the
  tiers use the genome-trained promoter score);
* ortho_rise: coverage rise in the orthogonal fragment library (wAlbB);
* confirmed: the reads are thinned to `--fraction`, the standard mode is run
  on the thinned data, and a thinned call counts as confirmed if the
  full-data standard mode calls it (+-2 nt) with >= 10 reads.

usage: standard_mode.py --data-dir DIR --dataset walbb --out results/standard_walbb
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

TIERS = ("high", "medium", "low")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--dataset", default="walbb")
    parser.add_argument("--manifest", default=os.path.join(HERE, "datasets.toml"))
    parser.add_argument("--cache-dir")
    parser.add_argument("--fraction", type=float, default=0.2)
    parser.add_argument("--no-genome", action="store_true")
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
    genome = None if args.no_genome else {chrom: seq}
    background = harness.MatchedBackground(seq)
    cov = None
    if ds.get("orthogonal"):
        cov = harness.coverage_track(path(ds["orthogonal"]), chrom, ds.get("orthogonal_mate", 1),
                                     mapq, cache)

    calls, design, _, _ = ct.standard_calls(enriched, controls, exclude, genome)
    print("design", design, file=sys.stderr)
    rng = np.random.default_rng(args.seed)
    thinned = [harness.thin(lib, args.fraction, rng) for lib in enriched]
    tcontrols = [harness.thin(lib, args.fraction, rng) for lib in controls]
    tcalls, _, _, _ = ct.standard_calls(thinned, tcontrols, exclude, genome)
    full = {(t.strand, t.pos + d) for t in calls if t.reads >= 10 for d in range(-2, 3)}

    rows = []
    label = args.dataset + ("_nogenome" if args.no_genome else "")
    for tier in TIERS + ("all",):
        group = [t for t in calls if tier in ("all", t.confidence)]
        tgroup = [t for t in tcalls if tier in ("all", t.confidence)]
        sites = harness.as_sites(group)
        prom = harness.promoter_rates(sites, seq, background)
        row = dict(dataset=label, method="%s_%s" % (design["method"], tier), n_tss=len(group),
                   share=len(group) / len(calls) if calls else float("nan"),
                   core10=prom["core10"], core10_bg=prom.get("core10_bg"),
                   n_thinned=len(tgroup),
                   confirmed=(sum((t.strand, t.pos) in full for t in tgroup) / len(tgroup)
                              if tgroup else float("nan")))
        if cov is not None:
            row["ortho_rise"] = harness.coverage_rise(sites, cov)
            row["ortho_rise_bg"] = harness.coverage_rise(harness.shifted(sites), cov)
        rows.append(row)
        print(row, file=sys.stderr)
    harness.write_rows(rows, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
