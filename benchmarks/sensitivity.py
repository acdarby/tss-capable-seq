#!/usr/bin/env python3
"""Parameter sensitivity of the local-background test (plan section 5).

For each setting, reports the TSS called at FDR 5% (clustered, 6 nt), their
-10 core rate (harness.promoter_elements), and that of the top N sites
ranked by p-value.

usage: sensitivity.py --data-dir DIR --dataset walbb --out results/sensitivity_walbb
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
from matched import top_sites  # noqa: E402

SETTINGS = [dict()] + [dict(window=w) for w in (50, 200)] + \
    [dict(trim=t) for t in (0, 10)] + [dict(min_reads=m) for m in (2, 5)] + \
    [dict(fdr=f) for f in (0.01, 0.1)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--dataset", default="walbb")
    parser.add_argument("--manifest", default=os.path.join(HERE, "datasets.toml"))
    parser.add_argument("--cache-dir")
    parser.add_argument("--top", type=int, default=4000)
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
    controls = [load(p) for p in ds.get("control", [])] or None
    exclude = ct.GeneIndex(ct.read_genes(path(ds["annotation"]), ("rRNA", "tRNA")))
    seq = harness.read_fasta(path(ds["genome"]))[chrom]
    rows = []
    for setting in SETTINGS:
        kw = dict(setting)
        fdr = kw.pop("fdr", 0.05)
        tested = ct.call_local(enriched, controls, fdr=1.0, exclude=exclude, **kw)
        called = ct.cluster_tss([t for t in tested if t.qvalue <= fdr], 6)
        ranked = [t.key for t in sorted(tested, key=lambda t: (t.pvalue, -t.enriched_rpm))]
        top = top_sites(ranked, args.top)
        label = ", ".join("%s=%s" % kv for kv in setting.items()) or "defaults"
        rows.append(dict(dataset=args.dataset, method=label, n_tss=len(called),
                         core10_called=harness.promoter_rates(harness.as_sites(called), seq)["core10"],
                         **{"core10_top%d" % args.top:
                            harness.promoter_rates([(st, p) for _, st, p in top], seq)["core10"]}))
        print(rows[-1], file=sys.stderr)
    harness.write_rows(rows, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
