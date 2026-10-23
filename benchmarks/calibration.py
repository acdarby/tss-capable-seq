#!/usr/bin/env python3
"""Calibration tests for statistical TSS calling (docs/benchmark_plan.md, section 5).

A test whose false-discovery rate is honest should call (almost) nothing when
there is nothing to find. Three null situations:

1. simulated null: enriched and control libraries drawn from the same rates
   (TSS present but not enriched), with or without extra variation between
   libraries; every call is false;
2. split library: one real library split at random into two halves, one used
   as "enriched" and the other as "control";
3. replicate against replicate: one enriched replicate used as "enriched" and
   another as "control"; this null keeps the real variation between libraries;
4. label swap: the real control used as "enriched" and the enriched library as
   "control"; anything called should be processed ends (rRNA/tRNA, depleted
   in Cappable-seq), not promoters.

Tests 3 and 4 need two enriched replicates and a control library respectively.

usage: calibration.py --data-dir DIR --dataset ecoli --dispersions 0,0.05 --out results/calibration
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
import simulate  # noqa: E402


def pvalue_summary(enriched, controls, exclude, min_fold=1.0, dispersion=0.0):
    """Share of tested positions with p below 0.05 and 0.01, and calls at FDR 5% and 1%."""
    kw = dict(min_fold=min_fold, exclude=exclude, dispersion=dispersion)
    calls5 = ct.call_binomial(enriched, controls, fdr=0.05, **kw)
    calls1 = ct.call_binomial(enriched, controls, fdr=0.01, **kw)
    # Re-run the per-position tests to get all p-values (calls only keep significant ones).
    allp = ct.call_binomial(enriched, controls, fdr=1.0, min_replicates=0, **kw)
    p = np.array([t.pvalue for t in allp])
    return {"tested": len(p),
            "p_below_0.05": float((p < 0.05).mean()) if len(p) else float("nan"),
            "p_below_0.01": float((p < 0.01).mean()) if len(p) else float("nan"),
            "calls_fdr5": len(calls5), "calls_fdr1": len(calls1)}, calls5


def split(signal, rng):
    """Split a library into two random halves (binomial thinning with p = 0.5)."""
    a, b = {}, {}
    for key, d in signal.items():
        pos = np.fromiter(d.keys(), dtype=np.int64, count=len(d))
        cnt = np.fromiter(d.values(), dtype=np.int64, count=len(d))
        half = rng.binomial(cnt, 0.5)
        rest = cnt - half
        a[key] = {int(p): float(v) for p, v in zip(pos, half) if v}
        b[key] = {int(p): float(v) for p, v in zip(pos, rest) if v}
    return a, b


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--dataset", default="ecoli")
    parser.add_argument("--manifest", default=os.path.join(HERE, "datasets.toml"))
    parser.add_argument("--cache-dir")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--dispersions", default="0,0.05",
                        help="Beta-binomial dispersions to test, comma-separated (0 = binomial)")
    parser.add_argument("--out", required=True, help="Output prefix (.tsv)")
    args = parser.parse_args(argv)
    rng = np.random.default_rng(args.seed)
    rows = []

    # Real data (tests 2-4).
    with open(args.manifest, "rb") as handle:
        ds = [d for d in tomllib.load(handle)["dataset"] if d["name"] == args.dataset][0]
    path = lambda p: os.path.join(args.data_dir, p)  # noqa: E731
    cache = args.cache_dir or args.out + "_cache"
    mate, mapq, chrom = ds.get("mate", 1), ds.get("min_mapq", 0), ds["chromosome"]
    enr = harness.load_signal(path(ds["enriched"][0]), mate, mapq, cache)
    enr2 = (harness.load_signal(path(ds["enriched"][1]), mate, mapq, cache)
            if len(ds["enriched"]) > 1 else None)
    ctl = (harness.load_signal(path(ds["control"][0]), mate, mapq, cache)
           if ds.get("control") else None)
    if enr2 is not None:
        print("Dispersion estimated from enriched replicates: %s"
              % ct.estimate_dispersion([enr, enr2]), file=sys.stderr)
    rna = ct.read_genes(path(ds["annotation"]), ("rRNA", "tRNA"))
    exclude = ct.GeneIndex(rna)
    seq = harness.read_fasta(path(ds["genome"]))[chrom]

    for disp in [float(x) for x in args.dispersions.split(",")]:
        # 1. Simulated null: TSS equally strong in enriched and control, no processed
        #    ends, no rRNA; two depths, with and without variation between libraries.
        for depth, null_disp in ((1.0, 0.0), (0.25, 0.0), (1.0, 0.05), (1.0, 0.1)):
            sim = simulate.simulate(length=1_000_000, n_tss=1000, replicates=1, depth=depth,
                                    tss_mu=4.0, background=0.2, with_control=True,
                                    dispersion=null_disp, control_tss_factor=1.0,
                                    n_processed=0, rrna_fraction=0.0, seed=args.seed)
            res, _ = pvalue_summary(sim["enriched"], sim["controls"], None, dispersion=disp)
            rows.append(dict(test="simulated null (no enrichment), depth %g, true dispersion %g"
                                  % (depth, null_disp), dispersion=disp, **res))
            print(rows[-1], file=sys.stderr)

        for label, lib in (("enriched", enr), ("control", ctl)):
            if lib is None:
                continue
            a, b = split(lib, rng)
            res, _ = pvalue_summary([a], [b], exclude, dispersion=disp)
            rows.append(dict(test="split %s library into halves" % label, dispersion=disp, **res))
            print(rows[-1], file=sys.stderr)

        if enr2 is not None:
            res, _ = pvalue_summary([enr], [enr2], exclude, dispersion=disp)
            rows.append(dict(test="enriched replicate 1 against replicate 2", dispersion=disp,
                             **res))
            print(rows[-1], file=sys.stderr)

        if ctl is None:
            continue
        res, calls = pvalue_summary([ctl], [enr], exclude, dispersion=disp)
        calls = ct.cluster_tss(calls, 6)
        in_rna = sum(1 for t in calls if exclude.overlapping(t.seqid, t.pos - 5, t.pos + 5))
        sites = [(t.strand, t.pos) for t in calls if t.seqid == chrom]
        box = harness.box_rates(sites, seq)[0]
        res.update(clustered=len(calls), share_in_rRNA_tRNA=in_rna / len(calls) if calls else float("nan"),
                   box10=box)
        rows.append(dict(test="label swap (control as enriched)", dispersion=disp, **res))
        print(rows[-1], file=sys.stderr)

        # Reference for the label swap: the real orientation.
        res, calls = pvalue_summary([enr], [ctl], exclude, dispersion=disp)
        calls = ct.cluster_tss(calls, 6)
        in_rna = sum(1 for t in calls if exclude.overlapping(t.seqid, t.pos - 5, t.pos + 5))
        sites = [(t.strand, t.pos) for t in calls if t.seqid == chrom]
        res.update(clustered=len(calls), share_in_rRNA_tRNA=in_rna / len(calls) if calls else float("nan"),
                   box10=harness.box_rates(sites, seq)[0])
        rows.append(dict(test="real orientation (for comparison)", dispersion=disp, **res))
        print(rows[-1], file=sys.stderr)

    harness.write_rows([dict(dataset=r.pop("test"),
                             method="betabinom" if r["dispersion"] else "binomial", **r)
                        for r in rows], args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
