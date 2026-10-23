#!/usr/bin/env python3
"""Benchmark of shoulder removal (approach G).

A call is a shoulder if a stronger call on the same strand lies within W nt
upstream and it has less than r times that call's RPM
(cappable_tss.remove_shoulders). For a grid of W and r this reports how many
calls are removed and the evidence for removed and kept calls:

* core10: sigma70 -10 core (harness.promoter_elements); core10_bg for
  composition-matched random positions;
* core10_same_rpm: kept calls with the RPM closest to each removed call, so
  that removed shoulders are compared with equally weak calls that are not
  shoulders;
* model_pos: share with a positive genome-trained promoter score
  (cappable_tss.promoter_scores);
* ref_recall_0 (E. coli Table S1, exact base) before and after;
* ortho_rise (orthogonal fragment library, wAlbB) of removed and kept.

usage: shoulders.py --data-dir DIR --dataset ecoli --methods ettwiller2016,local --out PREFIX
"""
import argparse
import bisect
import os
import sys
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
import cappable_tss as ct  # noqa: E402
import harness  # noqa: E402


def share(values):
    values = [v for v in values if v is not None]
    return sum(values) / len(values) if values else float("nan")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--dataset", default="ecoli")
    parser.add_argument("--manifest", default=os.path.join(HERE, "datasets.toml"))
    parser.add_argument("--cache-dir")
    parser.add_argument("--methods", default="ettwiller2016,local")
    parser.add_argument("--windows", default="10,20,30,50")
    parser.add_argument("--ratios", default="0.05,0.1,0.2,0.3")
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
    genome = {chrom: ct.read_fasta(path(ds["genome"]))[chrom]}
    seq = genome[chrom]
    background = harness.MatchedBackground(seq)
    reference = None
    if ds.get("reference"):
        reference = [(st, p) for s, st, p in ct.read_sites(path(ds["reference"]))
                     if s in (chrom, "gi|48994873|gb|U00096.2|")]
    cov = None
    if ds.get("orthogonal"):
        cov = harness.coverage_track(path(ds["orthogonal"]), chrom, ds.get("orthogonal_mate", 1),
                                     mapq, cache)

    def core(t):
        e = harness.promoter_elements(seq, t.strand, t.pos)
        return e["core10"] if e else None

    rows = []
    for method in args.methods.split(","):
        calls = [t for t in harness.METHODS[method](enriched, controls, {"exclude": exclude})
                 if t.seqid == chrom]
        ranked = sorted(calls, key=lambda t: -t.enriched_rpm)
        scores = ct.promoter_scores(genome, [t.key for t in ranked[:len(ranked) // 3]],
                                    [t.key for t in calls])
        base = dict(dataset=args.dataset, method=method, window="", ratio="",
                    n_tss=len(calls), core10_all=share(core(t) for t in calls))
        if reference:
            base["ref_recall_0"] = harness.match(reference, harness.as_sites(calls), 0) / len(reference)
        rows.append(base)
        for window in [int(x) for x in args.windows.split(",")]:
            ct.mark_shoulders(calls, window)
            for ratio in [float(x) for x in args.ratios.split(",")]:
                removed = [t for t in calls if t.peak_ratio < ratio]
                kept = [t for t in calls if t.peak_ratio >= ratio]
                kept_rpm = sorted((t.enriched_rpm, i) for i, t in enumerate(kept))
                rpms = [r for r, _ in kept_rpm]
                same = []
                for t in removed:
                    j = min(bisect.bisect_left(rpms, t.enriched_rpm), len(rpms) - 1)
                    same.append(kept[kept_rpm[j][1]])
                row = dict(dataset=args.dataset, method=method, window=window, ratio=ratio,
                           n_tss=len(kept), removed=len(removed),
                           core10_removed=share(core(t) for t in removed),
                           core10_removed_bg=harness.promoter_rates(
                               harness.as_sites(removed), seq, background).get("core10_bg"),
                           core10_same_rpm=share(core(t) for t in same),
                           core10_kept=share(core(t) for t in kept),
                           model_pos_removed=share(scores.get(t.key, 0) > 0 for t in removed),
                           model_pos_same_rpm=share(scores.get(t.key, 0) > 0 for t in same))
                if reference:
                    row["ref_recall_0"] = harness.match(reference, harness.as_sites(kept), 0) / len(reference)
                if cov is not None:
                    row["ortho_removed"] = harness.coverage_rise(harness.as_sites(removed), cov)
                    row["ortho_same_rpm"] = harness.coverage_rise(harness.as_sites(same), cov)
                rows.append(row)
                print(row, file=sys.stderr)
    harness.write_rows(rows, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
