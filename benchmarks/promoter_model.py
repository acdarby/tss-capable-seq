#!/usr/bin/env python3
"""Benchmark of the genome-trained promoter score (approach H).

For one dataset, with the local-background calls (cappable_tss.call_local):

1. recognition: AUROC of the cross-validated promoter score for strong (>= 20
   pooled reads) and weak (< 10) calls against composition-matched random
   positions, next to the fixed sigma70 score of harness.promoter_elements;
2. rescue at low depth: the reads are thinned to `--fraction`, calls and
   model come from the thinned data only, and weak candidates there (2-5
   pooled reads, the highest position within +-5 nt) are scored; the truth is
   whether the full data call them with >= 10 reads. AUROC of read count,
   promoter score and both (logistic regression, 5-fold);
3. by score quartile of the full-data calls: share confirmed by the
   orthogonal fragment library (coverage rise, wAlbB) where available;
4. the strongest learned weights.

Needs numpy and scikit-learn (benchmarks only; the tool itself does not).

usage: promoter_model.py --data-dir DIR --dataset walbb --out results/promoter_model_walbb
"""
import argparse
import math
import os
import sys
import tomllib

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_predict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
import cappable_tss as ct  # noqa: E402
import harness  # noqa: E402


def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    return roc_auc_score([1] * len(pos) + [0] * len(neg), list(pos) + list(neg))


def fixed_score(seq, strand, pos):
    e = harness.promoter_elements(seq, strand, pos)
    if not e:
        return None
    return 2 * e["core10"] + e["box10"] + e["box35"] + 0.5 * e["ext10"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--dataset", default="walbb")
    parser.add_argument("--manifest", default=os.path.join(HERE, "datasets.toml"))
    parser.add_argument("--cache-dir")
    parser.add_argument("--fraction", type=float, default=0.2)
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
    controls = [load(p) for p in ds.get("control", [])] or None
    exclude = ct.GeneIndex(ct.read_genes(path(ds["annotation"]), ("rRNA", "tRNA")))
    genome = {chrom: ct.read_fasta(path(ds["genome"]))[chrom]}
    seq = genome[chrom]
    rows = []

    def calls_and_scores(libs, extra_sites=()):
        calls = ct.cluster_tss(ct.call_local(libs, controls, exclude=exclude), 6)
        ranked = sorted(calls, key=lambda t: -t.enriched_rpm)
        positives = [t.key for t in ranked[:len(ranked) // 3]]
        scores = ct.promoter_scores(genome, positives,
                                    [t.key for t in calls] + list(extra_sites), seed=args.seed)
        return calls, scores, positives

    # 1. Recognition on the full data.
    pooled = ct.pool_signals(enriched)
    calls, _, positives = calls_and_scores(enriched)
    count = {t.key: pooled[(t.seqid, t.strand)][t.pos] for t in calls}
    groups = {"strong": [t.key for t in calls if count[t.key] >= 20],
              "weak": [t.key for t in calls if count[t.key] < 10]}
    background = {g: ct.matched_background(genome, keys, per_site=1, seed=args.seed + 7)
                  for g, keys in groups.items()}
    extra = [k for v in background.values() for k in v]
    scores = ct.promoter_scores(genome, positives, [t.key for t in calls] + extra, seed=args.seed)
    for g, keys in groups.items():
        fx_pos = [fixed_score(seq, s, p) for _, s, p in keys]
        fx_neg = [fixed_score(seq, s, p) for _, s, p in background[g]]
        rows.append(dict(dataset=args.dataset, method="recognition_" + g, n=len(keys),
                         auroc_model=auroc([scores[k] for k in keys if k in scores],
                                           [scores[k] for k in background[g] if k in scores]),
                         auroc_fixed=auroc([x for x in fx_pos if x is not None],
                                           [x for x in fx_neg if x is not None])))
        print(rows[-1], file=sys.stderr)

    # 3. Orthogonal support by score quartile (full-data calls).
    if ds.get("orthogonal"):
        cov = harness.coverage_track(path(ds["orthogonal"]), chrom, ds.get("orthogonal_mate", 1),
                                     mapq, cache)
        scored = sorted((scores[t.key], (t.strand, t.pos)) for t in calls if t.key in scores)
        q = len(scored) // 4
        for i in range(4):
            part = [site for _, site in scored[i * q:(i + 1) * q if i < 3 else None]]
            rows.append(dict(dataset=args.dataset, method="ortho_rise_score_q%d" % (i + 1),
                             n=len(part), ortho_rise=harness.coverage_rise(part, cov)))
            print(rows[-1], file=sys.stderr)

    # 2. Rescue at low depth.
    rng = np.random.default_rng(args.seed)
    thinned = [harness.thin(lib, args.fraction, rng) for lib in enriched]
    tpooled = ct.pool_signals(thinned)
    cand = []
    for (seqid, strand), d in tpooled.items():
        for pos, x in d.items():
            if 2 <= x <= 5 and all(d.get(q, 0) <= x for q in range(pos - 5, pos + 6)):
                cand.append((seqid, strand, pos, x))
    _, tscores, _ = calls_and_scores(thinned, [c[:3] for c in cand])
    full = {(t.strand, t.pos + d) for t in calls if count[t.key] >= 10 for d in range(-2, 3)}
    cand = [c for c in cand if c[:3] in tscores]
    y = [int((s, p) in full) for _, s, p, _ in cand]
    reads = np.array([[math.log(x)] for *_, x in cand])
    both = np.array([[math.log(c[3]), tscores[c[:3]]] for c in cand])
    combined = cross_val_predict(LogisticRegression(), both, y, cv=5, method="decision_function")
    rows.append(dict(dataset=args.dataset, method="rescue_%g" % args.fraction, n=len(cand),
                     confirmed=float(np.mean(y)), auroc_reads=roc_auc_score(y, reads[:, 0]),
                     auroc_model=roc_auc_score(y, both[:, 1]),
                     auroc_both=roc_auc_score(y, combined)))
    print(rows[-1], file=sys.stderr)
    # Confirmation by score quartile among the weak candidates.
    order = np.argsort(both[:, 1])
    q = len(order) // 4
    for i in range(4):
        idx = order[i * q:(i + 1) * q if i < 3 else None]
        rows.append(dict(dataset=args.dataset, method="rescue_score_q%d" % (i + 1), n=len(idx),
                         confirmed=float(np.mean([y[j] for j in idx]))))
    harness.write_rows(rows, args.out)

    # 4. What the model learned (trained on all positives).
    negatives = ct.matched_background(genome, positives, seed=args.seed)
    data = [(ct.promoter_features(w), lab) for lab, grp in ((1, positives), (0, negatives))
            for w in (ct.upstream_window(seq, s, p) for _, s, p in grp) if w]
    model = ct.LogisticModel(4 * ct.PROMOTER_UP + 16 * (ct.PROMOTER_UP - 1)).fit(
        [f for f, _ in data], [lab for _, lab in data])
    mono = [(model.w[4 * i + j], i - ct.PROMOTER_UP, "ACGT"[j])
            for i in range(ct.PROMOTER_UP) for j in range(4)]
    with open(args.out + ".weights.txt", "w") as handle:
        handle.write("Strongest single-base weights (position relative to the TSS):\n")
        for w, pos, base in sorted(mono, key=lambda t: -abs(t[0]))[:20]:
            handle.write("%+4d %s %+.2f\n" % (pos, base, w))
        handle.write("\nConsensus of positive weights, -60..-1 (upper case > 0.15):\n")
        cons = []
        for i in range(ct.PROMOTER_UP):
            w, j = max((model.w[4 * i + j], j) for j in range(4))
            cons.append("ACGT"[j] if w > 0.15 else "acgt"[j] if w > 0.05 else ".")
        handle.write("".join(cons) + "\n")
    print("wrote %s.weights.txt" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
