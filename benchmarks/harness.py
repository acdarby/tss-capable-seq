#!/usr/bin/env python3
"""Benchmark harness for TSS-calling methods (see docs/benchmark_plan.md).

Runs every method on every dataset in datasets.toml (real data) and on
simulated libraries, and writes one table of measures.

usage:
  harness.py real --data-dir DIR [--datasets ecoli,walbb] [--methods ...]
                  [--downsample 0.5,0.25] --out results/real
  harness.py sim  [--depths 1,0.25,0.05] [--methods ...] --out results/sim

Measures on real data (no single ground truth, so several lines of evidence):
  n_tss, class shares, −10/−35 box rates (and a shifted-position background),
  agreement with a reference TSS set, leave-one-replicate-out support,
  recovery after downsampling, and coverage rise in orthogonal fragment-library
  data. On simulated data: precision and recall against the planted TSS.
"""
import argparse
import collections
import math
import os
import pickle
import statistics
import sys
import tomllib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
import cappable_tss as ct  # noqa: E402
import simulate  # noqa: E402

CLASSES = ("Primary", "Secondary", "Internal", "Antisense", "Orphan")


# ---------------------------------------------------------------------------
# Methods: each takes (enriched signals, control signals) and returns TSS
# objects; n replicates are passed implicitly as len(enriched).
# ---------------------------------------------------------------------------

def method_ettwiller2016(enriched, controls, ctx=None):
    """Published Cappable-seq pipeline; with replicates, TSS in all of them."""
    preset = dict(ct.PRESETS["ettwiller2016"])
    distance = preset.pop("cluster_distance")
    calls = ct.call_positions(enriched, controls, min_replicates=len(enriched), **preset)
    return ct.cluster_tss(calls, distance)


def method_suliman2025(enriched, controls, ctx=None):
    """Suliman et al. 2025: >= 10 CPM, 10 nt clusters, all replicates, no control."""
    calls = ct.call_positions(enriched, [], min_rpm=10.0, min_replicates=len(enriched))
    return ct.cluster_tss(calls, 10)


def _binomial(enriched, controls, ctx, min_fold, dispersion=0.0):
    if not controls:
        return None                       # needs a control library
    if dispersion == "auto":              # as `call --dispersion auto`
        dispersion = ct.estimate_dispersion(enriched)
        if dispersion is None:
            dispersion = ct.DEFAULT_DISPERSION
    calls = ct.call_binomial(enriched, controls, fdr=0.05, min_fold=min_fold,
                             exclude=(ctx or {}).get("exclude"), dispersion=dispersion)
    return ct.cluster_tss(calls, 6)


def method_binomial(enriched, controls, ctx=None):
    """Binomial test against the control, FDR 5%, any enrichment, rRNA/tRNA excluded."""
    return _binomial(enriched, controls, ctx, 1.0)


def method_binomial_fold2(enriched, controls, ctx=None):
    """As binomial, but the null hypothesis is enrichment of at most 2-fold."""
    return _binomial(enriched, controls, ctx, 2.0)


def method_betabinom(enriched, controls, ctx=None):
    """Beta-binomial test (dispersion estimated from replicates, else 0.05), FDR 5%."""
    return _binomial(enriched, controls, ctx, 1.0, "auto")


def method_betabinom_fold2(enriched, controls, ctx=None):
    """As betabinom, but the null hypothesis is enrichment of at most 2-fold."""
    return _binomial(enriched, controls, ctx, 2.0, "auto")


def method_local(enriched, controls, ctx=None):
    """5'-end pile-ups above the local background (no control needed), FDR 5%."""
    calls = ct.call_local(enriched, controls or None, exclude=(ctx or {}).get("exclude"))
    return ct.cluster_tss(calls, 6)


METHODS = {
    "ettwiller2016": method_ettwiller2016,
    "suliman2025": method_suliman2025,
    "binomial": method_binomial,
    "binomial_fold2": method_binomial_fold2,
    "betabinom": method_betabinom,
    "betabinom_fold2": method_betabinom_fold2,
    "local": method_local,
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_signal(path, mate, min_mapq, cache_dir):
    """5'-end signal of a BAM, cached as a pickle (reading BAMs is slow)."""
    key = "%s.m%d.q%d.pkl" % (os.path.basename(path), mate, min_mapq)
    cache = os.path.join(cache_dir, key) if cache_dir else None
    if cache and os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime(path):
        with open(cache, "rb") as handle:
            return pickle.load(handle)
    sig = ct.five_prime_from_bam(path, mate=mate, min_mapq=min_mapq)
    sig = {k: dict(v) for k, v in sig.items()}
    if cache:
        os.makedirs(cache_dir, exist_ok=True)
        with open(cache, "wb") as handle:
            pickle.dump(sig, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return sig


def read_fasta(path):
    seqs, name, parts = {}, None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip().upper())
    seqs[name] = "".join(parts)
    return seqs


def thin(signal, fraction, rng):
    """Binomial downsampling of a signal to `fraction` of its reads."""
    out = {}
    for key, d in signal.items():
        pos = np.fromiter(d.keys(), dtype=np.int64, count=len(d))
        cnt = np.fromiter(d.values(), dtype=np.int64, count=len(d))
        kept = rng.binomial(cnt, fraction)
        nz = kept > 0
        out[key] = dict(zip(pos[nz].tolist(), kept[nz].astype(float).tolist()))
    return out


# ---------------------------------------------------------------------------
# Measures
# ---------------------------------------------------------------------------

COMP = str.maketrans("ACGT", "TGCA")


def upstream(seq, strand, pos, a, b):
    """Bases from −a to −b upstream of a 1-based TSS, in transcript direction."""
    if strand == "+":
        return seq[pos - 1 - a:pos - 1 - b + 1]
    return seq[pos + b - 1:pos + a][::-1].translate(COMP)


def has_box(window, motif):
    k = len(motif)
    return any(sum(x != y for x, y in zip(window[i:i + k], motif)) <= 1
               for i in range(len(window) - k + 1))


def box_rates(sites, seq):
    """Share of sites with a −10 (TATAAT, −15..−5) and −35 (TTGACA, −40..−28) box."""
    ok = [(st, p) for st, p in sites if 60 < p < len(seq) - 160]
    if not ok:
        return float("nan"), float("nan")
    m10 = sum(has_box(upstream(seq, st, p, 15, 5), "TATAAT") for st, p in ok) / len(ok)
    m35 = sum(has_box(upstream(seq, st, p, 40, 28), "TTGACA") for st, p in ok) / len(ok)
    return m10, m35


def shifted(sites, by=100):
    return [(st, p + by if st == "+" else p - by) for st, p in sites]


def match(query, reference, window):
    """Number of query sites within `window` nt of a reference site (same strand)."""
    ref = collections.defaultdict(list)
    for st, p in reference:
        ref[st].append(p)
    for v in ref.values():
        v.sort()
    import bisect
    hits = 0
    for st, p in query:
        v = ref.get(st, [])
        i = bisect.bisect_left(v, p - window)
        if i < len(v) and v[i] <= p + window:
            hits += 1
    return hits


def coverage_track(path, seqid, mate, min_mapq, cache_dir):
    """Strand-specific coverage of one mate on one sequence (orthogonal data)."""
    key = "%s.cov.%s.m%d.q%d.pkl" % (os.path.basename(path), seqid, mate, min_mapq)
    cache = os.path.join(cache_dir, key) if cache_dir else None
    if cache and os.path.exists(cache):
        with open(cache, "rb") as handle:
            return pickle.load(handle)
    import pysam
    bam = pysam.AlignmentFile(path)
    length = bam.get_reference_length(seqid)
    cov = {"+": np.zeros(length + 2, dtype=np.float32), "-": np.zeros(length + 2, dtype=np.float32)}
    for read in bam.fetch(seqid):
        if read.is_secondary or read.is_supplementary or read.mapping_quality < min_mapq:
            continue
        if (mate == 1 and read.is_read2) or (mate == 2 and read.is_read1):
            continue
        strand = "-" if read.is_reverse else "+"
        for a, b in read.get_blocks():
            cov[strand][a + 1:b + 1] += 1
    if cache:
        os.makedirs(cache_dir, exist_ok=True)
        with open(cache, "wb") as handle:
            pickle.dump(cov, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return cov


def coverage_rise(sites, cov, min_ratio=1.5, width=50):
    """Share of sites where coverage over +1..+width exceeds −width..−1 by min_ratio."""
    n = hits = 0
    for st, p in sites:
        c = cov[st]
        if not width < p < len(c) - width - 1:
            continue
        if st == "+":
            down, up = c[p:p + width].mean(), c[p - width:p].mean()
        else:
            down, up = c[p - width + 1:p + 1].mean(), c[p + 1:p + width + 1].mean()
        n += 1
        hits += (down + 1) / (up + 1) >= min_ratio
    return hits / n if n else float("nan")


def on_chrom(tss_list, chrom):
    return [t for t in tss_list if t.seqid == chrom]


def as_sites(tss_list):
    return [(t.strand, t.pos) for t in tss_list]


# ---------------------------------------------------------------------------
# Runs
# ---------------------------------------------------------------------------

def evaluate_real(ds, data_dir, methods, downsample, cache_dir, seed=0):
    path = lambda p: os.path.join(data_dir, p)  # noqa: E731
    mapq, mate, chrom = ds.get("min_mapq", 0), ds.get("mate", 1), ds["chromosome"]
    enriched = [load_signal(path(p), mate, mapq, cache_dir) for p in ds["enriched"]]
    controls = [load_signal(path(p), mate, mapq, cache_dir) for p in ds.get("control", [])]
    seq = read_fasta(path(ds["genome"]))[chrom]
    genes = ct.read_genes(path(ds["annotation"]), tuple(ds.get("feature_types", ["CDS"])))
    ctx = {"exclude": ct.GeneIndex(ct.read_genes(path(ds["annotation"]), ("rRNA", "tRNA")))}
    reference = None
    if ds.get("reference"):
        reference = [(st, p) for s, st, p in ct.read_sites(path(ds["reference"]))
                     if s in (chrom, "gi|48994873|gb|U00096.2|")]
    cov = None
    if ds.get("orthogonal"):
        cov = coverage_track(path(ds["orthogonal"]), chrom, ds.get("orthogonal_mate", 2),
                             20, cache_dir)
    rows = []
    for name in methods:
        run = METHODS[name]
        # Controls are used only when there is one per enriched library.
        ctl = controls if len(controls) == len(enriched) else []
        tss = run(enriched, ctl, ctx)
        if tss is None:
            continue
        ct.classify(tss, genes)
        chrom_tss = on_chrom(tss, chrom)
        sites = as_sites(chrom_tss)
        out = collections.OrderedDict(dataset=ds["name"], method=name, n_tss=len(sites))
        counts = collections.Counter(c for t in chrom_tss for c in t.classes)
        for c in CLASSES:
            out["share_" + c] = counts[c] / len(sites) if sites else float("nan")
        m10, m35 = box_rates(sites, seq)
        out["box10_all"], out["box35_all"] = m10, m35
        for c in ("Primary", "Internal", "Antisense"):
            out["box10_" + c] = box_rates([(t.strand, t.pos) for t in chrom_tss if c in t.classes], seq)[0]
        out["box10_background"] = box_rates(shifted(sites), seq)[0]
        if reference:
            for w in (0, 5):
                out["ref_precision_%d" % w] = match(sites, reference, w) / len(sites) if sites else float("nan")
                out["ref_recall_%d" % w] = match(reference, sites, w) / len(reference)
        if cov is not None:
            out["ortho_rise"] = coverage_rise(sites, cov)
            out["ortho_rise_background"] = coverage_rise(shifted(sites), cov)
        if len(enriched) >= 2:
            support = []
            for i in range(len(enriched)):
                rest = enriched[:i] + enriched[i + 1:]
                rest_ctl = ctl[:i] + ctl[i + 1:] if ctl else []
                train = as_sites(on_chrom(run(rest, rest_ctl, ctx), chrom))
                held = as_sites(on_chrom(run([enriched[i]], [ctl[i]] if ctl else [], ctx), chrom))
                support.append(match(train, held, 2) / len(train) if train else float("nan"))
            out["loo_support"] = statistics.mean(support)
        rng = np.random.default_rng(seed)
        for frac in downsample:
            thin_e = [thin(s, frac, rng) for s in enriched]
            thin_c = [thin(s, frac, rng) for s in ctl]
            sub = as_sites(on_chrom(run(thin_e, thin_c, ctx), chrom))
            out["down%g_n" % frac] = len(sub)
            out["down%g_recovery" % frac] = match(sites, sub, 2) / len(sites) if sites else float("nan")
        rows.append(out)
        print("  %-10s %-14s %d TSS" % (ds["name"], name, len(sites)), file=sys.stderr)
    return rows


def evaluate_sim(methods, depths, seed=0, label="sim", controls=(False, True), **sim_args):
    rows = []
    for with_control in controls:
        kind = "%s_%s" % (label, "with_control" if with_control else "enriched_only")
        for depth in depths:
            sim = simulate.simulate(depth=depth, with_control=with_control, seed=seed,
                                    rrna_fraction=0.7 if with_control else 0.0, **sim_args)
            truth = [(st, p) for _, st, p in sim["truth"]]
            strengths = sorted(sim["strength"].values())
            cut1, cut2 = strengths[len(strengths) // 3], strengths[2 * len(strengths) // 3]
            tiers = {"weak": [], "medium": [], "strong": []}
            for (_, st, p), s in sim["strength"].items():
                tiers["weak" if s < cut1 else "medium" if s < cut2 else "strong"].append((st, p))
            processed = [(st, p) for _, st, p in sim["processed"]]
            ctx = {"exclude": ct.GeneIndex([ct.Gene(simulate.SEQID, sim["rrna_block"][0],
                                                    sim["rrna_block"][1], "+", "rRNA")])}
            for name in methods:
                tss = METHODS[name](sim["enriched"], sim["controls"], ctx)
                if tss is None:
                    continue
                sites = as_sites(tss)
                tp = match(sites, truth, 3)
                out = collections.OrderedDict(dataset="%s_depth%g" % (kind, depth), method=name,
                                              n_tss=len(sites))
                out["precision"] = tp / len(sites) if sites else float("nan")
                out["recall"] = match(truth, sites, 3) / len(truth)
                for tier, tsites in tiers.items():
                    out["recall_" + tier] = match(tsites, sites, 3) / len(tsites)
                out["processed_called"] = match(processed, sites, 3) / len(processed) if processed else float("nan")
                rows.append(out)
                print("  %-28s %-14s %d TSS" % (out["dataset"], name, len(sites)), file=sys.stderr)
    return rows


def write_rows(rows, prefix):
    os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(prefix + ".tsv", "w") as out:
        out.write("\t".join(keys) + "\n")
        for r in rows:
            out.write("\t".join(_fmt(r.get(k, "")) for k in keys) + "\n")
    # Markdown: one table, measures as rows, dataset/method as columns.
    cols = ["%s / %s" % (r["dataset"], r["method"]) for r in rows]
    with open(prefix + ".md", "w") as out:
        out.write("| measure | " + " | ".join(cols) + " |\n")
        out.write("|---|" + "---|" * len(cols) + "\n")
        for k in keys[2:]:
            out.write("| %s | %s |\n" % (k, " | ".join(_fmt(r.get(k, "")) for r in rows)))
    print("wrote %s.tsv and %s.md" % (prefix, prefix))


def _fmt(v):
    if isinstance(v, float):
        return "NA" if math.isnan(v) else ("%.3f" % v)
    return str(v)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="mode", required=True)
    real = sub.add_parser("real", help="Benchmark on the datasets in datasets.toml")
    real.add_argument("--data-dir", required=True)
    real.add_argument("--manifest", default=os.path.join(HERE, "datasets.toml"))
    real.add_argument("--datasets", help="Comma-separated dataset names (default: all)")
    real.add_argument("--methods", default=",".join(METHODS))
    real.add_argument("--downsample", default="0.5,0.25",
                      help="Comma-separated read fractions ('' to skip)")
    real.add_argument("--cache-dir", help="Where to cache 5'-end signals (default: <out>_cache)")
    real.add_argument("--out", required=True, help="Output prefix (.tsv and .md)")
    sim = sub.add_parser("sim", help="Benchmark on simulated libraries")
    sim.add_argument("--depths", default="1,0.25,0.05")
    sim.add_argument("--methods", default=",".join(METHODS))
    sim.add_argument("--length", type=int, default=1_000_000)
    sim.add_argument("--n-tss", type=int, default=1000)
    sim.add_argument("--seed", type=int, default=0)
    sim.add_argument("--profile", help="Calibrate the simulation from one dataset in the "
                     "manifest (needs --data-dir): genome length, TSS number and strength, "
                     "background rate of its first enriched library")
    sim.add_argument("--data-dir")
    sim.add_argument("--manifest", default=os.path.join(HERE, "datasets.toml"))
    sim.add_argument("--replicates", type=int, default=3)
    sim.add_argument("--position-bias", type=float, default=0.0,
                     help="Variance of a fixed per-position background multiplier "
                          "shared by all libraries (default 0; real data about 2)")
    sim.add_argument("--dispersion", type=float, default=0.1,
                     help="Variance of the per-library rate multipliers (default 0.1; "
                          "replicates of the public data give 0.01-0.05)")
    sim.add_argument("--cache-dir")
    sim.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    methods = [m for m in args.methods.split(",") if m]
    unknown = [m for m in methods if m not in METHODS]
    if unknown:
        parser.error("unknown method(s): %s (have %s)" % (", ".join(unknown), ", ".join(METHODS)))

    if args.mode == "sim":
        depths = [float(d) for d in args.depths.split(",")]
        params, label = dict(length=args.length, n_tss=args.n_tss), "sim"
        if args.profile:
            with open(args.manifest, "rb") as handle:
                ds = [d for d in tomllib.load(handle)["dataset"] if d["name"] == args.profile][0]
            p = lambda x: os.path.join(args.data_dir, x)  # noqa: E731
            sig = load_signal(p(ds["enriched"][0]), ds.get("mate", 1), ds.get("min_mapq", 0),
                              args.cache_dir or args.out + "_cache")
            chrom = ds["chromosome"]
            sig = {k: v for k, v in sig.items() if k[0] == chrom}
            called = method_ettwiller2016([sig], [], None)
            length = len(read_fasta(p(ds["genome"]))[chrom])
            params = simulate.profile_from_signal(sig, [t.key for t in called], length)
            label = "sim_" + args.profile
            print("profile %s: %s" % (args.profile, params), file=sys.stderr)
        if args.dispersion != 0.1:
            label += "_disp%g" % args.dispersion
        if args.position_bias:
            label += "_bias%g" % args.position_bias
        rows = evaluate_sim(methods, depths, seed=args.seed, label=label,
                            replicates=args.replicates, dispersion=args.dispersion,
                            position_bias=args.position_bias, **params)
    else:
        with open(args.manifest, "rb") as handle:
            datasets = tomllib.load(handle)["dataset"]
        if args.datasets:
            wanted = args.datasets.split(",")
            datasets = [d for d in datasets if d["name"] in wanted]
        down = [float(x) for x in args.downsample.split(",") if x]
        cache = args.cache_dir or args.out + "_cache"
        rows = []
        for ds in datasets:
            rows += evaluate_real(ds, args.data_dir, methods, down, cache)
    write_rows(rows, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
