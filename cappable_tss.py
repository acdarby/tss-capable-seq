#!/usr/bin/env python3
"""ANNOgesic-style TSS annotation for Cappable-seq data.

ANNOgesic (Yu et al. 2018, doi:10.1093/gigascience/giy096) predicts
transcription start sites (TSS) with TSSpredator, which was built for
dRNA-seq: it compares a TEX-treated library (TEX+) against an untreated one
(TEX-) and looks for coverage "steps".  Cappable-seq (Ettwiller et al. 2016,
doi:10.1186/s12864-016-2539-z) works differently: 5'-triphosphorylated
(primary) transcripts are capped with a biotinylated GTP analogue and pulled
down, so the 5' end of each read in the enriched library marks a TSS at
single-base resolution.  A TSS is therefore a position with many read
5' ends in the enriched library and relatively few in the non-enriched
control.

This module re-implements the TSS part of the ANNOgesic workflow for that
kind of data:

``call``
    Count read 5' ends per position and strand, normalise to reads per
    million (RPM), keep positions that pass an RPM and log2 enrichment
    threshold in enough replicates, merge nearby positions into one TSS and
    classify each TSS against a genome annotation as Primary, Secondary,
    Internal, Antisense or Orphan, using the categories of Sharma et al.
    2010 (doi:10.1038/nature08756) that ANNOgesic and TSSpredator use.
    Writes an ANNOgesic-style GFF3 file and a TSV table.

``annogesic-libs``
    Convert Cappable-seq alignments into wiggle files plus the
    ``--tex_notex_libs`` strings needed to run the original ANNOgesic
    ``tss_ps`` subcommand, with the enriched library standing in for TEX+
    and the control for TEX-.

Only the Python standard library is needed.  BAM input additionally needs
``pysam``; BED and wiggle input do not.
"""

import argparse
import bisect
import gzip
import math
import os
import sys
from collections import defaultdict

STRANDS = ("+", "-")
CLASS_ORDER = ("Primary", "Secondary", "Internal", "Antisense", "Orphan")


# ---------------------------------------------------------------------------
# Reading signal
# ---------------------------------------------------------------------------

def open_text(path):
    """Open a plain or gzip-compressed text file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def new_signal():
    """Signal container: {(seqid, strand): {position (1-based): value}}."""
    return defaultdict(lambda: defaultdict(float))


def five_prime_from_bed(path):
    """Count read 5' ends from a BED file of alignments (needs column 6)."""
    signal = new_signal()
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[5] not in STRANDS:
                raise ValueError(
                    "%s: BED input needs a strand in column 6: %r" % (path, line))
            seqid, start, end, strand = (
                fields[0], int(fields[1]), int(fields[2]), fields[5])
            # BED is 0-based half-open; positions here are 1-based.
            pos = start + 1 if strand == "+" else end
            signal[(seqid, strand)][pos] += 1
    return signal


def _usable_alignment(read, mate, min_mapq):
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    if read.is_qcfail or read.mapping_quality < min_mapq:
        return False
    if read.is_paired:
        return read.is_read1 if mate == 1 else read.is_read2
    return True


def _pysam():
    try:
        import pysam
    except ImportError:
        sys.exit("BAM input needs pysam (pip install pysam); "
                 "alternatively convert to BED with `bedtools bamtobed`.")
    return pysam


def five_prime_from_bam(path, mate=1, min_mapq=0):
    """Count read 5' ends from a BAM file.

    For paired-end data only the mate given by ``mate`` is counted; its 5' end
    must correspond to the RNA 5' end (read 1 in most Cappable-seq libraries,
    but check your library kit).  The strand is the strand that mate aligns to.
    """
    pysam = _pysam()
    signal = new_signal()
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not _usable_alignment(read, mate, min_mapq):
                continue
            if read.is_reverse:
                signal[(read.reference_name, "-")][read.reference_end] += 1
            else:
                signal[(read.reference_name, "+")][read.reference_start + 1] += 1
    return signal


def read_wig(path, strand, signal=None):
    """Read a variableStep/fixedStep wiggle file for one strand.

    Values are taken as absolute values because ANNOgesic/READemption wiggle
    files store reverse-strand signal as negative numbers.
    """
    signal = new_signal() if signal is None else signal
    chrom, step, span, pos = None, 1, 1, None
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            if line.startswith(("variableStep", "fixedStep")):
                params = dict(item.split("=", 1) for item in line.split()[1:])
                chrom = params["chrom"]
                span = int(params.get("span", 1))
                fixed = line.startswith("fixedStep")
                if fixed:
                    pos = int(params["start"])
                    step = int(params.get("step", 1))
                continue
            if chrom is None:
                raise ValueError("%s: data line before a wiggle header" % path)
            if fixed:
                value = abs(float(line))
                start = pos
                pos += step
            else:
                start_text, value_text = line.split()[:2]
                start, value = int(start_text), abs(float(value_text))
            if value:
                for offset in range(span):
                    signal[(chrom, strand)][start + offset] += value
    return signal


def load_library(spec, mate=1, min_mapq=0):
    """Load one library.

    ``spec`` is either a BAM/BED path, or ``forward.wig,reverse.wig`` for a
    pair of strand-specific wiggle files of 5'-end counts.
    """
    if "," in spec:
        forward, reverse = spec.split(",", 1)
        signal = read_wig(forward, "+")
        return read_wig(reverse, "-", signal)
    if spec.endswith(".bam"):
        return five_prime_from_bam(spec, mate=mate, min_mapq=min_mapq)
    if spec.endswith((".bed", ".bed.gz")):
        return five_prime_from_bed(spec)
    raise ValueError("Unrecognised library %r: give a .bam, .bed(.gz) or "
                     "'forward.wig,reverse.wig'" % spec)


def total_signal(signal):
    return sum(sum(positions.values()) for positions in signal.values())


# ---------------------------------------------------------------------------
# TSS calling
# ---------------------------------------------------------------------------

class TSS(object):
    """A called transcription start site."""

    def __init__(self, seqid, pos, strand, enriched_rpm, control_rpm,
                 log2fc, replicates):
        self.seqid = seqid
        self.pos = pos
        self.strand = strand
        self.enriched_rpm = enriched_rpm   # mean over replicates
        self.control_rpm = control_rpm     # mean over replicates
        self.log2fc = log2fc               # mean over replicates
        self.replicates = replicates       # replicates passing the cut-offs
        self.pvalue = float("nan")         # statistical calling only
        self.qvalue = float("nan")
        self.consistency_p = float("nan")  # replicate agreement (local calling)
        self.promoter_score = float("nan")  # genome-trained promoter score
        self.peak_ratio = float("nan")      # RPM / strongest call upstream (mark_shoulders)
        self.classes = set()
        self.genes = defaultdict(list)     # class -> [(gene name, UTR length)]

    @property
    def key(self):
        return (self.seqid, self.strand, self.pos)


def log2_enrichment(enriched_rpm, control_rpm, pseudocount):
    top, bottom = enriched_rpm + pseudocount, control_rpm + pseudocount
    if top <= 0:
        return float("-inf")
    if bottom <= 0:
        return float("inf")
    return math.log2(top / bottom)


# Settings that reproduce the published Cappable-seq pipeline of Ettwiller
# et al. 2016 (scripts at https://github.com/Ettwiller/TSS):
#   bam2firstbasegtf.pl --cutoff 1.5   RRS >= 1.5, RRS rounded to 3 decimals
#   filter_tss.pl --cutoff 0           log2(RRS_cap / RRS_control) >= 0, no
#                                      pseudocount; a position with no control
#                                      reads gets the RRS of one control read
#   cluster_tss.pl --cutoff 5          bedtools merge -s -d 5 on 1-nt
#                                      features, i.e. starts <= 6 nt apart are
#                                      chained, the one with most reads kept
PRESETS = {
    "default": dict(min_rpm=1.0, min_log2fc=1.0, pseudocount=1.0,
                    control_floor=0, rpm_decimals=None, cluster_distance=5),
    "ettwiller2016": dict(min_rpm=1.5, min_log2fc=0.0, pseudocount=0.0,
                          control_floor=1, rpm_decimals=3, cluster_distance=6),
}


def call_positions(enriched, controls, min_rpm=1.0, min_log2fc=1.0,
                   min_replicates=1, pseudocount=1.0, control_floor=0,
                   rpm_decimals=None):
    """Return candidate TSS positions before clustering.

    ``enriched`` is a list of signals, one per replicate; ``controls`` is a
    list of the same length (paired by replicate) or empty if the experiment
    has no non-enriched control, in which case only ``min_rpm`` applies.
    ``control_floor`` is the minimum control read count used in the ratio
    (Ettwiller et al. use 1), and ``rpm_decimals`` rounds RPM values before
    they are compared with ``min_rpm`` (Ettwiller et al. use 3).
    """
    if controls and len(controls) != len(enriched):
        raise ValueError("Give one control per enriched replicate "
                         "(%d enriched, %d control)" % (len(enriched), len(controls)))
    if min_replicates > len(enriched):
        raise ValueError("min_replicates (%d) is larger than the number of "
                         "replicates (%d)" % (min_replicates, len(enriched)))

    scales = [1e6 / (total_signal(lib) or 1.0) for lib in enriched]
    control_scales = [1e6 / (total_signal(lib) or 1.0) for lib in controls]

    keys = set()
    for lib in enriched:
        for (seqid, strand), positions in lib.items():
            keys.update((seqid, strand, pos) for pos in positions)

    calls = []
    for seqid, strand, pos in keys:
        e_values, c_values, fcs, passed = [], [], [], 0
        for i, lib in enumerate(enriched):
            e_rpm = lib.get((seqid, strand), {}).get(pos, 0.0) * scales[i]
            if rpm_decimals is not None:
                e_rpm = round(e_rpm, rpm_decimals)
            if controls:
                c_count = controls[i].get((seqid, strand), {}).get(pos, 0.0)
                c_rpm = max(c_count, control_floor) * control_scales[i]
                if rpm_decimals is not None:
                    c_rpm = round(c_rpm, rpm_decimals)
                fc = log2_enrichment(e_rpm, c_rpm, pseudocount)
                ok = e_rpm >= min_rpm and fc >= min_log2fc
            else:
                c_rpm, fc = 0.0, float("nan")
                ok = e_rpm >= min_rpm
            e_values.append(e_rpm)
            c_values.append(c_rpm)
            fcs.append(fc)
            passed += ok
        if passed >= min_replicates:
            n = float(len(enriched))
            calls.append(TSS(seqid, pos, strand, sum(e_values) / n,
                             sum(c_values) / n, sum(fcs) / n, passed))
    calls.sort(key=lambda t: t.key)
    return calls


def cluster_tss(calls, distance=5):
    """Merge calls on the same strand no more than ``distance`` nt apart.

    Within each cluster the position with the highest enriched RPM is kept,
    as many bacterial promoters give a few adjacent start positions.
    ``distance=0`` disables clustering.
    """
    if distance <= 0:
        return list(calls)
    merged, cluster = [], []
    for tss in sorted(calls, key=lambda t: t.key):
        if cluster and (tss.seqid, tss.strand) == (cluster[-1].seqid, cluster[-1].strand) \
                and tss.pos - cluster[-1].pos <= distance:
            cluster.append(tss)
            continue
        if cluster:
            merged.append(max(cluster, key=lambda t: t.enriched_rpm))
        cluster = [tss]
    if cluster:
        merged.append(max(cluster, key=lambda t: t.enriched_rpm))
    return merged


def mark_shoulders(calls, window=30):
    """Set tss.peak_ratio for each call: its enriched RPM divided by that of
    the strongest call on the same strand within `window` nt upstream (in the
    direction of transcription), or 1.0 if there is none or it is weaker.

    A low ratio means the call sits in the downstream shoulder of a stronger
    TSS: reads from just after the start, or degradation of that transcript.
    """
    by_strand = defaultdict(list)
    for tss in calls:
        by_strand[(tss.seqid, tss.strand)].append(tss)
    for (_, strand), group in by_strand.items():
        group.sort(key=lambda t: t.pos)
        positions = [t.pos for t in group]
        for i, tss in enumerate(group):
            if strand == "+":
                lo, hi = bisect.bisect_left(positions, tss.pos - window), i
            else:
                lo, hi = i + 1, bisect.bisect_right(positions, tss.pos + window)
            top = max((group[j].enriched_rpm for j in range(lo, hi)), default=0.0)
            tss.peak_ratio = (min(1.0, tss.enriched_rpm / top)
                              if top > tss.enriched_rpm else 1.0)
    return calls


def remove_shoulders(calls, window=30, ratio=0.1):
    """Calls whose peak_ratio (mark_shoulders) is at least `ratio`."""
    mark_shoulders(calls, window)
    return [t for t in calls if t.peak_ratio >= ratio]


# ---------------------------------------------------------------------------
# Statistical calling
# ---------------------------------------------------------------------------

def _beta_cf(a, b, x, max_iter=10000, eps=1e-15):
    """Continued fraction for the incomplete beta function (modified Lentz)."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betainc(a, b, x):
    """Regularized incomplete beta function I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_front = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                 + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(log_front) * _beta_cf(a, b, x) / a
    return 1.0 - math.exp(log_front) * _beta_cf(b, a, 1.0 - x) / b


def binom_sf(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    return betainc(k, n - k + 1, p)


def betabinom_sf(k, n, p, dispersion):
    """P(X >= k) for X ~ BetaBinomial(n, a, b) with mean p.

    The beta-binomial adds variation of the success probability between
    libraries to binomial sampling. With independent gamma-distributed rate
    multipliers of variance `dispersion` (d) in each of two libraries, the log
    odds of the split of reads between them has variance about 2d whatever the
    library sizes. The Beta is matched to that: its logit has variance
    trigamma(a) + trigamma(b), about 1/a + 1/b, so a = 1 / (2d(1 - p)) and
    b = 1 / (2dp), which keeps the mean at p and is exact (Beta(1/d, 1/d)) for
    p = 0.5. `dispersion` 0 gives the binomial test.
    """
    if dispersion <= 0:
        return binom_sf(k, n, p)
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    a, b = 1.0 / (2.0 * dispersion * (1.0 - p)), 1.0 / (2.0 * dispersion * p)

    def log_pmf(x):
        return (math.lgamma(n + 1) - math.lgamma(x + 1) - math.lgamma(n - x + 1)
                + math.lgamma(x + a) + math.lgamma(n - x + b) - math.lgamma(n + a + b)
                - math.lgamma(a) - math.lgamma(b) + math.lgamma(a + b))

    def ratio(x):                      # pmf(x + 1) / pmf(x)
        return (n - x) / (x + 1.0) * (x + a) / (n - x - 1.0 + b)

    mean = n * p
    if k > mean:                       # sum the upper tail from k upwards
        term = math.exp(log_pmf(k))
        total_p, x = term, k
        while x < n and term > total_p * 1e-14:
            term *= ratio(x)
            x += 1
            total_p += term
        return min(1.0, total_p)
    term = math.exp(log_pmf(k - 1))    # 1 - lower tail, summed downwards
    lower, x = term, k - 1
    while x > 0 and term > lower * 1e-14:
        term /= ratio(x - 1)
        x -= 1
        lower += term
    return max(0.0, 1.0 - lower)


def estimate_dispersion(libraries, min_count=50):
    """Between-library dispersion from replicate libraries.

    For each pair of replicates, the variance of log(count share) differences at
    positions with at least `min_count` reads in both is about 2 * dispersion
    plus a Poisson part (1/x + 1/y); the median over pairs is returned, or
    None with fewer than two libraries or too few positions.
    """
    if len(libraries) < 2:
        return None
    sizes = [total_signal(lib) or 1.0 for lib in libraries]
    estimates = []
    for i in range(len(libraries)):
        for j in range(i + 1, len(libraries)):
            diffs, poisson = [], []
            for key, da in libraries[i].items():
                db = libraries[j].get(key, {})
                for pos, x in da.items():
                    y = db.get(pos, 0.0)
                    if x >= min_count and y >= min_count:
                        diffs.append(math.log(x / sizes[i]) - math.log(y / sizes[j]))
                        poisson.append(1.0 / x + 1.0 / y)
            if len(diffs) >= 50:
                mean = sum(diffs) / len(diffs)
                var = sum((d - mean) ** 2 for d in diffs) / (len(diffs) - 1)
                estimates.append(max(0.0, (var - sum(poisson) / len(poisson)) / 2.0))
    if not estimates:
        return None
    estimates.sort()
    return estimates[len(estimates) // 2]


DEFAULT_DISPERSION = 0.05


def chi2_sf_even(x, dof):
    """Upper tail of a chi-squared distribution with an even number of degrees."""
    half, term, total = x / 2.0, 1.0, 1.0
    for i in range(1, dof // 2):
        term *= half / i
        total += term
    return min(1.0, math.exp(-half) * total)


def chi2_sf(x, dof):
    """Upper tail of a chi-squared distribution with an integer number of degrees."""
    if x <= 0:
        return 1.0
    if dof % 2 == 0:
        return chi2_sf_even(x, dof)
    half = x / 2.0
    total = math.erfc(math.sqrt(half))          # one degree of freedom
    term = math.sqrt(2.0 * x / math.pi) * math.exp(-half)
    for k in range(3, dof + 1, 2):              # sf(k) = sf(k - 2) + term
        total += term
        term *= x / k
    return min(1.0, total)


def consistency_pvalue(counts, sizes, dispersion):
    """Do replicate counts at a site agree, given library sizes?

    Under agreement replicate i has mean E_i = n * N_i / sum(N) (n = pooled
    count) and variance E_i + dispersion * E_i^2, the between-library
    variation estimated with estimate_dispersion. Returns the chi-squared
    p-value of sum((x_i - E_i)^2 / var_i) with (replicates - 1) degrees of
    freedom; small values mark sites carried by one library (e.g. a stack of
    PCR duplicates). NaN with fewer than two replicates.
    """
    if len(counts) < 2:
        return float("nan")
    n, total = sum(counts), sum(sizes)
    if n <= 0 or total <= 0:
        return 1.0
    stat = 0.0
    for x, size in zip(counts, sizes):
        e = n * size / total
        stat += (x - e) ** 2 / (e + dispersion * e * e)
    return chi2_sf(stat, len(counts) - 1)


def fisher_combine(pvalues):
    """Fisher's method: combined p-value of independent tests."""
    pvalues = [max(p, 1e-300) for p in pvalues]
    if len(pvalues) == 1:
        return pvalues[0]
    return chi2_sf_even(-2.0 * sum(math.log(p) for p in pvalues), 2 * len(pvalues))


def bh_qvalues(pvalues, n_tests=None):
    """Benjamini-Hochberg adjusted p-values (q-values), in input order.

    `n_tests` is the number of hypotheses if larger than len(pvalues): when
    only some positions are tested (say, those with at least 3 reads) but the
    p-values do not condition on that selection, the untested positions still
    count as tests (their p-values are all larger than the tested ones).
    """
    n = len(pvalues)
    total = max(n, n_tests or 0)
    order = sorted(range(n), key=lambda i: pvalues[i], reverse=True)
    q = [0.0] * n
    running = 1.0
    for rank, i in enumerate(order):
        running = min(running, pvalues[i] * total / (n - rank))
        q[i] = running
    return q


def library_size(signal, exclude=None):
    """Reads in a library, optionally leaving out reads that start in the
    regions of `exclude` (a GeneIndex, e.g. of rRNA and tRNA genes)."""
    if exclude is None:
        return total_signal(signal)
    size = 0.0
    for (seqid, _), positions in signal.items():
        for pos, value in positions.items():
            if not exclude.overlapping(seqid, pos, pos):
                size += value
    return size


def call_binomial(enriched, controls, fdr=0.05, min_fold=1.0, min_reads=3,
                  min_replicates=None, exclude=None, dispersion=0.0):
    """Call TSS by testing enriched against control counts at each position.

    For one replicate, with enriched count E, control count C and library
    sizes Ne and Nc, the null hypothesis is that the position is enriched by
    at most `min_fold`: under it E ~ Binomial(E + C, p0) with
    p0 = f*Ne / (f*Ne + Nc). The one-sided p-value is P(X >= E). Replicates
    are combined with Fisher's method, positions with at least `min_reads`
    reads in some replicate are tested, and positions with a Benjamini-Hochberg
    q-value <= `fdr` whose enrichment (point estimate) exceeds `min_fold` in at
    least `min_replicates` replicates (default: all) are called. Library sizes
    leave out reads in `exclude` (e.g. rRNA and tRNA genes), so an
    rRNA-dominated control does not make every other position look enriched.
    With `dispersion` > 0 the beta-binomial replaces the binomial, allowing for
    variation between libraries beyond sampling noise.
    """
    if not controls or len(controls) != len(enriched):
        raise ValueError("Statistical calling needs one control per enriched library")
    reps = len(enriched)
    min_replicates = reps if min_replicates is None else min_replicates
    e_size = [library_size(lib, exclude) or 1.0 for lib in enriched]
    c_size = [library_size(lib, exclude) or 1.0 for lib in controls]
    p0 = [min_fold * e / (min_fold * e + c) for e, c in zip(e_size, c_size)]

    keys = set()
    for lib in enriched:
        for (seqid, strand), positions in lib.items():
            keys.update((seqid, strand, pos) for pos, v in positions.items() if v >= min_reads)

    tested = []
    for seqid, strand, pos in keys:
        pvals, e_rpm, c_rpm, fcs, passed = [], [], [], [], 0
        for i in range(reps):
            e = enriched[i].get((seqid, strand), {}).get(pos, 0.0)
            c = controls[i].get((seqid, strand), {}).get(pos, 0.0)
            pvals.append(betabinom_sf(int(round(e)), int(round(e + c)), p0[i], dispersion))
            er, cr = e / e_size[i] * 1e6, c / c_size[i] * 1e6
            fc = math.log2((er + 1e-9) / (cr + 1e-9)) if er or cr else float("nan")
            e_rpm.append(er)
            c_rpm.append(cr)
            fcs.append(fc)
            passed += e * c_size[i] > min_fold * c * e_size[i]
        tested.append((seqid, strand, pos, fisher_combine(pvals), e_rpm, c_rpm, fcs, passed))

    qvals = bh_qvalues([t[3] for t in tested])
    calls = []
    for (seqid, strand, pos, p, e_rpm, c_rpm, fcs, passed), q in zip(tested, qvals):
        if q <= fdr and passed >= min_replicates:
            finite = [f for f in fcs if f == f and abs(f) != float("inf")]
            tss = TSS(seqid, pos, strand, sum(e_rpm) / reps, sum(c_rpm) / reps,
                      sum(finite) / len(finite) if finite else float("inf"), passed)
            tss.pvalue, tss.qvalue = p, q
            calls.append(tss)
    calls.sort(key=lambda t: t.key)
    return calls


def poisson_sf(k, mean):
    """P(X >= k) for X ~ Poisson(mean)."""
    if k <= 0:
        return 1.0
    if mean <= 0.0:
        return 0.0
    log_pmf = lambda x: x * math.log(mean) - mean - math.lgamma(x + 1)  # noqa: E731
    if k > mean:                       # upper tail, terms decrease from k
        term = math.exp(log_pmf(k))
        total, x = term, k
        while term > total * 1e-16:
            x += 1
            term *= mean / x
            total += term
        return min(1.0, total)
    term = math.exp(log_pmf(k - 1))    # 1 - lower tail, terms decrease from k - 1
    lower, x = term, k - 1
    while x > 0 and term > lower * 1e-16:
        term *= x / mean
        x -= 1
        lower += term
    return max(0.0, 1.0 - lower)


def nbinom_sf(k, mean, dispersion):
    """P(X >= k) for a negative binomial with this mean and variance
    mean + dispersion * mean^2 (dispersion 0 gives the Poisson)."""
    if dispersion <= 0:
        return poisson_sf(k, mean)
    if k <= 0:
        return 1.0
    if mean <= 0.0:
        return 0.0
    size = 1.0 / dispersion
    return betainc(k, size, mean / (size + mean))


class LocalBackground(object):
    """Reads per position around a site, from its neighbourhood on the same strand.

    The window is `window` nt either side, leaving out `gap` nt either side of
    the site (where its own reads spread). The `trim` busiest positions of the
    window are also left out, so a neighbouring TSS does not count as
    background. One pseudo-read keeps the rate above zero in empty windows.
    """

    def __init__(self, signal, window=100, gap=2, trim=5):
        self.window, self.gap, self.trim = window, gap, trim
        self.width = 2 * (window - gap) - trim
        self.index = {}
        for key, positions in signal.items():
            pos = sorted(positions)
            self.index[key] = (pos, [positions[p] for p in pos])

    def rate(self, seqid, strand, pos):
        positions, counts = self.index.get((seqid, strand), ((), ()))
        lo = bisect.bisect_left(positions, pos - self.window)
        hi = bisect.bisect_right(positions, pos + self.window)
        values = [counts[j] for j in range(lo, hi) if abs(positions[j] - pos) > self.gap]
        total = sum(values)
        if self.trim and values:
            total -= sum(sorted(values)[-self.trim:])
        return (total + 1.0) / self.width


# Dispersion of 5'-end counts around the local background. Read 5' ends are
# spread very unevenly, also between neighbouring positions and in the same
# way in every replicate, and no single negative binomial fits both the many
# empty positions and the heavy tail (maximum-likelihood estimates on public
# data run to 40-1000). On the benchmarks the ranking of TSS hardly changes
# between 0.5 and 10, and with replicates 2 stops where calls stop looking
# like promoters (benchmarks/README.md); q-values are nominal.
DEFAULT_LOCAL_DISPERSION = 2.0


def pool_signals(signals):
    """Sum of several signals."""
    pooled = new_signal()
    for signal in signals:
        for key, positions in signal.items():
            target = pooled[key]
            for pos, value in positions.items():
                target[pos] = target.get(pos, 0.0) + value
    return pooled


def call_local(enriched, controls=None, fdr=0.05, min_reads=3, min_replicates=1,
               window=100, gap=2, trim=5, dispersion=DEFAULT_LOCAL_DISPERSION,
               exclude=None):
    """Call TSS as 5'-end pile-ups above the local background (no control needed).

    Replicates are pooled and tested once: the uneven spread of read 5' ends
    is reproducible between replicates, so per-replicate tests would not be
    independent. A position is tested if the pooled count is at least
    `min_reads` and it has reads in at least `min_replicates` replicates
    (default 1: on public biological triplicates, weak sites seen in only
    some replicates were as promoter-like as those seen in all; missing
    reads there are sampling, not artefacts). Each call gets
    `consistency_p` (consistency_pvalue, with the between-replicate
    dispersion from estimate_dispersion or DEFAULT_DISPERSION). Under the null the count is negative binomial with the
    local background rate (LocalBackground) as mean and `dispersion`
    (see DEFAULT_LOCAL_DISPERSION). Positions with a
    Benjamini-Hochberg q-value <= `fdr` are called; every position of the
    genome (two strands) counts as a test. With controls, a call must
    also not be depleted in the enriched library (RPM at least that of the
    control, a missing control counting as one read, as in Ettwiller et al.
    2016), which removes processed ends.
    """
    reps = len(enriched)
    min_replicates = reps if min_replicates is None else min_replicates
    pooled = pool_signals(enriched)
    background = LocalBackground(pooled, window, gap, trim)
    between = (estimate_dispersion(enriched) or DEFAULT_DISPERSION) if reps > 1 else 0.0
    e_size = [library_size(lib, exclude) or 1.0 for lib in enriched]
    ctl = pool_signals(controls) if controls else None
    c_size = (library_size(ctl, exclude) or 1.0) if ctl is not None else None
    pooled_size = sum(e_size)

    tested = []
    for (seqid, strand), positions in pooled.items():
        for pos, x in positions.items():
            if x < min_reads:
                continue
            present = sum(1 for lib in enriched if lib.get((seqid, strand), {}).get(pos, 0.0) > 0)
            if present < min_replicates:
                continue
            rate = background.rate(seqid, strand, pos)
            tested.append((seqid, strand, pos, x, present, rate,
                           nbinom_sf(int(round(x)), rate, dispersion)))
    # Every position of the genome is a hypothesis: the p-values do not
    # condition on the read count that selected the tested positions.
    n_tests = 0
    for seqid in {key[0] for key in pooled}:
        n_tests += len(STRANDS) * max(max(pooled.get((seqid, s), {0: 0})) for s in STRANDS)
    qvals = bh_qvalues([t[6] for t in tested], n_tests)
    calls = []
    for (seqid, strand, pos, x, present, rate, p), q in zip(tested, qvals):
        if q > fdr:
            continue
        e_rpm = x / pooled_size * 1e6
        c_rpm = log2fc = float("nan")
        if ctl is not None:
            c = ctl.get((seqid, strand), {}).get(pos, 0.0)
            c_rpm = c / c_size * 1e6
            if e_rpm < max(c, 1.0) / c_size * 1e6:
                continue                 # depleted after enrichment: processed end
            log2fc = math.log2(e_rpm / c_rpm) if c else float("inf")
        tss = TSS(seqid, pos, strand, e_rpm, c_rpm, log2fc, present)
        tss.pvalue, tss.qvalue = p, q
        tss.background = rate
        tss.consistency_p = consistency_pvalue(
            [lib.get((seqid, strand), {}).get(pos, 0.0) for lib in enriched], e_size, between)
        calls.append(tss)
    calls.sort(key=lambda t: t.key)
    return calls


# ---------------------------------------------------------------------------
# Genome-trained promoter score (approach H)
# ---------------------------------------------------------------------------

PROMOTER_UP = 60          # bases -60..-1; downstream bases are left out so the
                          # model cannot learn 5'-ligation preferences
PROMOTER_BLOCK = 50000    # genome blocks for cross-validation
_BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    """{sequence name: upper-case sequence} from a (gzipped) FASTA file."""
    seqs, name, parts = {}, None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip().upper())
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def upstream_window(seq, strand, pos, length=PROMOTER_UP):
    """Bases -length..-1 of a 1-based TSS in transcript direction (None if not ACGT)."""
    if strand == "+":
        window = seq[max(0, pos - 1 - length):pos - 1]
    else:
        window = seq[pos:pos + length][::-1].translate(_COMPLEMENT)
    if len(window) != length or any(b not in _BASE for b in window):
        return None
    return window


def promoter_features(window):
    """Active binary features: the base at each position and each pair of
    neighbouring bases (so that, e.g., TA steps can be weighted)."""
    feats = [4 * i + _BASE[b] for i, b in enumerate(window)]
    offset = 4 * len(window)
    feats += [offset + 16 * i + 4 * _BASE[a] + _BASE[b]
              for i, (a, b) in enumerate(zip(window, window[1:]))]
    return feats


class LogisticModel(object):
    """L2-regularised logistic regression on binary (sparse) features.

    Full-batch gradient descent with Adam steps on the mean log-loss plus
    l2 / 2 * |w|^2 (the intercept is not penalised).
    """

    def __init__(self, n_features, l2=1e-3):
        self.w = [0.0] * n_features
        self.b = 0.0
        self.l2 = l2

    def decision(self, feats):
        return self.b + sum(self.w[j] for j in feats)

    def fit(self, rows, labels, iterations=150, rate=0.05):
        n, k = len(rows), len(self.w)
        m, v = [0.0] * (k + 1), [0.0] * (k + 1)
        beta1, beta2, eps = 0.9, 0.999, 1e-8
        for it in range(1, iterations + 1):
            grad = [0.0] * k
            gb = 0.0
            for feats, y in zip(rows, labels):
                z = self.decision(feats)
                p = 1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0
                err = p - y
                gb += err
                for j in feats:
                    grad[j] += err
            grad = [g / n + self.l2 * w for g, w in zip(grad, self.w)]
            grad.append(gb / n)
            corr1, corr2 = 1 - beta1 ** it, 1 - beta2 ** it
            for j, g in enumerate(grad):
                m[j] = beta1 * m[j] + (1 - beta1) * g
                v[j] = beta2 * v[j] + (1 - beta2) * g * g
                step = rate * (m[j] / corr1) / (math.sqrt(v[j] / corr2) + eps)
                if j < k:
                    self.w[j] -= step
                else:
                    self.b -= step
        return self


def _block_fold(pos, block=PROMOTER_BLOCK):
    return (pos // block) % 2


def matched_background(genome, sites, per_site=3, n_sample=100000, seed=0):
    """Random positions whose -60..-1 window has the same AT count as each site.

    AT-rich genomes contain promoter-like words by chance, and the regions
    around TSS are often more AT-rich than genes; matching composition makes
    the model learn the promoter rather than base composition.
    """
    import random
    rng = random.Random(seed)
    names = [n for n in genome if len(genome[n]) > 2 * PROMOTER_UP + 10]
    weights = [len(genome[n]) for n in names]
    bins = defaultdict(list)
    for _ in range(n_sample):
        name = rng.choices(names, weights)[0]
        pos = rng.randint(PROMOTER_UP + 2, len(genome[name]) - PROMOTER_UP - 2)
        strand = rng.choice(STRANDS)
        window = upstream_window(genome[name], strand, pos)
        if window:
            bins[sum(b in "AT" for b in window)].append((name, strand, pos))
    out = []
    for seqid, strand, pos in sites:
        window = upstream_window(genome.get(seqid, ""), strand, pos)
        if not window:
            continue
        at = sum(b in "AT" for b in window)
        for d in (0, 1, -1, 2, -2, 3, -3):
            pool = bins.get(at + d)
            if pool:
                out += [pool[rng.randrange(len(pool))] for _ in range(per_site)]
                break
    return out


def promoter_scores(genome, positives, sites, seed=0, max_positives=4000):
    """Genome-trained promoter score for `sites`, cross-validated by genome block.

    `positives` are confident TSS of this genome (seqid, strand, pos); negatives
    are composition-matched random positions (three per positive). The genome
    is cut into 50-kb blocks in two alternating halves; a model trained on one
    half scores the sites in the other, so no site is scored by a model that
    saw it. Returns {site: log-odds score} (sites without a clean window are
    left out).
    """
    import random
    rng = random.Random(seed)
    positives = list(positives)
    if len(positives) > max_positives:
        positives = rng.sample(positives, max_positives)
    negatives = matched_background(genome, positives, seed=seed)
    n_features = 4 * PROMOTER_UP + 16 * (PROMOTER_UP - 1)
    data = []
    for label, group in ((1, positives), (0, negatives)):
        for seqid, strand, pos in group:
            window = upstream_window(genome.get(seqid, ""), strand, pos)
            if window:
                data.append((_block_fold(pos), promoter_features(window), label))
    models = {}
    for fold in (0, 1):
        rows = [(f, y) for blk, f, y in data if blk != fold]
        models[fold] = LogisticModel(n_features).fit([f for f, _ in rows], [y for _, y in rows])
    scores = {}
    for seqid, strand, pos in sites:
        window = upstream_window(genome.get(seqid, ""), strand, pos)
        if window:
            scores[(seqid, strand, pos)] = models[_block_fold(pos)].decision(
                promoter_features(window))
    return scores


def add_promoter_scores(calls, genome, strong_share=1 / 3.0, min_positives=200, stream=None):
    """Score each call (tss.promoter_score) with a model trained on the
    strongest `strong_share` of the calls themselves (by enriched RPM)."""
    ranked = sorted(calls, key=lambda t: -t.enriched_rpm)
    positives = [t.key for t in ranked[:int(len(ranked) * strong_share)]
                 if t.seqid in genome]
    if len(positives) < min_positives:
        if stream is not None:
            stream.write("Promoter score skipped: %d confident TSS, need %d\n"
                         % (len(positives), min_positives))
        return False
    scores = promoter_scores(genome, positives, [t.key for t in calls])
    for t in calls:
        t.promoter_score = scores.get(t.key, float("nan"))
    if stream is not None:
        stream.write("Promoter score: trained on %d of this genome's strongest TSS "
                     "(cross-validated by 50-kb blocks)\n" % min(len(positives), 4000))
    return True


# ---------------------------------------------------------------------------
# Annotation and classification
# ---------------------------------------------------------------------------

class Gene(object):
    def __init__(self, seqid, start, end, strand, name):
        self.seqid, self.start, self.end = seqid, start, end
        self.strand, self.name = strand, name


def parse_attributes(text):
    attributes = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key.strip()] = value.strip()
    return attributes


def read_genes(path, feature_types=("gene",)):
    """Read features of the given types from a GFF3 file."""
    genes = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in feature_types:
                continue
            if fields[6] not in STRANDS:
                continue
            attrs = parse_attributes(fields[8])
            name = (attrs.get("locus_tag") or attrs.get("Name")
                    or attrs.get("ID") or "%s:%s-%s" % (fields[0], fields[3], fields[4]))
            genes.append(Gene(fields[0], int(fields[3]), int(fields[4]),
                              fields[6], name))
    return genes


class GeneIndex(object):
    """Find genes overlapping an interval on one sequence."""

    def __init__(self, genes):
        self.by_seq = defaultdict(list)
        for gene in genes:
            self.by_seq[gene.seqid].append(gene)
        self.starts, self.max_len = {}, {}
        for seqid, items in self.by_seq.items():
            items.sort(key=lambda g: g.start)
            self.starts[seqid] = [g.start for g in items]
            self.max_len[seqid] = max(g.end - g.start for g in items)

    def overlapping(self, seqid, lo, hi):
        """Genes on ``seqid`` that overlap [lo, hi] (1-based, inclusive)."""
        if seqid not in self.by_seq:
            return []
        starts, items = self.starts[seqid], self.by_seq[seqid]
        first = bisect.bisect_left(starts, lo - self.max_len[seqid])
        last = bisect.bisect_right(starts, hi)
        return [g for g in items[first:last] if g.end >= lo and g.start <= hi]


def classify(tss_list, genes, utr_length=300, antisense_distance=100):
    """Assign TSS classes following Sharma et al. 2010.

    * Primary   - within ``utr_length`` nt upstream of a gene on the same
                  strand, and the strongest such TSS for that gene.
    * Secondary - within ``utr_length`` nt upstream of a gene but not the
                  strongest.
    * Internal  - inside a gene, on the same strand.
    * Antisense - inside a gene or within ``antisense_distance`` nt of it, on
                  the opposite strand.
    * Orphan    - none of the above.

    A TSS can belong to several classes (reported joined by '&', as in
    ANNOgesic).  Gene coordinates are taken from the start codon end of the
    feature (start for '+', end for '-'), so pass CDS or gene features rather
    than transcripts that already include a UTR.
    """
    index = GeneIndex(genes)
    upstream = defaultdict(list)    # gene -> [(tss, utr length)]

    for tss in tss_list:
        lo, hi = tss.pos - utr_length, tss.pos + utr_length
        nearby = index.overlapping(tss.seqid, lo - antisense_distance,
                                   hi + antisense_distance)
        for gene in nearby:
            if gene.strand == tss.strand:
                if tss.strand == "+":
                    utr = gene.start - tss.pos
                    inside = gene.start < tss.pos <= gene.end
                else:
                    utr = tss.pos - gene.end
                    inside = gene.start <= tss.pos < gene.end
                if 0 <= utr <= utr_length:
                    upstream[gene].append((tss, utr))
                elif inside:
                    tss.classes.add("Internal")
                    tss.genes["Internal"].append((gene.name, None))
            elif (gene.start - antisense_distance <= tss.pos
                    <= gene.end + antisense_distance):
                tss.classes.add("Antisense")
                tss.genes["Antisense"].append((gene.name, None))

    for gene, candidates in upstream.items():
        best = max(candidates, key=lambda item: (item[0].enriched_rpm, -item[1]))
        for tss, utr in candidates:
            label = "Primary" if tss is best[0] else "Secondary"
            tss.classes.add(label)
            tss.genes[label].append((gene.name, utr))

    for tss in tss_list:
        if not tss.classes:
            tss.classes.add("Orphan")
    return tss_list


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _class_text(tss):
    # "NA" when no annotation was given, so nothing was classified.
    return "&".join(c for c in CLASS_ORDER if c in tss.classes) or "NA"


def _gene_text(tss):
    names = []
    for label in CLASS_ORDER:
        for name, _ in tss.genes.get(label, []):
            if name not in names:
                names.append(name)
    if names:
        return ",".join(names)
    return "orphan" if "Orphan" in tss.classes else "NA"


def _utr_text(tss):
    parts = []
    for label in CLASS_ORDER:
        for _, utr in tss.genes.get(label, []):
            parts.append("%s_%s" % (label, "NA" if utr is None else utr))
    if parts:
        return "&".join(parts)
    return "Orphan_NA" if "Orphan" in tss.classes else "NA"


def _fmt_p(value):
    return "NA" if value != value else "%.3g" % value


def _fmt(value):
    return "NA" if value != value else "%.3f" % value   # NaN check


def write_gff(tss_list, path, source="ANNOgesic-Cappable"):
    with open(path, "w") as out:
        out.write("##gff-version 3\n")
        for i, tss in enumerate(tss_list):
            attrs = [
                "ID=%s_tss%d" % (tss.seqid, i),
                "Name=TSS:%d_%s" % (tss.pos, "f" if tss.strand == "+" else "r"),
                "type=%s" % _class_text(tss),
                "associated_gene=%s" % _gene_text(tss),
                "UTR_length=%s" % _utr_text(tss),
                "enriched_rpm=%s" % _fmt(tss.enriched_rpm),
                "control_rpm=%s" % _fmt(tss.control_rpm),
                "log2_enrichment=%s" % _fmt(tss.log2fc),
                "replicates=%d" % tss.replicates,
            ]
            if tss.qvalue == tss.qvalue:
                attrs += ["pvalue=%.3g" % tss.pvalue, "qvalue=%.3g" % tss.qvalue]
            if tss.consistency_p == tss.consistency_p:
                attrs.append("consistency_p=%.3g" % tss.consistency_p)
            if tss.promoter_score == tss.promoter_score:
                attrs.append("promoter_score=%.2f" % tss.promoter_score)
            if tss.peak_ratio == tss.peak_ratio:
                attrs.append("peak_ratio=%.3g" % tss.peak_ratio)
            out.write("\t".join([tss.seqid, source, "TSS", str(tss.pos),
                                 str(tss.pos), _fmt(tss.enriched_rpm),
                                 tss.strand, ".", ";".join(attrs)]) + "\n")


def write_table(tss_list, path):
    header = ["seqid", "position", "strand", "type", "associated_gene",
              "UTR_length", "enriched_rpm", "control_rpm", "log2_enrichment",
              "replicates", "pvalue", "qvalue", "consistency_p", "promoter_score",
              "peak_ratio"]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for tss in tss_list:
            out.write("\t".join([
                tss.seqid, str(tss.pos), tss.strand, _class_text(tss),
                _gene_text(tss), _utr_text(tss), _fmt(tss.enriched_rpm),
                _fmt(tss.control_rpm), _fmt(tss.log2fc), str(tss.replicates),
                _fmt_p(tss.pvalue), _fmt_p(tss.qvalue), _fmt_p(tss.consistency_p),
                _fmt(tss.promoter_score), _fmt(tss.peak_ratio),
            ]) + "\n")


def summarise(tss_list):
    counts = defaultdict(int)
    for tss in tss_list:
        for label in tss.classes:
            counts[label] += 1
    lines = ["TSS called: %d" % len(tss_list)]
    lines += ["  %-10s %d" % (label, counts[label]) for label in CLASS_ORDER]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# ANNOgesic export
# ---------------------------------------------------------------------------

def coverage_from_bed(path):
    """Per-base read coverage from a BED file, by strand."""
    diffs = defaultdict(lambda: defaultdict(int))
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[5] not in STRANDS:
                raise ValueError("%s: BED input needs a strand in column 6" % path)
            key = (fields[0], fields[5])
            diffs[key][int(fields[1]) + 1] += 1
            diffs[key][int(fields[2]) + 1] -= 1
    return _sweep(diffs)


def coverage_from_bam(path, mate=1, min_mapq=0):
    """Per-base read coverage from a BAM file, by strand of the chosen mate."""
    pysam = _pysam()
    diffs = defaultdict(lambda: defaultdict(int))
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not _usable_alignment(read, mate, min_mapq):
                continue
            key = (read.reference_name, "-" if read.is_reverse else "+")
            for start, end in read.get_blocks():
                diffs[key][start + 1] += 1
                diffs[key][end + 1] -= 1
    return _sweep(diffs)


def _sweep(diffs):
    signal = new_signal()
    for key, changes in diffs.items():
        depth, previous = 0, None
        for pos in sorted(changes):
            if depth and previous is not None:
                for p in range(previous, pos):
                    signal[key][p] = depth
            depth += changes[pos]
            previous = pos
    return signal


def write_wig(signal, strand, path, name, negative=False):
    sign = -1 if negative else 1
    with open(path, "w") as out:
        out.write('track type=wiggle_0 name="%s"\n' % name)
        for (seqid, s), positions in sorted(signal.items()):
            if s != strand or not positions:
                continue
            out.write("variableStep chrom=%s span=1\n" % seqid)
            for pos in sorted(positions):
                value = positions[pos] * sign
                text = ("%d" % value) if float(value).is_integer() else ("%.4f" % value)
                out.write("%d %s\n" % (pos, text))


def export_annogesic(enriched, controls, outdir, signal="coverage", mate=1,
                     min_mapq=0, condition=1, negative_reverse=True):
    """Write wiggle files and ANNOgesic ``--tex_notex_libs`` strings.

    Enriched libraries are labelled ``tex`` and controls ``notex``, so the
    TEX+/TEX- comparison in ANNOgesic/TSSpredator becomes an
    enriched/control comparison.  Returns the list of library strings.
    """
    if len(controls) != len(enriched):
        raise ValueError("ANNOgesic needs one control per enriched library")
    os.makedirs(outdir, exist_ok=True)
    libs = []
    replicate_ids = "abcdefghijklmnopqrstuvwxyz"
    for kind, paths in (("tex", enriched), ("notex", controls)):
        for i, path in enumerate(paths):
            if signal == "coverage":
                data = (coverage_from_bam(path, mate, min_mapq)
                        if path.endswith(".bam") else coverage_from_bed(path))
            else:
                data = load_library(path, mate=mate, min_mapq=min_mapq)
            stem = os.path.basename(path).split(".")[0]
            for strand, suffix in (("+", "forward"), ("-", "reverse")):
                wig = os.path.join(outdir, "%s_%s_%s.wig" % (stem, kind, suffix))
                write_wig(data, strand, wig, "%s_%s" % (stem, suffix),
                          negative=(strand == "-" and negative_reverse))
                libs.append("%s:%s:%d:%s:%s" % (
                    os.path.abspath(wig), kind, condition, replicate_ids[i], strand))
    return libs


# ---------------------------------------------------------------------------
# Comparison with a reference TSS set
# ---------------------------------------------------------------------------

STRAND_NAMES = {"+": "+", "-": "-", "forward": "+", "reverse": "-",
                "f": "+", "r": "-", "1": "+", "-1": "-", "+1": "+"}


def read_sites(path, columns=None, rename=None):
    """Read TSS positions as a set of (seqid, strand, position).

    GFF/GTF files are read from columns 1, 4 and 7.  For other tables give
    ``columns`` as 1-based (seqid, position, strand) column numbers; the
    separator is a tab, or a comma for .csv files, and lines whose position
    column is not a number (headers) are skipped.  ``rename`` maps sequence
    names in the file to the names used in the other set.
    """
    rename = rename or {}
    is_gff = columns is None
    if is_gff:
        columns = (1, 4, 7)
    sep = "," if path.endswith((".csv", ".csv.gz")) else "\t"
    sites = set()
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = [f.strip().strip('"') for f in line.rstrip("\n").split(sep)]
            try:
                seqid = fields[columns[0] - 1]
                pos = int(float(fields[columns[1] - 1]))
                strand = STRAND_NAMES[fields[columns[2] - 1].lower()]
            except (IndexError, ValueError, KeyError):
                if is_gff:
                    raise ValueError("%s: not a GFF/GTF line: %r" % (path, line))
                continue
            sites.add((rename.get(seqid, seqid), strand, pos))
    return sites


def _nearest(sorted_positions, pos):
    """Nearest position in a sorted list, or None."""
    i = bisect.bisect_left(sorted_positions, pos)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(sorted_positions):
            if best is None or abs(sorted_positions[j] - pos) < abs(best - pos):
                best = sorted_positions[j]
    return best


def compare_sites(query, reference, windows=(0, 5, 25)):
    """Compare two TSS sets.

    For every site the nearest site of the other set on the same sequence and
    strand is found.  Offsets are query minus reference in the direction of
    transcription (positive: the query TSS is downstream).  Returns a summary
    dict and per-query rows (seqid, strand, position, offset or None).
    """
    def index(sites):
        by_key = defaultdict(list)
        for seqid, strand, pos in sites:
            by_key[(seqid, strand)].append(pos)
        for positions in by_key.values():
            positions.sort()
        return by_key

    ref_index, query_index = index(reference), index(query)

    def offsets(sites, other_index):
        result = []
        for seqid, strand, pos in sorted(sites):
            near = _nearest(other_index.get((seqid, strand), []), pos)
            if near is None:
                result.append((seqid, strand, pos, None))
            else:
                diff = pos - near
                result.append((seqid, strand, pos, diff if strand == "+" else -diff))
        return result

    query_rows = offsets(query, ref_index)
    ref_rows = offsets(reference, query_index)
    summary = {"query": len(query), "reference": len(reference),
               "shared_seqids": sorted({s for s, _, _ in query} & {s for s, _, _ in reference}),
               "windows": []}
    for w in windows:
        q_hit = sum(1 for r in query_rows if r[3] is not None and abs(r[3]) <= w)
        r_hit = sum(1 for r in ref_rows if r[3] is not None and abs(r[3]) <= w)
        summary["windows"].append({
            "window": w,
            "query_matched": q_hit,
            "reference_recovered": r_hit,
            "precision": q_hit / float(len(query)) if query else float("nan"),
            "recall": r_hit / float(len(reference)) if reference else float("nan"),
        })
    return summary, query_rows


def format_comparison(summary):
    lines = ["query TSS:     %d" % summary["query"],
             "reference TSS: %d" % summary["reference"]]
    if not summary["shared_seqids"]:
        lines.append("WARNING: no sequence names in common; use --rename")
    lines.append("window  query_in_ref  ref_recovered  precision  recall")
    for w in summary["windows"]:
        lines.append("%6d  %12d  %13d  %9.3f  %6.3f" % (
            w["window"], w["query_matched"], w["reference_recovered"],
            w["precision"], w["recall"]))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Library quality checks
# ---------------------------------------------------------------------------

# Minimum share of reads at sharp peaks for a library to look like 5'-end
# data. Measured on public data (benchmarks/README.md): 5'-end enriched
# libraries 77-97%, the E. coli 5'-end control 42%, fragment (random-primed)
# libraries 22-26%.
QC_MIN_PEAK_SHARE = {"enriched": 0.50, "control": 0.30}
# Minimum share of read 5' ends inside genes that lie on the gene's strand.
# 5'-end libraries: 66-96%; the antisense mate of a stranded library: 2%.
QC_MIN_SENSE_SHARE = 0.30


def peak_share(signal, window=20, fold=5.0, min_reads=2):
    """Share of read 5' ends at sharp peaks.

    A position is a peak if it has at least `min_reads` reads and at least
    `fold` times the mean of the surrounding +-`window` nt (itself excluded).
    In 5'-end libraries most reads start at TSS or processed ends; in fragment
    libraries read starts are spread along transcripts.
    """
    total = peaks = 0.0
    for positions in signal.values():
        if not positions:
            continue
        ordered = sorted(positions)
        values = [positions[p] for p in ordered]
        prefix = [0.0]
        for v in values:
            prefix.append(prefix[-1] + v)
        lo = hi = 0
        for i, pos in enumerate(ordered):
            while ordered[lo] < pos - window:
                lo += 1
            while hi < len(ordered) and ordered[hi] <= pos + window:
                hi += 1
            neighbours = (prefix[hi] - prefix[lo] - values[i]) / (2.0 * window)
            total += values[i]
            if values[i] >= min_reads and values[i] >= fold * max(neighbours, 1.0 / (2 * window)):
                peaks += values[i]
    return peaks / total if total else float("nan")


def sense_share(signal, genes):
    """Share of read 5' ends inside genes that are on the gene's strand."""
    index = GeneIndex(genes)
    sense = anti = 0.0
    for (seqid, strand), positions in signal.items():
        for pos, value in positions.items():
            hits = index.overlapping(seqid, pos, pos)
            if not hits:
                continue
            if any(g.strand == strand for g in hits):
                sense += value
            else:
                anti += value
    return sense / (sense + anti) if sense + anti else float("nan")


def library_qc(signal, genes=None):
    """Summary statistics used to check that a library is 5'-end data."""
    reads = total_signal(signal)
    return {"reads": reads,
            "positions": sum(len(v) for v in signal.values()),
            "peak_share": peak_share(signal),
            "sense_share": sense_share(signal, genes) if genes else float("nan")}


def qc_warnings(qc, role):
    """Plain-language warnings for one library ('enriched' or 'control')."""
    warnings = []
    if qc["reads"] == 0:
        return ["no reads were counted (check --mate, --min-mapq and sequence names)"]
    if qc["peak_share"] == qc["peak_share"] and qc["peak_share"] < QC_MIN_PEAK_SHARE[role]:
        warnings.append(
            "only %.0f%% of reads start at sharp peaks (expected at least %.0f%% for a "
            "5'-end %s library); this looks like fragment (random-primed) RNA-seq, where "
            "read starts do not mark TSS" % (100 * qc["peak_share"],
                                             100 * QC_MIN_PEAK_SHARE[role], role))
    if qc["sense_share"] == qc["sense_share"] and qc["sense_share"] < QC_MIN_SENSE_SHARE:
        warnings.append(
            "only %.0f%% of reads inside genes are on the gene's strand; the counted read "
            "is probably antisense (try the other mate with --mate, or check the library "
            "strandedness)" % (100 * qc["sense_share"]))
    return warnings


def run_qc(paths, signals, roles, genes=None, out_path=None, stream=sys.stderr):
    """Check every library, print a table and warnings, optionally write a TSV.

    Returns the number of libraries with warnings.
    """
    rows, flagged = [], 0
    for path, signal, role in zip(paths, signals, roles):
        qc = library_qc(signal, genes)
        warnings = qc_warnings(qc, role)
        flagged += bool(warnings)
        rows.append((path, role, qc, warnings))
    stream.write("Library check (5'-end data expected):\n")
    stream.write("  %-8s %12s %10s %8s %8s  %s\n" % ("role", "reads", "positions", "peaks",
                                                     "sense", "library"))
    for path, role, qc, warnings in rows:
        stream.write("  %-8s %12.0f %10d %7.1f%% %7s  %s\n" % (
            role, qc["reads"], qc["positions"], 100 * qc["peak_share"],
            "NA" if qc["sense_share"] != qc["sense_share"] else "%.1f%%" % (100 * qc["sense_share"]),
            path))
        for w in warnings:
            stream.write("  WARNING (%s): %s\n" % (os.path.basename(path), w))
    if out_path:
        with open(out_path, "w") as out:
            out.write("library\trole\treads\tpositions\tpeak_share\tsense_share\twarnings\n")
            for path, role, qc, warnings in rows:
                out.write("%s\t%s\t%.0f\t%d\t%s\t%s\t%s\n" % (
                    path, role, qc["reads"], qc["positions"], _fmt(qc["peak_share"]),
                    _fmt(qc["sense_share"]), "; ".join(warnings) or "none"))
    return flagged


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description="ANNOgesic-style TSS annotation for Cappable-seq data.")
    sub = parser.add_subparsers(dest="command")
    sub.required = True

    def add_inputs(p):
        p.add_argument("--enriched", nargs="+", required=True,
                       help="Cappable-seq enriched libraries, one per replicate "
                            "(.bam, .bed[.gz] or forward.wig,reverse.wig)")
        p.add_argument("--control", nargs="*", default=[],
                       help="Non-enriched control libraries, in the same "
                            "replicate order as --enriched")
        p.add_argument("--mate", type=int, choices=(1, 2), default=1,
                       help="For paired-end BAM: mate whose 5' end is the RNA "
                            "5' end (default 1)")
        p.add_argument("--min-mapq", type=int, default=0,
                       help="Minimum mapping quality for BAM reads (default 0)")

    call = sub.add_parser("call", help="Call and classify TSS")
    add_inputs(call)
    call.add_argument("--annotation", help="GFF3 genome annotation")
    call.add_argument("--shoulder-ratio", type=float, default=0.0,
                      help="Drop a TSS with less than this share of the RPM of a "
                           "stronger TSS on the same strand within --shoulder-window nt "
                           "upstream (its downstream shoulder); 0 keeps all (default 0)")
    call.add_argument("--shoulder-window", type=int, default=30,
                      help="Window for --shoulder-ratio, nt upstream (default 30)")
    call.add_argument("--genome",
                      help="Genome FASTA: adds a promoter score per TSS from a sequence "
                           "model trained on this genome's strongest TSS (-60..-1, "
                           "cross-validated by 50-kb genome blocks)")
    call.add_argument("--feature-types", default="gene",
                      help="Comma-separated GFF feature types used as genes "
                           "(default: gene; CDS is a common alternative)")
    call.add_argument("--method", choices=("cutoff", "binomial", "local"), default="cutoff",
                      help="cutoff: fixed RPM and enrichment cut-offs (the presets); "
                           "binomial: test enriched against control counts at each "
                           "position with a false-discovery rate (needs --control); "
                           "local: test 5'-end pile-ups against the local background "
                           "(no control needed; a control removes depleted sites). "
                           "Default: cutoff")
    call.add_argument("--fdr", type=float, default=0.05,
                      help="binomial, local: false-discovery rate (default 0.05)")
    call.add_argument("--min-fold", type=float, default=1.0,
                      help="binomial: enrichment over the control that the test must "
                           "exceed (default 1, i.e. any enrichment)")
    call.add_argument("--min-reads", type=int, default=3,
                      help="binomial: test only positions with at least this many reads "
                           "in one enriched library; local: in the pooled replicates "
                           "(default 3)")
    call.add_argument("--window", type=int, default=100,
                      help="local: background window either side of a site (default 100)")
    call.add_argument("--gap", type=int, default=2,
                      help="local: nt either side of a site left out of its background "
                           "(default 2)")
    call.add_argument("--trim", type=int, default=5,
                      help="local: busiest window positions left out of the background, "
                           "so a neighbouring TSS does not count (default 5)")
    call.add_argument("--local-dispersion", type=float, default=DEFAULT_LOCAL_DISPERSION,
                      help="local: negative-binomial dispersion of counts around the "
                           "background (default %g)" % DEFAULT_LOCAL_DISPERSION)
    call.add_argument("--dispersion", default="auto",
                      help="binomial: variation between libraries beyond sampling noise "
                           "(beta-binomial test); 'auto' estimates it from enriched "
                           "replicates, or uses %g with one replicate; 0 gives the plain "
                           "binomial test (default auto)" % DEFAULT_DISPERSION)
    call.add_argument("--normalise", choices=("total", "exclude-rna"), default="exclude-rna",
                      help="binomial, local: library sizes from all reads, or leaving out reads "
                           "in rRNA and tRNA genes of --annotation (default exclude-rna; "
                           "falls back to total without an annotation)")
    call.add_argument("--preset", choices=sorted(PRESETS), default="default",
                      help="Starting set of thresholds; 'ettwiller2016' "
                           "reproduces Ettwiller et al. 2016. Options given "
                           "explicitly override the preset (default: default)")
    call.add_argument("--min-rpm", type=float,
                      help="Minimum 5'-end reads per million in the enriched "
                           "library (default 1.0)")
    call.add_argument("--min-log2fc", type=float,
                      help="Minimum log2(enriched/control) (default 1.0)")
    call.add_argument("--pseudocount", type=float,
                      help="RPM pseudocount for the enrichment ratio (default 1.0)")
    call.add_argument("--control-floor", type=float,
                      help="Minimum control read count used in the ratio "
                           "(default 0)")
    call.add_argument("--min-replicates", type=int,
                      help="Replicates in which a position must pass (cutoff, binomial) "
                           "or have reads (local) (default: 1 for cutoff and local, all "
                           "for binomial)")
    call.add_argument("--cluster-distance", type=int,
                      help="Merge TSS on the same strand no more than this "
                           "many nt apart; 0 disables (default 5)")
    call.add_argument("--utr-length", type=int, default=300,
                      help="Maximum 5' UTR length for Primary/Secondary (default 300)")
    call.add_argument("--antisense-distance", type=int, default=100,
                      help="Antisense window either side of a gene (default 100)")
    call.add_argument("--out-prefix", required=True,
                      help="Writes <prefix>.gff, <prefix>.tsv and <prefix>.qc.tsv")
    call.add_argument("--strict-qc", action="store_true",
                      help="Stop if a library fails the 5'-end data check "
                           "(default: warn and continue)")

    qc = sub.add_parser("qc", help="Check that libraries look like 5'-end data")
    add_inputs(qc)
    qc.add_argument("--annotation", help="GFF3 annotation (adds the strand check)")
    qc.add_argument("--feature-types", default="CDS",
                    help="GFF feature types used as genes (default: CDS)")
    qc.add_argument("--out", help="Write the checks to this TSV")

    export = sub.add_parser(
        "annogesic-libs",
        help="Write wiggle files and --tex_notex_libs strings for ANNOgesic tss_ps")
    add_inputs(export)
    export.add_argument("--signal", choices=("coverage", "5prime"),
                        default="coverage",
                        help="Wiggle content: read coverage (what TSSpredator "
                             "was designed for) or 5'-end counts (default coverage)")
    export.add_argument("--condition", type=int, default=1,
                        help="ANNOgesic condition id (default 1)")
    export.add_argument("--positive-reverse", action="store_true",
                        help="Write reverse-strand values as positive numbers")
    export.add_argument("--outdir", required=True)

    compare = sub.add_parser(
        "compare", help="Compare TSS calls with a reference TSS set")
    compare.add_argument("--query", required=True,
                         help="TSS to evaluate (GFF/GTF, e.g. <prefix>.gff)")
    compare.add_argument("--reference", required=True,
                         help="Reference TSS (GFF/GTF, or a table with --ref-columns)")
    compare.add_argument("--ref-columns",
                         help="1-based seqid,position,strand columns of a "
                              "non-GFF reference table, e.g. 1,2,3")
    compare.add_argument("--rename", action="append", default=[],
                         help="OLD=NEW sequence name in the reference "
                              "(can be repeated)")
    compare.add_argument("--windows", default="0,5,25",
                         help="Comma-separated match windows in nt (default 0,5,25)")
    compare.add_argument("--out", help="Write per-query offsets to this TSV")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.command == "annogesic-libs":
        libs = export_annogesic(
            args.enriched, args.control, args.outdir, signal=args.signal,
            mate=args.mate, min_mapq=args.min_mapq, condition=args.condition,
            negative_reverse=not args.positive_reverse)
        print("Use with: annogesic tss_ps ... --tex_notex_libs \\")
        print("    " + " \\\n    ".join(libs))
        return 0

    if args.command == "compare":
        rename = dict(item.split("=", 1) for item in args.rename)
        columns = (tuple(int(c) for c in args.ref_columns.split(","))
                   if args.ref_columns else None)
        query = read_sites(args.query)
        reference = read_sites(args.reference, columns, rename)
        windows = [int(w) for w in args.windows.split(",")]
        summary, rows = compare_sites(query, reference, windows)
        print(format_comparison(summary))
        if args.out:
            with open(args.out, "w") as out:
                out.write("seqid\tstrand\tposition\toffset_to_nearest_reference\n")
                for seqid, strand, pos, offset in rows:
                    out.write("%s\t%s\t%d\t%s\n" % (
                        seqid, strand, pos, "NA" if offset is None else offset))
        return 0

    if args.command == "qc":
        genes = (read_genes(args.annotation, tuple(args.feature_types.split(",")))
                 if args.annotation else None)
        paths = args.enriched + args.control
        signals = [load_library(p, args.mate, args.min_mapq) for p in paths]
        roles = ["enriched"] * len(args.enriched) + ["control"] * len(args.control)
        flagged = run_qc(paths, signals, roles, genes, args.out, stream=sys.stdout)
        return 2 if flagged else 0

    settings = dict(PRESETS[args.preset])
    for name in settings:
        value = getattr(args, name, None)
        if value is not None:
            settings[name] = value

    enriched = [load_library(p, args.mate, args.min_mapq) for p in args.enriched]
    controls = [load_library(p, args.mate, args.min_mapq) for p in args.control]
    genes = (read_genes(args.annotation, tuple(args.feature_types.split(",")))
             if args.annotation else None)
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    flagged = run_qc(args.enriched + args.control, enriched + controls,
                     ["enriched"] * len(enriched) + ["control"] * len(controls),
                     genes, args.out_prefix + ".qc.tsv")
    if flagged and args.strict_qc:
        sys.stderr.write("Stopping: %d librar%s failed the check (--strict-qc).\n"
                         % (flagged, "y" if flagged == 1 else "ies"))
        return 2
    min_replicates = args.min_replicates
    exclude = None
    if args.method in ("binomial", "local") and args.normalise == "exclude-rna":
        if args.annotation:
            exclude = GeneIndex(read_genes(args.annotation, ("rRNA", "tRNA")))
        else:
            sys.stderr.write("No --annotation: library sizes include rRNA/tRNA reads\n")
    if args.method == "local":
        if len(enriched) < 2:
            sys.stderr.write("Warning: --method local with one library cannot require reads "
                             "in several replicates; isolated read stacks (e.g. PCR "
                             "duplicates) are then called and the FDR is too optimistic\n")
        calls = call_local(enriched, controls or None, fdr=args.fdr, min_reads=args.min_reads,
                           min_replicates=min_replicates or 1,
                           window=args.window, gap=args.gap, trim=args.trim,
                           dispersion=args.local_dispersion, exclude=exclude)
    elif args.method == "binomial":
        if len(controls) != len(enriched):
            sys.stderr.write("--method binomial needs one --control per --enriched library\n")
            return 1
        if args.dispersion == "auto":
            dispersion = estimate_dispersion(enriched)
            source = "estimated from replicates"
            if dispersion is None:
                dispersion, source = DEFAULT_DISPERSION, "default (one replicate)"
        else:
            dispersion, source = float(args.dispersion), "given"
        sys.stderr.write("Dispersion %.4f (%s)\n" % (dispersion, source))
        calls = call_binomial(enriched, controls, fdr=args.fdr, min_fold=args.min_fold,
                              min_reads=args.min_reads, dispersion=dispersion,
                              min_replicates=min_replicates or len(enriched),
                              exclude=exclude)
    else:
        calls = call_positions(enriched, controls, min_rpm=settings["min_rpm"],
                               min_log2fc=settings["min_log2fc"],
                               min_replicates=min_replicates or 1,
                               pseudocount=settings["pseudocount"],
                               control_floor=settings["control_floor"],
                               rpm_decimals=settings["rpm_decimals"])
    calls = cluster_tss(calls, settings["cluster_distance"])
    mark_shoulders(calls, args.shoulder_window)
    if args.shoulder_ratio > 0:
        before = len(calls)
        calls = [t for t in calls if t.peak_ratio >= args.shoulder_ratio]
        sys.stderr.write("Shoulders removed: %d\n" % (before - len(calls)))
    if args.genome:
        add_promoter_scores(calls, read_fasta(args.genome), stream=sys.stderr)
    if genes is not None:
        classify(calls, genes, args.utr_length, args.antisense_distance)
    write_gff(calls, args.out_prefix + ".gff")
    write_table(calls, args.out_prefix + ".tsv")
    print(summarise(calls))
    return 0


if __name__ == "__main__":
    sys.exit(main())
