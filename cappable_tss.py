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
    ``tss_ptss`` subcommand, with the enriched library standing in for TEX+
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
        self.classes = set()
        self.genes = defaultdict(list)     # class -> [(gene name, UTR length)]

    @property
    def key(self):
        return (self.seqid, self.strand, self.pos)


def log2_enrichment(enriched_rpm, control_rpm, pseudocount):
    return math.log2((enriched_rpm + pseudocount) / (control_rpm + pseudocount))


def call_positions(enriched, controls, min_rpm=1.0, min_log2fc=1.0,
                   min_replicates=1, pseudocount=1.0):
    """Return candidate TSS positions before clustering.

    ``enriched`` is a list of signals, one per replicate; ``controls`` is a
    list of the same length (paired by replicate) or empty if the experiment
    has no non-enriched control, in which case only ``min_rpm`` applies.
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
            if controls:
                c_rpm = (controls[i].get((seqid, strand), {}).get(pos, 0.0)
                         * control_scales[i])
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
    """Merge calls on the same strand closer than ``distance`` nt.

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
            out.write("\t".join([tss.seqid, source, "TSS", str(tss.pos),
                                 str(tss.pos), _fmt(tss.enriched_rpm),
                                 tss.strand, ".", ";".join(attrs)]) + "\n")


def write_table(tss_list, path):
    header = ["seqid", "position", "strand", "type", "associated_gene",
              "UTR_length", "enriched_rpm", "control_rpm", "log2_enrichment",
              "replicates"]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for tss in tss_list:
            out.write("\t".join([
                tss.seqid, str(tss.pos), tss.strand, _class_text(tss),
                _gene_text(tss), _utr_text(tss), _fmt(tss.enriched_rpm),
                _fmt(tss.control_rpm), _fmt(tss.log2fc), str(tss.replicates),
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
    call.add_argument("--feature-types", default="gene",
                      help="Comma-separated GFF feature types used as genes "
                           "(default: gene; CDS is a common alternative)")
    call.add_argument("--min-rpm", type=float, default=1.0,
                      help="Minimum 5'-end reads per million in the enriched "
                           "library (default 1.0)")
    call.add_argument("--min-log2fc", type=float, default=1.0,
                      help="Minimum log2(enriched/control) (default 1.0)")
    call.add_argument("--pseudocount", type=float, default=1.0,
                      help="RPM pseudocount for the enrichment ratio (default 1.0)")
    call.add_argument("--min-replicates", type=int, default=1,
                      help="Replicates in which a position must pass (default 1)")
    call.add_argument("--cluster-distance", type=int, default=5,
                      help="Merge TSS on the same strand within this many nt; "
                           "0 disables (default 5)")
    call.add_argument("--utr-length", type=int, default=300,
                      help="Maximum 5' UTR length for Primary/Secondary (default 300)")
    call.add_argument("--antisense-distance", type=int, default=100,
                      help="Antisense window either side of a gene (default 100)")
    call.add_argument("--out-prefix", required=True,
                      help="Writes <prefix>.gff and <prefix>.tsv")

    export = sub.add_parser(
        "annogesic-libs",
        help="Write wiggle files and --tex_notex_libs strings for ANNOgesic tss_ptss")
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
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.command == "annogesic-libs":
        libs = export_annogesic(
            args.enriched, args.control, args.outdir, signal=args.signal,
            mate=args.mate, min_mapq=args.min_mapq, condition=args.condition,
            negative_reverse=not args.positive_reverse)
        print("Use with: annogesic tss_ptss ... --tex_notex_libs \\")
        print("    " + " \\\n    ".join(libs))
        return 0

    enriched = [load_library(p, args.mate, args.min_mapq) for p in args.enriched]
    controls = [load_library(p, args.mate, args.min_mapq) for p in args.control]
    calls = call_positions(enriched, controls, min_rpm=args.min_rpm,
                           min_log2fc=args.min_log2fc,
                           min_replicates=args.min_replicates,
                           pseudocount=args.pseudocount)
    calls = cluster_tss(calls, args.cluster_distance)
    if args.annotation:
        genes = read_genes(args.annotation, tuple(args.feature_types.split(",")))
        classify(calls, genes, args.utr_length, args.antisense_distance)
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    write_gff(calls, args.out_prefix + ".gff")
    write_table(calls, args.out_prefix + ".tsv")
    print(summarise(calls))
    return 0


if __name__ == "__main__":
    sys.exit(main())
