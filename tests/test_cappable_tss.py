import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import cappable_tss as ct  # noqa: E402


def write(directory, name, text):
    path = os.path.join(directory, name)
    with open(path, "w") as handle:
        handle.write(text)
    return path


def bed_reads(seqid, pos, strand, n, length=30):
    """n BED lines whose 5' end is at 1-based position pos."""
    if strand == "+":
        start, end = pos - 1, pos - 1 + length
    else:
        start, end = pos - length, pos
    return "".join("%s\t%d\t%d\tr\t0\t%s\n" % (seqid, start, end, strand)
                   for _ in range(n))


GFF = """##gff-version 3
chr\t.\tgene\t200\t400\t.\t+\t.\tID=g1;locus_tag=geneA
chr\t.\tgene\t600\t900\t.\t-\t.\tID=g2;locus_tag=geneB
"""


class FivePrimeTests(unittest.TestCase):
    def test_bed_five_prime_ends(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = write(tmp, "a.bed", bed_reads("chr", 150, "+", 3)
                        + bed_reads("chr", 950, "-", 2))
            signal = ct.five_prime_from_bed(bed)
        self.assertEqual(signal[("chr", "+")][150], 3)
        self.assertEqual(signal[("chr", "-")][950], 2)

    def test_wig_reads_negative_reverse_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            fwd = write(tmp, "f.wig", "track x\nvariableStep chrom=chr\n10 5\n")
            rev = write(tmp, "r.wig", "fixedStep chrom=chr start=20 step=1\n-4\n0\n-1\n")
            signal = ct.load_library(fwd + "," + rev)
        self.assertEqual(signal[("chr", "+")][10], 5)
        self.assertEqual(signal[("chr", "-")][20], 4)
        self.assertEqual(signal[("chr", "-")][22], 1)
        self.assertNotIn(21, signal[("chr", "-")])


class CallingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = self.tmp.name
        # Enriched: TSS at 150(+, primary for geneA), 152(+, weaker
        # neighbour, clustered away), 120(+, secondary), 950(-, primary for
        # geneB), 300(+, internal), 450(-, antisense: within 100 nt of geneA),
        # 2000(+, orphan) and 700(+, equally strong in the control).
        enriched = (bed_reads("chr", 150, "+", 50) + bed_reads("chr", 152, "+", 10)
                    + bed_reads("chr", 120, "+", 20)
                    + bed_reads("chr", 950, "-", 40) + bed_reads("chr", 300, "+", 30)
                    + bed_reads("chr", 450, "-", 30) + bed_reads("chr", 2000, "+", 30)
                    + bed_reads("chr", 700, "+", 30))
        # Control: position 700 is as strong in the control -> processed end.
        control = (bed_reads("chr", 150, "+", 2) + bed_reads("chr", 700, "+", 30)
                   + bed_reads("chr", 5000, "+", 140))
        self.enriched = ct.five_prime_from_bed(write(d, "e.bed", enriched))
        self.control = ct.five_prime_from_bed(write(d, "c.bed", control))
        self.genes = ct.read_genes(write(d, "a.gff", GFF))

    def tearDown(self):
        self.tmp.cleanup()

    def called(self, **kwargs):
        calls = ct.call_positions([self.enriched], [self.control], **kwargs)
        calls = ct.cluster_tss(calls, 5)
        ct.classify(calls, self.genes)
        return {(t.pos, t.strand): t for t in calls}

    def test_control_removes_non_enriched_ends(self):
        calls = self.called()
        self.assertNotIn((700, "+"), calls)
        self.assertIn((150, "+"), calls)

    def test_clustering_keeps_strongest(self):
        calls = self.called()
        self.assertNotIn((152, "+"), calls)

    def test_classes(self):
        calls = self.called()
        self.assertEqual(calls[(150, "+")].classes, {"Primary"})
        self.assertEqual(calls[(120, "+")].classes, {"Secondary"})
        self.assertEqual(calls[(300, "+")].classes, {"Internal"})
        self.assertEqual(calls[(950, "-")].classes, {"Primary"})
        self.assertEqual(calls[(450, "-")].classes, {"Antisense"})
        self.assertEqual(calls[(2000, "+")].classes, {"Orphan"})
        self.assertEqual(calls[(150, "+")].genes["Primary"], [("geneA", 50)])
        self.assertEqual(calls[(950, "-")].genes["Primary"], [("geneB", 50)])

    def test_replicate_requirement(self):
        calls = ct.call_positions([self.enriched, self.control],
                                  [self.control, self.control], min_replicates=2)
        self.assertEqual(calls, [])

    def test_no_control(self):
        calls = ct.call_positions([self.enriched], [])
        self.assertIn(700, [t.pos for t in calls])

    def test_cli_outputs(self):
        d = self.tmp.name
        e = write(d, "e2.bed", bed_reads("chr", 150, "+", 50))
        c = write(d, "c2.bed", bed_reads("chr", 999, "+", 50))
        a = write(d, "a2.gff", GFF)
        prefix = os.path.join(d, "results", "out")   # directory not yet made
        ct.main(["call", "--enriched", e, "--control", c, "--annotation", a,
                 "--out-prefix", prefix])
        with open(prefix + ".gff") as handle:
            lines = [l for l in handle if not l.startswith("#")]
        self.assertEqual(len(lines), 1)
        self.assertIn("type=Primary", lines[0])
        self.assertIn("associated_gene=geneA", lines[0])
        self.assertIn("UTR_length=Primary_50", lines[0])
        self.assertTrue(os.path.exists(prefix + ".tsv"))


class EttwillerPresetTests(unittest.TestCase):
    """Behaviour of filter_tss.pl / bam2firstbasegtf.pl (Ettwiller et al.)."""

    def signal(self, counts):
        sig = ct.new_signal()
        for pos, n in counts.items():
            sig[("chr", "+")][pos] = n
        return sig

    def call(self, enriched, control):
        preset = dict(ct.PRESETS["ettwiller2016"])
        preset.pop("cluster_distance")
        return {t.pos: t for t in ct.call_positions(
            [self.signal(enriched)], [self.signal(control)], **preset)}

    def test_missing_control_counts_as_one_read(self):
        # 1e6 reads in each library; 10 enriched reads = 10 RPM, control 0
        calls = self.call({1: 10, 2: 999990}, {5: 1000000})
        self.assertAlmostEqual(calls[1].control_rpm, 1.0)
        self.assertAlmostEqual(calls[1].log2fc, 3.321928, places=5)

    def test_equal_enrichment_passes_at_zero(self):
        calls = self.call({1: 10, 2: 999990}, {1: 10, 5: 999990})
        self.assertIn(1, calls)
        self.assertAlmostEqual(calls[1].log2fc, 0.0)

    def test_depleted_position_removed(self):
        calls = self.call({1: 10, 2: 999990}, {1: 11, 5: 999989})
        self.assertNotIn(1, calls)

    def test_rpm_rounded_before_cutoff(self):
        # 1499.6 reads in 1e9 -> 1.4996 RPM, rounds to 1.500 and passes
        calls = self.call({1: 1499.6, 2: 1e9 - 1499.6}, {5: 1e9})
        self.assertIn(1, calls)

    def test_cluster_distance_matches_bedtools_merge_d5(self):
        def tss(pos, rpm):
            return ct.TSS("chr", pos, "+", rpm, 0, 0, 1)
        merged = ct.cluster_tss([tss(100, 5), tss(106, 9), tss(113, 1)],
                                ct.PRESETS["ettwiller2016"]["cluster_distance"])
        self.assertEqual([t.pos for t in merged], [106, 113])


class CompareTests(unittest.TestCase):
    def test_windows_and_offsets(self):
        query = {("chr", "+", 100), ("chr", "+", 205), ("chr", "-", 500),
                 ("chr", "+", 900)}
        reference = {("chr", "+", 100), ("chr", "+", 200), ("chr", "-", 503),
                     ("chr", "-", 700)}
        summary, rows = ct.compare_sites(query, reference, windows=(0, 5))
        w0, w5 = summary["windows"]
        self.assertEqual((w0["query_matched"], w0["reference_recovered"]), (1, 1))
        self.assertEqual((w5["query_matched"], w5["reference_recovered"]), (3, 3))
        offsets = {(r[1], r[2]): r[3] for r in rows}
        self.assertEqual(offsets[("+", 205)], 5)    # 5 nt downstream
        self.assertEqual(offsets[("-", 500)], 3)    # downstream on - strand

    def test_read_table_with_columns_and_rename(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(tmp, "ref.csv", "chrom,pos,strand\nNC_1,10,forward\nNC_1,20,-\n")
            sites = ct.read_sites(path, columns=(1, 2, 3), rename={"NC_1": "chr"})
        self.assertEqual(sites, {("chr", "+", 10), ("chr", "-", 20)})


class ExportTests(unittest.TestCase):
    def test_annogesic_libs(self):
        with tempfile.TemporaryDirectory() as tmp:
            e = write(tmp, "enr.bed", "chr\t9\t12\tr\t0\t+\nchr\t20\t22\tr\t0\t-\n")
            c = write(tmp, "ctl.bed", "chr\t9\t10\tr\t0\t+\n")
            out = os.path.join(tmp, "wig")
            libs = ct.export_annogesic([e], [c], out)
            self.assertEqual(len(libs), 4)
            self.assertTrue(libs[0].endswith(":tex:1:a:+"))
            self.assertTrue(libs[3].endswith(":notex:1:a:-"))
            fwd = ct.read_wig(os.path.join(out, "enr_tex_forward.wig"), "+")
            with open(os.path.join(out, "enr_tex_reverse.wig")) as handle:
                rev_text = handle.read()
        self.assertEqual(sorted(fwd[("chr", "+")]), [10, 11, 12])
        self.assertIn("21 -1", rev_text)
        self.assertIn("22 -1", rev_text)


class LibraryQcTests(unittest.TestCase):
    def five_prime_like(self):
        # Sharp peaks: 100 reads at isolated positions, little background.
        sig = {("chr", "+"): {}}
        for p in range(100, 10000, 500):
            sig[("chr", "+")][p] = 100.0
            sig[("chr", "+")][p + 1] = 5.0
        for p in range(50, 10000, 997):
            sig[("chr", "+")][p] = 1.0
        return sig

    def fragment_like(self):
        # Read starts spread evenly along transcripts.
        return {("chr", "+"): {p: 3.0 for p in range(100, 10000)}}

    def test_peak_share_separates_library_types(self):
        self.assertGreater(ct.peak_share(self.five_prime_like()), 0.9)
        self.assertLess(ct.peak_share(self.fragment_like()), 0.05)

    def test_peak_share_ignores_single_reads(self):
        self.assertEqual(ct.peak_share({("chr", "+"): {100: 1.0, 5000: 1.0}}), 0.0)

    def test_sense_share(self):
        genes = [ct.Gene("chr", 1000, 2000, "+", "g1"), ct.Gene("chr", 3000, 4000, "-", "g2")]
        # 9 sense + 1 antisense read on each strand.
        sig = {("chr", "+"): {1500: 9.0, 3500: 1.0}, ("chr", "-"): {3500: 9.0, 1500: 1.0}}
        self.assertAlmostEqual(ct.sense_share(sig, genes), 0.9)

    def test_warnings(self):
        good = {"reads": 1e6, "positions": 1000, "peak_share": 0.8, "sense_share": 0.8}
        self.assertEqual(ct.qc_warnings(good, "enriched"), [])
        frag = dict(good, peak_share=0.25)
        self.assertEqual(len(ct.qc_warnings(frag, "enriched")), 1)
        self.assertEqual(ct.qc_warnings(dict(good, peak_share=0.4), "control"), [])
        anti = dict(good, sense_share=0.02)
        self.assertIn("antisense", ct.qc_warnings(anti, "enriched")[0])
        self.assertIn("no reads", ct.qc_warnings(dict(good, reads=0), "enriched")[0])

    def test_call_writes_qc_and_strict_stops(self):
        with tempfile.TemporaryDirectory() as d:
            frag = write(d, "frag.bed", "".join(
                "chr\t%d\t%d\tr\t0\t+\n" % (p, p + 50) for p in range(100, 3000) for _ in range(2)))
            prefix = os.path.join(d, "out")
            self.assertEqual(ct.main(["call", "--enriched", frag, "--out-prefix", prefix]), 0)
            with open(prefix + ".qc.tsv") as handle:
                self.assertIn("fragment", handle.read())
            self.assertEqual(ct.main(["call", "--enriched", frag, "--out-prefix", prefix,
                                      "--strict-qc"]), 2)
            self.assertEqual(ct.main(["qc", "--enriched", frag]), 2)


try:
    import pysam
except ImportError:
    pysam = None


@unittest.skipIf(pysam is None, "pysam not installed")
class BamTests(unittest.TestCase):
    def test_paired_end_uses_chosen_mate(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "x.bam")
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr", "LN": 1000}]}
            with pysam.AlignmentFile(path, "wb", header=header) as bam:
                for name, flag, start in (("p1", 0x1 | 0x40, 99),
                                          ("p1", 0x1 | 0x80 | 0x10, 199),
                                          ("s1", 0x10, 299),
                                          ("s2", 0x4, 0)):
                    read = pysam.AlignedSegment()
                    read.query_name = name
                    read.flag = flag
                    read.reference_id = 0 if not flag & 0x4 else -1
                    read.reference_start = start
                    read.mapping_quality = 30
                    read.cigartuples = [(0, 20)]
                    read.query_sequence = "A" * 20
                    bam.write(read)
            mate1 = ct.five_prime_from_bam(path, mate=1)
            mate2 = ct.five_prime_from_bam(path, mate=2)
        # read 1 forward at 0-based 99 -> 1-based 100
        self.assertEqual(dict(mate1[("chr", "+")]), {100: 1})
        # single-end reverse read at 299..318 -> 5' end 319
        self.assertEqual(dict(mate1[("chr", "-")]), {319: 1})
        # read 2 reverse at 199..218 -> 5' end 219, plus the single-end read
        self.assertEqual(dict(mate2[("chr", "-")]), {219: 1, 319: 1})


if __name__ == "__main__":
    unittest.main()
