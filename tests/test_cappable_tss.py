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
        prefix = os.path.join(d, "out")
        ct.main(["call", "--enriched", e, "--control", c, "--annotation", a,
                 "--out-prefix", prefix])
        with open(prefix + ".gff") as handle:
            lines = [l for l in handle if not l.startswith("#")]
        self.assertEqual(len(lines), 1)
        self.assertIn("type=Primary", lines[0])
        self.assertIn("associated_gene=geneA", lines[0])
        self.assertIn("UTR_length=Primary_50", lines[0])
        self.assertTrue(os.path.exists(prefix + ".tsv"))


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
