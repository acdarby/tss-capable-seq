import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "benchmarks"))

try:
    import numpy as np
    import harness
    import simulate
except ImportError:  # numpy is only needed for the benchmarks
    np = None


@unittest.skipIf(np is None, "numpy not installed")
class SimulateTests(unittest.TestCase):
    def test_reproducible_and_planted(self):
        a = simulate.simulate(length=50_000, n_tss=40, replicates=2, seed=1)
        b = simulate.simulate(length=50_000, n_tss=40, replicates=2, seed=1)
        self.assertEqual(a["truth"], b["truth"])
        self.assertEqual(a["enriched"], b["enriched"])
        self.assertEqual(len(a["truth"]), 40)
        self.assertEqual(len(a["enriched"]), 2)
        self.assertEqual(a["controls"], [])

    def test_strong_tss_have_reads_at_the_site(self):
        sim = simulate.simulate(length=50_000, n_tss=20, tss_mu=6, tss_sigma=0.1,
                                background=0.0, dispersion=0.0, seed=2)
        sig = sim["enriched"][0]
        for seqid, strand, pos in sim["truth"]:
            self.assertGreater(sig[(seqid, strand)].get(pos, 0), 50)

    def test_control_has_processed_ends_and_rrna(self):
        sim = simulate.simulate(length=50_000, n_tss=20, n_processed=10, with_control=True,
                                rrna_fraction=0.5, dispersion=0.0, seed=3)
        ctl = sim["controls"][0]
        total = sum(sum(d.values()) for d in ctl.values())
        block = sum(v for p, v in ctl[("sim", "+")].items() if 25_000 < p <= 30_000)
        self.assertGreater(block / total, 0.4)
        for seqid, strand, pos in sim["processed"]:
            self.assertGreater(ctl[(seqid, strand)].get(pos, 0), 100)

    def test_profile_from_signal(self):
        sig = {("c", "+"): {100: 50.0, 101: 10.0, 500: 1.0}, ("c", "-"): {300: 20.0}}
        prof = simulate.profile_from_signal(sig, [("c", "+", 100), ("c", "-", 300)], 1000)
        self.assertEqual(prof["n_tss"], 2)
        self.assertAlmostEqual(prof["background"], 1.0 / 2000)
        self.assertAlmostEqual(prof["tss_mu"], (np.log(60) + np.log(20)) / 2)


@unittest.skipIf(np is None, "numpy not installed")
class MeasureTests(unittest.TestCase):
    def test_match_window_and_strand(self):
        ref = [("+", 100), ("-", 200)]
        self.assertEqual(harness.match([("+", 102), ("-", 100), ("+", 300)], ref, 2), 1)
        self.assertEqual(harness.match([("+", 103)], ref, 2), 0)

    def test_thin_keeps_expected_fraction(self):
        rng = np.random.default_rng(0)
        sig = {("c", "+"): {i: 100.0 for i in range(1, 1001)}}
        out = harness.thin(sig, 0.25, rng)
        total = sum(out[("c", "+")].values())
        self.assertAlmostEqual(total / 100000, 0.25, delta=0.01)

    def test_box_detection(self):
        # TSS at 101 on the + strand with TATAAT at -12..-7.
        seq = "G" * 88 + "TATAAT" + "G" * 400
        self.assertTrue(harness.has_box(harness.upstream(seq, "+", 101, 15, 5), "TATAAT"))
        self.assertFalse(harness.has_box(harness.upstream(seq, "+", 301, 15, 5), "TATAAT"))

    def test_coverage_rise(self):
        cov = {"+": np.zeros(400, dtype=np.float32), "-": np.zeros(400, dtype=np.float32)}
        cov["+"][200:] = 10          # rise at position 200 on +
        cov["-"][:201] = 10          # rise at position 200 on - (transcribed leftwards)
        self.assertEqual(harness.coverage_rise([("+", 200)], cov), 1.0)
        self.assertEqual(harness.coverage_rise([("-", 200)], cov), 1.0)
        self.assertEqual(harness.coverage_rise([("+", 300)], cov), 0.0)


if __name__ == "__main__":
    unittest.main()
