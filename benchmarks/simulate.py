#!/usr/bin/env python3
"""Synthetic 5'-end Cappable-seq libraries with planted TSS.

The model is deliberately simple and every part of it is a parameter:

* background: read 5' ends from degradation and fragmentation, Poisson with a
  rate that changes along the genome in segments (expression varies by region);
* TSS: planted at random positions with strengths drawn from a log-normal
  distribution; reads start at the TSS or within a few nucleotides of it
  (imprecise initiation);
* processed ends (control only, or depleted in the enriched library): strong
  in the non-enriched control, weak after enrichment;
* rRNA-like pile-up (control only): a block that takes a large share of control
  reads, as in the E. coli control (76% rRNA);
* replicates: each site's rate is multiplied per replicate by a gamma variable
  with mean 1 and variance `dispersion` (biological variability), then counts
  are Poisson.

`simulate()` returns libraries in the signal format used by cappable_tss.py
({(seqid, strand): {position: count}}) and the set of true TSS.
"""
import math

import numpy as np

SEQID = "sim"


def profile_from_signal(signal, tss_sites, length, window=3):
    """Simulation parameters that mimic one real library.

    `signal` is a 5'-end signal ({(seqid, strand): {pos: count}}) of one
    replicate, `tss_sites` the (seqid, strand, pos) of TSS called on it and
    `length` the genome length. Strengths are the reads within `window` nt of
    each TSS; they are truncated by the calling threshold, so the fitted
    log-normal under-represents weak TSS (the simulation is then conservative
    about how many weak TSS exist).
    """
    near = set()
    strengths = []
    for seqid, strand, pos in tss_sites:
        d = signal.get((seqid, strand), {})
        total = 0.0
        for q in range(pos - window, pos + window + 1):
            if q in d:
                total += d[q]
                near.add((seqid, strand, q))
        if total > 0:
            strengths.append(total)
    all_reads = sum(sum(d.values()) for d in signal.values())
    tss_reads = sum(strengths)
    logs = [math.log(x) for x in strengths]
    mu = sum(logs) / len(logs)
    sigma = (sum((x - mu) ** 2 for x in logs) / max(1, len(logs) - 1)) ** 0.5
    return {"length": int(length), "n_tss": len(strengths), "tss_mu": mu,
            "tss_sigma": sigma, "background": (all_reads - tss_reads) / (2.0 * length)}


def _segments(rng, length, mean_len, mean_rate, sigma):
    """Piecewise-constant background rate along one strand."""
    rate = np.empty(length)
    pos = 0
    while pos < length:
        seg = max(50, int(rng.exponential(mean_len)))
        rate[pos:pos + seg] = mean_rate * rng.lognormal(-sigma ** 2 / 2, sigma)
        pos += seg
    return rate


def _to_signal(counts_by_strand):
    signal = {}
    for strand, counts in counts_by_strand.items():
        idx = np.nonzero(counts)[0]
        signal[(SEQID, strand)] = {int(i) + 1: float(counts[i]) for i in idx}
    return signal


def simulate(length=1_000_000, n_tss=1000, replicates=3, depth=1.0,
             tss_mu=3.0, tss_sigma=1.5, background=0.05, background_sigma=1.0,
             segment_length=1000, spread=(0.6, 0.15, 0.05), dispersion=0.1,
             with_control=False, control_tss_factor=0.1, n_processed=300,
             processed_strength=200.0, control_background=1.0,
             rrna_fraction=0.0, seed=0):
    """Simulate enriched (and optionally control) libraries.

    Expected counts are per library at `depth` 1.0; `depth` scales everything.
    `spread` gives the share of a TSS's reads at offsets 0, ±1, ±2 (the rest
    of the mass is split equally); `tss_mu`/`tss_sigma` are the log-normal
    parameters of TSS strength (median exp(tss_mu) reads at depth 1.0).
    `control_tss_factor` is the TSS rate in the control relative to the
    enriched library (enrichment = 1 / factor). `rrna_fraction` is the share
    of control reads placed in an rRNA-like block.

    Returns a dict with "enriched" and "controls" (lists of signals), "truth"
    (set of (seqid, strand, position)), "strength" (dict of expected reads per
    true TSS at depth 1.0) and "processed" (set of processed-end positions).
    """
    rng = np.random.default_rng(seed)
    strands = ("+", "-")

    # True TSS, at least 30 nt apart on the same strand, away from the ends.
    truth, strength = set(), {}
    for strand in strands:
        n = n_tss // 2 + (n_tss % 2 if strand == "+" else 0)
        positions = set()
        while len(positions) < n:
            p = int(rng.integers(200, length - 200))
            if all(abs(p - q) >= 30 for q in range(p - 30, p + 31) if q in positions):
                positions.add(p)
        for p in positions:
            truth.add((SEQID, strand, p))
            strength[(SEQID, strand, p)] = float(rng.lognormal(tss_mu, tss_sigma))

    processed = set()
    for strand in strands:
        for p in rng.choice(np.arange(200, length - 200), n_processed // 2, replace=False):
            if (SEQID, strand, int(p)) not in truth:
                processed.add((SEQID, strand, int(p)))

    rest = (1.0 - spread[0] - 2 * spread[1] - 2 * spread[2]) / 2.0
    shape = {0: spread[0], -1: spread[1], 1: spread[1], -2: spread[2], 2: spread[2],
             -3: rest / 2, 3: rest / 2}

    bg_rate = {s: _segments(rng, length, segment_length, background, background_sigma)
               for s in strands}

    def library(tss_factor, bg_factor, processed_rate, rrna):
        rate = {s: bg_rate[s] * bg_factor * depth for s in strands}
        for (_, s, p), st in strength.items():
            for off, frac in shape.items():
                q = p + off if s == "+" else p - off
                rate[s][q - 1] += st * tss_factor * frac * depth
        for (_, s, p) in processed:
            rate[s][p - 1] += processed_rate * depth
        if rrna > 0:
            block = slice(length // 2, length // 2 + 5000)
            total = sum(r.sum() for r in rate.values())
            rate["+"][block] += rrna / (1 - rrna) * total / 5000
        return rate

    def replicate(rate):
        counts = {}
        for s, r in rate.items():
            mult = rng.gamma(1 / dispersion, dispersion, size=r.size) if dispersion > 0 else 1.0
            counts[s] = rng.poisson(r * mult)
        return _to_signal(counts)

    enriched_rate = library(1.0, 1.0, processed_strength * 0.05, 0.0)
    out = {"enriched": [replicate(enriched_rate) for _ in range(replicates)],
           "controls": [], "truth": truth, "strength": strength, "processed": processed,
           "length": length, "rrna_block": (length // 2 + 1, length // 2 + 5000)}
    if with_control:
        control_rate = library(control_tss_factor, control_background, processed_strength,
                               rrna_fraction)
        out["controls"] = [replicate(control_rate) for _ in range(replicates)]
    return out
