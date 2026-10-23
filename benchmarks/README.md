# Benchmarks

Tools to judge TSS-calling methods on real and simulated Cappable-seq data
(plan: `../docs/benchmark_plan.md`).

| File | Purpose |
|---|---|
| `datasets.toml` | Dataset manifest: files (relative to a data directory), accessions, replicate structure, mapping quality, which mate carries the RNA 5' end, optional reference TSS and orthogonal data |
| `harness.py` | Runs every method on every dataset and writes one table of measures |
| `simulate.py` | Synthetic 5'-end libraries with planted TSS, optionally calibrated from a real library |
| `calibration.py` | Null tests for the statistical calls (false calls where there is nothing to find) |
| `matched.py` | Compares rankings (RPM, local test, control test) at the same number of TSS on independent evidence |
| `sensitivity.py` | Settings of the local-background test against call numbers and promoter rates |
| `ettwiller2016/` | Rebuilds the *E. coli* data and compares with the published TSS table |
| `frontiers2026/` | TSSpredator and ANNOgesic on wAlbB fragment-library data |
| `tools/` | Download SRA runs, list runs per project, annotate genomes, TSS statistics |
| `results/` | Results tables (TSV) |

```bash
# real data (BAMs and genomes laid out as in datasets.toml under DATA)
python benchmarks/harness.py real --data-dir DATA --out results/real
# simulation calibrated from the wAlbB triplicate, three sequencing depths
python benchmarks/harness.py sim --profile walbb --data-dir DATA --depths 1,0.25,0.05 --out results/sim
```

Methods are registered in `harness.METHODS`; a new approach is a function
`(enriched_signals, control_signals) -> list of TSS`.

## Measures

| Measure | Meaning |
|---|---|
| `n_tss`, `share_<class>` | TSS on the chromosome and their class shares |
| `box10_*`, `box35_*` | Share with a TATAAT-like −10 box (−15..−5) or TTGACA-like −35 box (−40..−28), one mismatch allowed; `box10_background` uses the same positions shifted 100 nt |
| `ref_precision_w`, `ref_recall_w` | Agreement with the reference TSS within `w` nt (*E. coli* Table S1) |
| `loo_support` | Leave-one-replicate-out: share of TSS called on the other replicates that the held-out replicate also calls (±2 nt) |
| `downF_n`, `downF_recovery` | TSS called after keeping a fraction F of reads, and the share of full-depth TSS recovered (±2 nt) |
| `ortho_rise` | Share of TSS with ≥ 1.5× coverage rise (50 nt downstream vs upstream) in orthogonal fragment-library data; `_background` for shifted positions |
| `precision`, `recall`, `recall_<tier>` (simulation) | Against planted TSS (±3 nt); tiers are thirds of TSS strength |
| `processed_called` (simulation) | Share of processed ends wrongly called as TSS |

## Baselines (October 2026)

`ettwiller2016` = published pipeline (≥ 1.5 RPM, enrichment ≥ 0 when a
control exists, TSS in all replicates); `suliman2025` = ≥ 10 CPM, 10 nt
clusters, all replicates, no control. Full tables in `results/`.

| | *E. coli* E / S | *Salmonella* E / S | wAlbB E / S | wMelPop-CLA E / S |
|---|---|---|---|---|
| TSS | 16,639 / 4,079 | 10,609 / 3,483 | 4,724 / 1,527 | 2,119 / 742 |
| −10 box, all (background) | 25% / 32% (2.5%) | 33% / 37% (2.5%) | 49% / 54% (8%) | 57% / 63% (7%) |
| −10 box, internal TSS | 19% / 25% | 30% / 34% | 45% / 50% | 55% / 62% |
| Table S1 recall, exact | 96% / 23% | – | – | – |
| Held-out replicate support | – | 88% / 94% | 94% / 97% | 89% / 95% |
| Recovered at 25% of reads | 89% / 96% | 88% / 95% | 88% / 95% | 76% / 84% |
| TSS called at 25% of reads | 19,185 / 4,109 | 10,025 / 3,427 | 4,384 / 1,470 | 1,833 / 649 |
| Orthogonal coverage rise (background) | – | – | 56% / 73% (9% / 6%) | – |

Simulation calibrated from wAlbB (5,189 planted TSS, three replicates):

| Depth | `ettwiller2016` precision / recall | `suliman2025` precision / recall |
|---|---|---|
| 1 | 0.94 / 0.95 | 1.00 / 0.68 |
| 0.25 | 0.76 / 0.94 | 1.00 / 0.64 |
| 0.05 | 0.99 / 0.71 | 1.00 / 0.55 |

Both fixed cut-offs depend on depth: at lower depth one read passes 1.5 RPM,
so the Ettwiller cut-off calls more noise (on real *E. coli* data, 25% of
the reads gives more TSS than 100%), while 10 CPM keeps precision but misses
most weak TSS (11% of the weakest third recovered at full depth). With a
control, the Ettwiller enrichment filter removes almost all processed ends
(1% called, against 85% without a control).

## Statistical calling against the control (October 2026)

`call --method binomial`: a beta-binomial test of enriched against control
counts, Fisher's method across replicates, Benjamini–Hochberg FDR 5%; library
sizes leave out rRNA and tRNA reads. Results tables: `results/calibration_*.tsv`,
`results/statistical_ecoli.tsv`, `results/matched_ecoli.tsv`,
`results/sim_walbb_disp*.tsv`.

**Dispersion.** The variance of library-to-library variation, estimated from
replicates (`cappable_tss.estimate_dispersion`): wAlbB 0.012, *Salmonella*
0.018, wMelPop-CLA 0.035, *E. coli* enriched against its technical replicate
0.046. With one replicate, `call` uses 0.05.

**Calibration.** Calls at FDR 5% where every call is false
(`calibration.py`):

| Null test | binomial (dispersion 0) | beta-binomial 0.02 | beta-binomial 0.05 |
|---|---|---|---|
| Simulated, no extra variation | 0 | 0 | 0 |
| Simulated, dispersion 0.05 | 87 | 1 | 0 |
| Simulated, dispersion 0.1 | 169 | 28 | 0 |
| One library split into random halves (all datasets) | 0 | 0 | 0 |
| *Salmonella* replicate 1 against replicate 2 | 624 | 1 | 0 |
| wAlbB replicate 1 against replicate 2 | 36 | 1 | 1 |
| wMelPop-CLA replicate 1 against replicate 2 | 53 | 11 | 3 |

The plain binomial test is valid only for sampling noise (split halves);
real replicates vary more, and the beta-binomial with a dispersion at or
above the estimate removes nearly all false calls. Swapping the *E. coli*
labels (control as "enriched") gives 18,000 clustered calls with a −10 box
rate of 3.6%, near the 2.5% background: positions really more abundant in
the control (processed ends), not promoters.

***E. coli*, real data** (`harness.py real`; one enriched library, so
dispersion 0.05):

| | `ettwiller2016` | binomial | beta-binomial |
|---|---|---|---|
| TSS | 16,639 | 11,412 | 7,259 |
| −10 box (background 2.7%) | 24.7% | 26.6% | 28.0% |
| Table S1 recall, exact | 96% | 64% | 39% |
| TSS called at 25% of reads | 19,185 | 4,170 | 2,506 |

**Matched comparison** (first version of `matched.py`, `results/matched_ecoli.tsv`): every position with ≥ 3 reads ranked by
each score and the top N compared (−10 box rate):

| Top N | RPM (Ettwiller filter) | fold enrichment | p-value, dispersion 0 | 0.02 | 0.05 |
|---|---|---|---|---|---|
| 1,000 | **35.9%** | 23.9% | 35.0% | 31.5% | 29.5% |
| 4,000 | **31.1%** | 23.9% | 30.4% | 28.9% | 27.6% |
| 7,000 | 27.7% | 23.7% | 28.2% | 27.7% | 27.0% |
| 11,000 | 25.3% | 22.8% | 25.6% | 25.4% | 25.3% |
| 16,000 | 23.5% | 21.8% | 23.4% | 23.5% | 23.2% |

The p-value ranking is no better than ranking by read count, and the higher
promoter rate of the statistical calls comes only from calling fewer TSS.
Fold enrichment over the control is the weakest ranking: the
non-enriched control also contains TSS 5' ends (median enrichment of Table S1
TSS with ≥ 20 reads is 4.6-fold; 13% are less than 2-fold enriched). The
TSS that the Ettwiller cut-off adds beyond rank 11,000 still have a −10 box
rate of 20%, against 2.5% background, so most are real: at FDR 5% the test
against the control stops too early.

**Power** (simulation calibrated from wAlbB with a control, recall at
precision ≥ 0.99 for all methods):

| Simulation | `ettwiller2016` | binomial | beta-binomial (dispersion estimated) |
|---|---|---|---|
| Dispersion 0.015 (as wAlbB), depth 1 | 0.96 | 0.78 | 0.76 |
| depth 0.25 | 0.89 | 0.47 | 0.45 |
| depth 0.05 | 0.59 | 0.12 | 0.11 |
| Dispersion 0.035 (as wMelPop), depth 1 | 0.95 | 0.78 | 0.73 |
| Dispersion 0.1, depth 1 | 0.94 | 0.75 | 0.56 |

At realistic dispersion the beta-binomial costs little compared with the
binomial, but both lose most weak TSS, especially at low depth, because a
few reads cannot show enrichment over a control that also holds TSS reads.
A null of at most 2-fold enrichment (`--min-fold 2`) loses almost
everything.

**Conclusion.** The test against the control is calibrated and robust to
depth (it calls fewer TSS with fewer reads, where the RPM cut-off calls
more), but enrichment over the control is a weak signal. It should be a
filter (not depleted, as in the Ettwiller pipeline), not the test that
decides. The evidence for a TSS is a 5'-end pile-up above the local
background, which is approach D, needed anyway for datasets without a
control. `call` keeps `--method cutoff` as the default.

## Test against the local background (approach D, October 2026)

`call --method local`: replicates pooled; a position with ≥ 3 pooled reads and
reads in every replicate is tested against the 5'-end rate of its
neighbourhood (±100 nt, leaving out ±2 nt and the 5 busiest positions, one
pseudo-read) with a negative binomial; BH FDR 5% over every position of the
genome; with a control, sites depleted after enrichment are removed (the
Ettwiller rule). Results: `results/matched_local_*`, `results/local_real_*`,
`results/sim_walbb_local_bias*.tsv`.

**What the background looks like.** Read 5' ends are spread very unevenly:
where the background is 2–8 reads per position, 44% of positions (E. coli,
after trimming) have no reads, against 2% for Poisson; the excess at a
position correlates 0.94 between wAlbB replicates, so the unevenness is
fixed (ligation preferences, processed ends), not noise. A Poisson test
would be far too lenient, and no single negative binomial fits (maximum
likelihood gives dispersions of 40–1000, which leave the test powerless).
The dispersion is therefore a fixed setting (2); the ranking hardly changes
between 0.5 and 10. Replicates are pooled rather than combined with
Fisher's method because the unevenness repeats in each replicate. Testing
only positions with ≥ 3 reads with unconditional p-values requires counting
all positions in the FDR; counting only the tested ones called almost every
tested position in a null test.

**Matched comparison** (`matched.py`, same candidates for each ranking; −10
box rate of the top N):

| Top N | E. coli RPM / local | wAlbB RPM / local | Salmonella RPM / local | wMelPop-CLA RPM / local |
|---|---|---|---|---|
| 1,000 | 35.9% / 35.0% | 55.0% / 55.9% | 38.7% / 38.6% | 60.4% / 59.7% |
| 4,000 | 31.1% / **35.8%** | 48.7% / 50.0% | 36.2% / 36.8% | 52.0% / 52.7% |
| 8,000 | 27.0% / **34.1%** | 45.0% / 47.5% | 34.2% / 35.0% | – |
| 16,000 | 23.5% / **30.5%** | 39.4% / 41.8% | 31.5% / 32.2% | – |

Ranking by the local test is better than ranking by read count in E. coli
(+1 purine also 86% against 79%), slightly better in wAlbB and *Salmonella*
and equal in wMelPop-CLA (where few candidates remain after requiring reads
in all three replicates). It also picks fewer sites at rRNA and tRNA genes
(0.1–1% against 0.3–2%). Count divided by background ranks worse than the
p-value.

**Where the FDR stops.** In wAlbB the local test stops at about 17,000 sites,
where the −10 box rate of further ranks falls from 26–30% to 17% and 11%
(background 8%); with replicates the stop is sensible. With one library
(E. coli) it stops at 52,000 and the box rate declines smoothly, because
isolated read stacks (e.g. PCR duplicates) cannot be removed by requiring
reads in several replicates; `call` warns in that case.

**Real data, called sets** (`harness.py real`):

| | E. coli E / local | Salmonella E / local | wAlbB E / local | wMelPop-CLA E / local |
|---|---|---|---|---|
| TSS | 16,639 / 51,064 | 10,609 / 37,071 | 4,724 / 16,707 | 2,119 / 3,924 |
| −10 box (background) | 25% / 18% (2.5%) | 33% / 26% (2.5%) | 49% / 41% (8%) | 57% / 53% (7%) |
| Held-out replicate support | – | 88% / 79% | 94% / 71% | 89% / 64% |
| TSS at 25% of reads | 19,185 / 18,397 | 10,025 / 18,755 | 4,384 / 9,295 | 1,833 / 1,809 |
| Orthogonal coverage rise (background) | – | – | 56% / 36% (9% / 12%) | – |

The local test calls 2–3.5 times more TSS. The additional sites are weak:
lower promoter rates (but 3–10 times background), less often found by the
held-out replicate alone and by a coverage rise in the fragment library,
and lost when reads are removed. Unlike the RPM cut-off, it calls fewer TSS
with fewer reads.

**Simulation** calibrated from wAlbB (dispersion 0.015, three replicates),
with a fixed per-position background bias of variance 2 as in real data
(`--position-bias 2`), precision / recall:

| | Depth | `ettwiller2016` | `suliman2025` | local |
|---|---|---|---|---|
| No control | 1 | 0.67 / 0.95 | 1.00 / 0.73 | 0.72 / 0.98 |
| | 0.25 | 0.47 / 0.94 | 1.00 / 0.67 | **0.87** / 0.93 |
| | 0.05 | 0.95 / 0.72 | 1.00 / 0.56 | 0.98 / 0.71 |
| With control | 1 | 0.88 / 0.95 | 1.00 / 0.73 | 0.93 / 0.98 |
| | 0.25 | 0.92 / 0.89 | 1.00 / 0.67 | 0.96 / 0.92 |
| | 0.05 | 0.99 / 0.64 | 1.00 / 0.56 | 1.00 / 0.69 |

Without a control the false calls of both are processed ends (planted at
5% of their control strength in the enriched library), which no method
without a control can tell from TSS.

**Parameter sensitivity** (`sensitivity.py`; TSS called, and −10 box rate of
the top 8,000 for wAlbB and E. coli, top 2,000 for wMelPop-CLA):

| Setting | wAlbB | wMelPop-CLA | E. coli |
|---|---|---|---|
| defaults | 16,710 (47.5%) | 3,924 (57.0%) | 51,064 (34.1%) |
| window ±50 / ±200 | 17,278 / 16,044 (47.4% / 47.2%) | 3,925 / 3,921 | 54,665 / 48,642 (34.1% / 33.9%) |
| trim 0 / 10 | 12,024 / 17,646 (46.6% / 47.5%) | 3,411 / 3,966 | 33,841 / 56,199 (34.1% / 33.9%) |
| minimum reads 2 / 5 | 16,710 / 15,583 | 3,924 / 3,455 | 51,385 / 41,002 |
| FDR 1% / 10% | 15,948 / 17,038 | 3,890 / 3,937 | 46,191 / 53,458 |

The ranking does not depend on the settings. The number of calls depends
mostly on trimming, and hardly on the FDR level: almost every position with
reads in all replicates and a quiet neighbourhood passes, so in practice the
call set is "≥ 3 pooled reads, in every replicate, standing out from the
trimmed background" and the q-value adds little to the stop.

**Conclusion.** The local test ranks TSS at least as well as read count in
every dataset (much better in E. coli), works without a control, and with
replicates stops near where calls stop looking like promoters. It calls
many more weak TSS than the RPM cut-off; whether those are wanted is a
choice of depth against confidence. The q-values are nominal (the
dispersion is a setting, not an estimate). `--method cutoff` stays the
default until the choice of stop is settled (plan section 7).
