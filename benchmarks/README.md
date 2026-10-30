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
| `promoter_model.py` | Genome-trained promoter score: recognition, rescue at low depth, learned weights |
| `shoulders.py` | Calls in the downstream shoulder of a stronger TSS: evidence for removed and kept calls |
| `standard_mode.py` | Validates the standard mode and its confidence tiers |
| `mixture/` | Two bacteria in one sample: *E. coli* mixed into *Salmonella* at 100, 10 and 2%, against each called alone |
| `design.py` | Design comparisons (replicates against one library, depth, control) for `docs/design_recommendations.md` |
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
| `prom_core10`, `prom_box10`, `prom_ext10`, `prom_box35`, `prom_purine` | Promoter score (from October 2026; see "Promoter score" below): σ70 −10 core, TATAAT, extended −10, −35 at the right spacing, purine at +1; `_bg` for composition-matched random positions |
| `box10_*`, `box35_*` | Earlier measure, kept for comparison: share with a TATAAT-like −10 box (−15..−5) or TTGACA-like −35 box (−40..−28), one mismatch allowed; `box10_background` uses the same positions shifted 100 nt |
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
reads in every replicate (as benchmarked here; the default is now reads in
one replicate, see approach F below) is tested against the 5'-end rate of its
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

## Replicates as evidence (approach F, October 2026)

The plan was to combine per-replicate p-values with Fisher's method. That is
not valid for the local test, because the uneven background repeats in every
replicate, so replicates are pooled (approach D). What replicates can still
show was tested on the three datasets with biological replicates:

* **Reads in every replicate are not needed.** Among weak calls, sites with
  reads in only some replicates are as promoter-like as sites with reads in
  all of them (−10 box rate, 3–5 pooled reads: wAlbB 32% in one replicate
  against 31% in all three; wMelPop-CLA 42% / 47%; *Salmonella*, 6–11 reads,
  26% / 23%). For weak sites, missing reads are sampling. `--method local`
  therefore needs reads in one replicate by default (`--min-replicates`).
  This adds 3,772 calls in wAlbB (−10 box 32%, background 8%), 1,855 in
  wMelPop-CLA (45%, background 7.5%) and 3,102 in *Salmonella* (20%,
  background 2.4%) (`results/replicates_real_*.tsv`). The exception is
  *E. coli* with its technical replicate, where sites found in one library
  only are weaker (12% against 17% at 6–11 reads), consistent with read
  stacks in a single library.
* **Rank reproducibility stops too early.** The share of sites that every
  replicate ranks highly on its own (the idea behind the irreproducible
  discovery rate, Li, Brown, Huang & Bickel 2011, *Ann. Appl. Stat.*
  5:1752–1779, doi:10.1214/11-AOAS466) falls from 97% to 15% between ranks
  6,000 and 16,000 in wAlbB, while the −10 box rate only goes from 47% to
  34% (background 8%). Single replicates are too shallow for weak sites, so
  a rule based on rank reproducibility stops at the depth limit, not where
  the evidence runs out.
* **Consistency between replicates** is reported per TSS (`consistency_p`:
  do the replicate counts fit the library sizes and the between-replicate
  dispersion?). Few sites are flagged (p < 0.001: 0–1.8%), but flagged sites
  look less like promoters (−10 box among sites with ≥ 30 reads: wAlbB 30%
  against 47%, *Salmonella* 20% against 33%; *E. coli* with its technical
  replicate 18% against 27%). Sites with few reads cannot be flagged (no
  power).

So with biological replicates the replicates agree within the measured
variation and the evidence lies in the pooled count; replicates are best
used to measure that variation and flag the rare site carried by one
library.

## Promoter score (October 2026)

The earlier measure counted TATAAT and TTGACA (≤ 1 mismatch) anywhere in
−15..−5 and −40..−28, separately. It is the textbook consensus but not the
current model of σ70 promoters, so a fixed score from the literature was
added (`harness.promoter_elements`; not learned from our calls, so it can
judge methods without circularity):

* **−10 core**: T, A and T at −12, −11 and −7, the hexamer starting at
  −14..−10. A(−11) and T(−7) are the bases σ reads, flipped out of the DNA
  (Feklistov & Darst 2011, *Cell* 147:1257, doi:10.1016/j.cell.2011.10.041).
* **−35 at the right distance**: TTGACA (≤ 2 mismatches) 16–18 nt upstream of
  a −10; some promoters lack a −35 and have a TG extended −10 instead
  (Shultzaberger et al. 2007, *Nucleic Acids Res* 35:771,
  doi:10.1093/nar/gkl956; Mitchell et al. 2003, *Nucleic Acids Res* 31:4689,
  doi:10.1093/nar/gkg694), reported as `ext10`.
* **+1 purine**.
* **Composition-matched background** (`harness.MatchedBackground`): random
  positions with the same AT share in −60..−1. AT-rich genomes contain more
  promoter-like words by chance; shifting sites by 100 nt does not correct
  for that.

Not included: the UP element (Ross et al. 1993, *Science* 262:1407,
doi:10.1126/science.8248780), the discriminator (Haugen et al. 2006, *Cell*
125:1069, doi:10.1016/j.cell.2006.04.034) and promoters of other σ factors
(review: Browning & Busby 2016, *Nat Rev Microbiol* 14:638,
doi:10.1038/nrmicro.2016.103); TSS of other σ factors count as misses.

**Check on our TSS** (`ettwiller2016` calls, background = shifted 100 nt):

| | *E. coli* TSS / background | wAlbB TSS / background |
|---|---|---|
| TATAAT ≤ 1 mismatch, −20..−5 (old style) | 28% / 6.6% | 58% / 20% |
| T·A···T core | 51% / 7.0% | 69% / 18% |
| −10 hexamer starts at −12 or −13 / −12 or −11 | 66% of hits | 55% of hits |
| −35 (≤ 2 mismatches) by spacer | peak at 16–17 nt | peak at 16–17 nt |
| Extended −10 TG | 10% / 5.7% | 15% / 8% |
| +1 purine | 81% / 50% | 88% / 53% |

The core finds more promoters than the hexamer and separates TSS from
background better (7× against 4× in *E. coli*); the −10 sits where the
structure places it, and the −35 is enriched at the expected spacing.

**Methods on the new score** (`results/promoter_real_*.tsv`; −10 core,
composition-matched background in brackets):

| | *E. coli* | *Salmonella* | wAlbB | wMelPop-CLA |
|---|---|---|---|---|
| `ettwiller2016` | 49.9% (7.4%) | 65.0% (8.7%) | 66.4% (17.3%) | 70.9% (17.7%) |
| `suliman2025` | 63.0% (8.5%) | 68.3% (8.5%) | 71.2% (19.3%) | 75.4% (21.2%) |
| beta-binomial (control) | 53.5% (8.1%) | – | – | – |
| local | 38.5% (6.7%) | 53.3% (7.5%) | 56.6% (16.8%) | 66.1% (16.5%) |

**Matched comparison** (`results/promoter_matched_*.tsv`; −10 core rate of
the top N; background 7–11% *E. coli* and *Salmonella*, 17–20% *Wolbachia*):

| Top N | *E. coli* RPM / local | *Salmonella* RPM / local | wAlbB RPM / local | wMelPop-CLA RPM / local |
|---|---|---|---|---|
| 1,000 | 70.5% / 74.7% | 70.3% / 70.8% | 72.5% / 72.9% | 73.7% / 74.1% |
| 4,000 | 59.6% / **70.9%** | 68.5% / 69.6% | 66.6% / 70.0% | 67.3% / 68.5% |
| 8,000 | 53.4% / **67.0%** | 66.4% / 67.8% | 62.5% / 66.3% | – |
| 16,000 | 47.1% / **60.8%** | 62.2% / 63.7% | 56.2% / 60.0% | – |

On the new score the conclusions are unchanged and clearer: the local test
ranks best everywhere (much better in *E. coli*), the spaced −35 agrees
(*E. coli* top 8,000: 17.9% against 15.2%), and in wAlbB the local test's FDR
stops at about 20,500 sites, after which the core rate falls to the
background (40%, 32%, 26%, 19% in successive bands of 2,000; background 20%).
The setting sensitivity is also unchanged (`results/promoter_sensitivity_*`).

## Promoter score trained on each genome (approach H, October 2026)

`call --genome`: positives are the strongest third of the calls,
negatives composition-matched random positions (three per positive, at most
4,000 positives); L2-regularised logistic regression on −60..−1 (bases and
neighbouring pairs); two-fold cross-validation over alternating 50-kb blocks,
so a TSS is never scored by a model that saw it. Standard library only.
Results: `results/promoter_model_*.tsv` and `.weights.txt`
(`promoter_model.py`, local-background calls).

**Recognition** (AUROC against composition-matched random positions;
fixed = the σ70 score above, combined as 2·core + TATAAT + spaced −35 +
0.5·extended −10):

| | E. coli | *Salmonella* | wAlbB | wMelPop-CLA |
|---|---|---|---|---|
| Strong calls (≥ 20 reads): genome model / fixed | 0.87 / 0.77 | 0.93 / 0.81 | 0.91 / 0.78 | 0.92 / 0.84 |
| Weak calls (< 10 reads): genome model / fixed | 0.75 / 0.62 | 0.89 / 0.72 | 0.86 / 0.70 | 0.88 / 0.78 |

Trained only on strong TSS, the model also recognises weak ones, and better
than the fixed score in every genome. A model that also used +1..+5 did
about as well in a prototype, so the gain is not 5'-ligation preference.

**What it learns.** The strongest positive single-base weights include A at
−11 and T at −7 in every genome (the bases σ reads; Feklistov & Darst 2011,
doi:10.1016/j.cell.2011.10.041). Consensus of positive weights:

| | −35 region | −10 region (−15..−5) | −1 |
|---|---|---|---|
| E. coli | ttggc at −36..−32 | gggtaAAATTt | C |
| *Salmonella* | tTgg | gGGTTAAATTTg | C |
| wAlbB | TtGgc | gGGTAAAATTt | C |
| wMelPop-CLA | TTGGc | GGGTAAAATTt | C |

**Rescue at low depth.** Reads thinned to 20%; calls and model from the
thinned data only; weak candidates there (2–5 pooled reads, local maximum)
are confirmed if the full data call them with ≥ 10 reads (AUROC; logistic
regression for the combination, 5-fold):

| | E. coli | *Salmonella* | wAlbB | wMelPop-CLA |
|---|---|---|---|---|
| Candidates (confirmed) | 42,418 (23%) | 14,766 (51%) | 10,904 (40%) | 2,230 (46%) |
| Read count | 0.718 | 0.765 | 0.741 | 0.763 |
| Promoter score | 0.718 | 0.580 | 0.663 | 0.580 |
| Both | **0.798** | 0.782 | **0.788** | 0.770 |
| Confirmed, lowest / highest score quarter | 9% / 44% | 41% / 60% | 20% / 56% | 36% / 54% |

**Orthogonal support** (wAlbB, all calls by score quarter): coverage rise in
the fragment library 29%, 31%, 32%, 39% (background 12%).

So the sequence score carries information that read counts do not, most in
the deepest-sampled genomes with many weak sites (E. coli, wAlbB): among weak
candidates, those in the top score quarter are confirmed 3–5 times as often
as those in the bottom quarter in E. coli and wAlbB. It is reported per TSS;
using it to filter or rescue calls is left to the standard mode (plan
section 8).

## Peak shape (approach G, October 2026)

**Position.** The reported position of a cluster (the strongest base) is
already precise: it is a purine 81–88% of the time against 28–49% one base
either side, and the −10 core starts at −12/−11 as expected, for
`ettwiller2016` and the local test on *E. coli* and wAlbB. Nothing to change.

**Shoulders.** Calls 7–30 nt downstream of a much stronger call on the same
strand (`cappable_tss.mark_shoulders`; `peak_ratio` = RPM / RPM of the
strongest call within 30 nt upstream). Removed calls compared with kept
calls of the same RPM (`shoulders.py`, `results/shoulders_*.tsv`; window
30 nt, ratio 0.1; −10 core rate):

| | Removed | Removed | Same-RPM kept | Background | Positive genome score, removed / same-RPM |
|---|---|---|---|---|---|
| E. coli `ettwiller2016` | 732 (4.4%) | 13% | 36% | 10% | 16% / 33% |
| E. coli local | 1,752 (3.4%) | 13% | 38% | 8% | 21% / 54% |
| wAlbB local | 698 (3.4%) | 44% | 44% | 20% | 44% / 59% |
| *Salmonella* `ettwiller2016` | 159 (1.5%) | 55% | 66% | 8% | 54% / 64% |
| *Salmonella* local | 1,034 (2.6%) | 40% | 48% | 8% | 50% / 67% |
| wMelPop-CLA `ettwiller2016` | 42 (2.0%) | 57% | 74% | 33% | 26% / 50% |

In *E. coli*, shoulders have background-level promoter evidence and removing
them costs only Table S1 sites that the same pipeline called (exact-base
recall 0.959 → 0.920). In *Salmonella* and both *Wolbachia* strains most
shoulders have promoters (tandem promoters); removing them would lose real
TSS. The data do not show why *E. coli* differs. Shoulder removal is
therefore an option (`--shoulder-ratio`, off by default), and `peak_ratio`
is reported for every TSS. The coverage-rise check cannot judge shoulders
(a site just downstream of a strong TSS always sees rising coverage) and is
not used here.

## Standard mode (October 2026)

`call --mode standard` (the default) chooses from the design
(`cappable_tss.choose_design`): the local-background test with ≥ 2 enriched
replicates, the Ettwiller cut-offs with one library; a control removes
depleted sites; enriched libraries failing the library check stop the run.
Each TSS gets a confidence (`assign_confidence`): **high** = ≥ 10 reads,
positive genome-trained promoter score (with `--genome`), consistent
replicates (consistency_p ≥ 0.001) and not a shoulder (peak_ratio ≥ 0.1);
**medium** = ≥ 10 reads or a positive score, and consistent; **low** = the
rest. The thresholds come from the earlier sections (10 reads: weak sites
were usually confirmed above it; 0.1: shoulder analysis).

Validation (`standard_mode.py`, `results/standard_mode_*.tsv`, with
`--genome`): *confirmed* = share of the TSS called from 20% of the reads
(in that tier) that the full data call with ≥ 10 reads; −10 core rate
(partly circular, as the tiers use the promoter score); coverage rise in the
fragment library (wAlbB):

| | Design | TSS: high / medium / low | Confirmed: high / medium / low | −10 core: high / medium / low (background) |
|---|---|---|---|---|
| E. coli | one library + control: cut-offs | 8,962 / 7,677 / 0 | 99% / 75% / 47% | 68% / 28% / – (7%) |
| *Salmonella* | 2 replicates: local | 13,831 / 17,752 / 8,590 | 99.5% / 89% / 79% | 67% / 51% / 36% (7%) |
| wAlbB | 3 replicates: local | 7,522 / 8,619 / 4,344 | 99% / 89% / 79% | 67% / 54% / 45% (17%) |
| wMelPop-CLA | 3 replicates: local | 1,780 / 2,543 / 1,456 | 100% / 88% / 79% | 73% / 65% / 60% (18%) |

wAlbB coverage rise: high 42%, medium 31%, low 22% (background 12–14%).
The tiers are ordered on every measure. In E. coli no full-data call is low
because the 1.5-RPM cut-off already implies ≥ 10 reads there. Without
`--genome` the tiers separate less (wAlbB: 11,768 high with a 62% −10 core
rate, 351 medium, 8,366 low), so `--genome` is recommended. "Low" is
relative: in the replicated datasets most low-confidence TSS are still
confirmed by deeper data.

## Several bacteria in one sample (multi-genome branch, October 2026)

`call --genome A.fa B.fa ...` splits each library by organism and calls each
separately. Test (`mixture/`, `results/mixture_ecoli_salmonella.tsv`):
*E. coli* reads mixed into two *Salmonella* replicates at 100%, 10% and 2%
of its reads, mapped to both genomes together, against each organism's same
reads mapped to its own genome and called alone. *Recovered*: share of the
TSS called alone that the mixture also calls (± 2 nt); *also alone*: share
of the mixture's TSS that were also called alone.

| *E. coli* share | MAPQ | *E. coli* reads, mixed / alone | *E. coli* TSS, mixed / alone | Recovered | Also alone | High-confidence also alone |
|---|---|---|---|---|---|---|
| 2% | 20 | 325,179 / 361,339 | 5,176 / 5,307 | 96.6% | 99.1% | 99.9% |
| 2% | 0 | 810,739 / 361,339 | 5,781 / 5,307 | 99.0% | 90.9% | 91.2% |
| 10% | 20 | 1.63 M / 1.81 M | 15,365 / 15,611 | 97.5% | 99.0% | 99.9% |
| 10% | 0 | 2.29 M / 1.81 M | 16,149 / 15,611 | 99.1% | 95.8% | 96.3% |
| 100% | 20 | 16.3 M / 18.1 M | 68,281 / 68,717 | 98.4% | 99.0% | 99.6% |

*Salmonella* (the abundant organism) at MAPQ 20: 98.3% recovered, 99.7%
also called alone, in every mixture.

* At MAPQ 20 the rare organism's calls are almost unchanged by the
  abundant one: with *E. coli* at 2%, 99.9% of its high-confidence TSS were
  also called alone.
* Without the filter (MAPQ 0), the abundant organism's reads leak into the
  rare one: *E. coli* at 2% gains about 450,000 reads, more than its own,
  and 526 calls not made alone (9%), 46% of them within 50 nt of rRNA or
  tRNA genes (0.4% of normal calls) and 139 rated high confidence. At
  MAPQ 20 there are 48 such calls (1 high).
* The filter costs about 10% of reads (those that fit both genomes) and
  2–3% of TSS, mostly in conserved regions.

So: map to all genomes and the host together, use `--min-mapq 20` (the tool
warns below 10), and call each organism separately.
