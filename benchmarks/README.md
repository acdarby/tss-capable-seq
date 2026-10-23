# Benchmarks

Tools to judge TSS-calling methods on real and simulated Cappable-seq data
(plan: `../docs/benchmark_plan.md`).

| File | Purpose |
|---|---|
| `datasets.toml` | Dataset manifest: files (relative to a data directory), accessions, replicate structure, mapping quality, which mate carries the RNA 5' end, optional reference TSS and orthogonal data |
| `harness.py` | Runs every method on every dataset and writes one table of measures |
| `simulate.py` | Synthetic 5'-end libraries with planted TSS, optionally calibrated from a real library |
| `ettwiller2016/` | Rebuilds the *E. coli* data and compares with the published TSS table |
| `frontiers2026/` | TSSpredator and ANNOgesic on wAlbB fragment-library data |
| `tools/` | Download SRA runs, list runs per project, annotate genomes, TSS statistics |
| `results/` | Baseline results (TSV) |

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
