# Plan: developing and benchmarking TSS-calling approaches

Goal: improve `cappable_tss.py call` beyond the fixed cut-offs of Ettwiller
et al. 2016 and Suliman et al. 2025, especially for symbiont genomes
(*Wolbachia*, *Buchnera*), and show the gain on real data before any change
becomes a default. `--preset ettwiller2016` stays as the fixed baseline.

## 1. Datasets

- [x] *E. coli* MG1655, Ettwiller et al. 2016 (PRJEB9717): enriched, technical
      replicate, non-enriched control; published TSS table (Table S1).
- [x] *Salmonella* SL1344 (PRJEB58986): 2 enriched replicates, no control.
- [x] *Wolbachia* wAlbB, Suliman et al. 2025 (PRJNA1145481): untreated 28 °C
      triplicate mapped; 33 further libraries (temperature, antibiotics) available.
- [x] *Wolbachia* wMelPop-CLA, Suliman et al. 2025: 28 °C 48 h triplicate
      mapped; low depth (1–2 M *Wolbachia* reads per library).
- [x] *Wolbachia* wAlbB PRJNA1392518 (Frontiers in Microbiology 2026,
      doi:10.3389/fmicb.2026.1860204): same Aa23 cell line as Suliman et al.,
      but whole primary transcripts captured, then heat-fragmented and
      sequenced with NEBNext Ultra II Directional (random-primed, stranded,
      read 2 = sense, 2 × 150 bp). Three experiments: unenriched / 1 / 2
      rounds of enrichment; unenriched / XRN-1 / XRN-1 + Cappable-seq /
      Cappable-seq; days 3–7 time course in triplicate. Not 5'-end data.
      Cappable-seq and unenriched libraries of the XRN-1 experiment mapped.
- [ ] Use PRJNA1392518 as an independent check of wAlbB TSS: coverage rises
      (enriched vs unenriched) at TSS called from the Suliman 5'-end data.
- [x] Run TSSpredator/ANNOgesic on PRJNA1392518 (its intended data type) and
      ANNOgesic transcript, terminator, UTR and operon modules
      (`benchmarks/frontiers2026/`). TSSpredator: 89% of its 250 TSS within
      25 nt of 5'-end TSS, but a median 12 nt downstream; 14% recall.
- [ ] *Buchnera* (own data): obtain reads, genome and annotation; record
      strain, host, library kit, replicates and whether a non-enriched control
      exists (see section 8).
- [ ] Suliman et al. unenriched libraries (Table 1 of the paper): locate or
      deposit; they would give *Wolbachia* a real control.
- [x] Write a dataset manifest (`benchmarks/datasets.toml`): accessions,
      genome and annotation accessions, host, library type, replicate
      structure, adapters, read layout, which mate carries the RNA 5' end.
- [ ] One script per dataset that rebuilds the BAMs from the manifest
      (as `benchmarks/ettwiller2016/run.sh` does now).

## 2. Evaluation measures

No dataset except *E. coli* has a published TSS list, so every approach is
scored on several independent lines of evidence:

- [ ] **Reference agreement** (*E. coli* only): recall and precision against
      Table S1 at 0, ±1, ±5 nt; report separately for strong (≥ 20 reads) and
      weak TSS, since Table S1 has a 1.5 RPM floor.
- [x] **Promoter evidence**: share of TSS with a −10 / −35 box, per class,
      against the same positions shifted 100 nt. Use a motif learned from each
      genome's own confident TSS (not *E. coli* TATAAT), because AT-rich
      genomes inflate fixed-word matches. Done as a fixed σ70 score from the
      literature (−10 core, spaced −35, extended −10, +1 purine) with a
      composition-matched background (benchmarks/README.md, "Promoter
      score"); the genome-learned motif is approach H.
- [ ] **Orthogonal support** (wAlbB): share of TSS with a coverage rise
      (≥ 1.5×, 50 nt downstream vs upstream) in the PRJNA1392518 enriched
      library, against shifted positions. First result: 73.1% of TSS called
      with the Suliman settings, 48.6% of the extra TSS called by the
      Ettwiller preset, 6.7% of shifted positions. (Enrichment over the
      unenriched library does not discriminate: 97% of shifted positions pass.)
- [ ] **Reproducibility**: leave-one-replicate-out — call on n−1 replicates
      and measure how many calls the held-out replicate supports.
- [ ] **Depth robustness**: downsample to 50%, 25% and 10% of reads and
      measure recovery of the full-depth calls; fixed RPM cut-offs should fail
      this, statistical calls should degrade gracefully.
- [ ] **Base-level precision**: offset distribution to the reference or to
      the replicate consensus.
- [ ] **Sanity checks**: +1 purine and −1/+1 YR frequency; TSS on sequences
      that should have none (e.g. the missing SL1344 plasmid pRSF1010); share
      of calls at known processing sites (tRNA/rRNA ends) as false positives.
- [x] **Simulation**: synthetic libraries with planted TSS of known strength
      on each real genome (background drawn from the real data), giving exact
      precision and recall per method and depth.
- [x] Benchmark harness (`benchmarks/harness.py`): one command that runs every
      method on every dataset and writes one summary table plus the report.

## 3. Baselines

- [x] `--preset ettwiller2016` (RRS ≥ 1.5, enrichment ≥ 0, 6 nt clusters).
- [x] Suliman et al. settings (MAPQ ≥ 20, ≥ 10 CPM, 10 nt clusters, all
      replicates).
- [x] TSSpredator via ANNOgesic 1.1.14, default parameters (*E. coli*).
- [ ] TSSpredator after ANNOgesic `optimize_tss_ps` on a manually checked set
      (a fair version of the coverage-step approach).

## 4. Approaches to develop

Each one is an option that stays off by default until it passes sections 5 and 6.

- [x] **A. Library-type check**: share of reads at the 1,000 busiest start
      positions, strand of reads inside genes, insert size; warn when a
      library does not look like 5'-end data (would have flagged PRJNA1392518).
- [x] **B. Normalisation**: size factors from reads outside rRNA and tRNA
      genes, or median-of-ratios over shared positions, so log2 enrichment 0
      means "no different from the control".
- [x] **C. Test against the control**: one-sided binomial test of enriched
      against control counts at each position, with an optional minimum fold
      change.
- [x] **D. Test against local background** (no control): Poisson or
      negative-binomial test of the position against read starts in the
      surrounding ±50–200 nt, excluding the cluster itself.
- [x] **E. False-discovery rate**: Benjamini–Hochberg q-values over tested
      positions (positions with ≥ 2–3 reads).
- [x] **F. Replicates as evidence**: combine per-replicate p-values (Fisher)
      with a direction check; report a reproducibility score per TSS.
      Result: Fisher is invalid for the local test (background repeats
      between replicates); pooling, reads in one replicate, and a per-TSS
      consistency p-value instead (benchmarks/README.md).
- [x] **G. Peak shape**: local-peak ratio and cluster spread, so the
      strongest position of an imprecise promoter is reported once and weak
      shoulders are not called separately. Result: positions already precise;
      `peak_ratio` per TSS and `--shoulder-ratio` (off): shoulders are
      background-like in E. coli but mostly promoters elsewhere.
- [x] **H. Genome-specific promoter score**: self-trained sequence model
      (labels from enriched vs depleted positions, trained per genome with
      block cross-validation) reported per TSS and optionally used to rescue
      borderline calls. Result: `call --genome` (labels: strongest third of
      calls against composition-matched random positions); learns −35/−10
      in every genome, beats the fixed σ70 score, improves prediction of
      which weak sites deeper sequencing confirms (benchmarks/README.md).
      Use for filtering/rescue: decide in the standard mode.
- [ ] **I. Negative-binomial model across replicates and conditions**
      (DESeq2/edgeR-like), only if F is not enough; adds heavier dependencies.

## 5. Testing the statistical approaches (C, D, E, F)

A statistical call is only useful if its p-values and false-discovery rates
mean what they say. Each test below must pass before the approach is judged
on biology in section 6.

**Correctness of the arithmetic**

- [x] Unit tests: binomial and Poisson tail probabilities, Fisher's combined
      p-value and Benjamini–Hochberg q-values against reference values
      computed independently (e.g. with SciPy in the test suite only).
- [ ] Edge cases: zero control reads, zero background, a single replicate,
      very deep positions (no underflow), ties in q-values.

**Calibration under the null (no false discoveries expected)**

- [x] (C done; D to do) Simulated null: libraries with no planted TSS. p-values must be close
      to uniform (QQ plot) and the FDR at 1%, 5%, 10% must hold.
- [x] Real-data null for test C: split one library into two halves by random
      read sampling (control vs control, and enriched vs enriched); calls at
      a given FDR should be about that fraction of tests, not thousands.
- [x] Real-data null for test C: swap enriched and control labels (*E. coli*);
      the "TSS" called should be processed ends (tRNA/rRNA), not promoters.
- [ ] Real-data null for test D: positions shifted 100 nt and randomly placed
      positions with the same read count; the call rate at a given FDR should
      be near the nominal level.
- [x] Dispersion check: across the Suliman triplicates, compare variance with
      mean per position; if the variance clearly exceeds the mean, the
      Poisson tests are anti-conservative and approach F/I must use a
      negative-binomial model instead.

**Power and robustness**

- [x] (C done; D to do) Simulation with planted TSS of known strength at several depths
      (including *Wolbachia*-like 1–2 M reads): recall at fixed FDR compared
      with the fixed cut-offs, as power curves by TSS strength.
- [x] (C done; D to do) Downsampling the real libraries (50%, 25%, 10%): statistical calls
      should lose weak TSS gradually; fixed RPM cut-offs change meaning with
      depth.
- [x] Parameter sensitivity: local-background window (±50, 100, 200 nt),
      minimum read count (2, 3, 5), FDR level; results should not hinge on one
      setting.

**Matched comparison with the fixed cut-offs**

- [x] For the same number of calls, compare statistical ranking with RPM
      ranking on the independent evidence: −10 box rate, orthogonal coverage
      rise in PRJNA1392518 (wAlbB), *E. coli* Table S1 recall, held-out
      replicate support. This separates "better ranking" from "more calls".
- [ ] For each dataset, report what the fixed cut-off and the FDR call each
      gain and lose, with the independent evidence for both sets.

## 6. How each approach is judged

For every approach and every dataset, compared with the `ettwiller2016` and
Suliman baselines:

- [ ] TSS gained: promoter-box rate close to that of TSS found by all methods.
- [ ] TSS lost: lower promoter-box rate and poorer reproducibility than kept TSS.
- [ ] *E. coli*: recall of strong Table S1 TSS (≥ 20 reads) not lower by more
      than 2 percentage points.
- [ ] Leave-one-replicate-out support not lower; depth robustness better.
- [ ] Simulation: precision and recall at least as good at every depth.
- [ ] Runtime and memory acceptable for a bacterial genome on a laptop.

An approach becomes a default only if it passes on *E. coli* and at least
two of the symbiont datasets, and the result is written into the report.

## 7. Development steps

- [x] Write the dataset manifest and the benchmark harness (sections 1–2);
      baselines in `benchmarks/README.md`.
- [x] Add simulation of planted TSS (section 2), calibrated from real libraries.
- [x] Implement A (library-type check): `qc` subcommand and checks in `call`;
      flags all fragment libraries (22–27% peak share) and passes all 5'-end
      libraries (77–97%; 42% for the E. coli control).
- [x] Implement B + C + E (control-based statistics); run the section 5
      tests, then benchmark on *E. coli*. Result (benchmarks/README.md):
      calibrated with a beta-binomial (dispersion 0.01–0.05 from replicates),
      but ranks no better than RPM at matched size and has low power at low
      depth; kept as `--method binomial`, not the default. Enrichment over the
      control is best used as a filter; D is the main test to develop next.
- [x] Implement D + E (local background); run the section 5 tests, then
      benchmark on *Salmonella*, wAlbB, wMelPop-CLA, *Buchnera*. Result
      (benchmarks/README.md): `--method local`; ranks at least as well as RPM
      everywhere (much better in E. coli), 2–3.5× more TSS, mostly weak; a
      negative binomial does not fit the uneven background, so the dispersion
      is a setting and q-values are nominal; needs replicates. Buchnera to do.
      Open: the real-data null for D (no clean null exists in 5'-end data:
      the background is reproducible) and the choice of stop.
- [x] Implement F (replicates); leave-one-out benchmark on the triplicates.
- [x] Implement G (peak shape); benchmark base-level precision.
- [x] Turn `selftrain.py` into H (`learn` subcommand); benchmark promoter
      scores per genome, including *Buchnera*. Done as `call --genome`
      (no separate subcommand); *Buchnera* still to do.
- [ ] Decide new defaults from sections 5 and 6; keep `--preset ettwiller2016`.
- [ ] Unit tests for every new statistic (exact values on small inputs).
- [ ] Update README, benchmark READMEs and the report after each step.

## 8. Deliverables

The end products of this plan:

- [x] **Ranked experimental-design recommendations** (docs/design_recommendations.md): which designs give
      reliable TSS, ranked on the benchmarks (control or not; number of
      biological replicates; technical against biological replicates;
      sequencing depth, using downsampling; library type, 5'-ligation against
      fragment libraries; spike-ins or UMIs if data exist). Each
      recommendation states the evidence (dataset, measure) and what is
      lost without it. Evidence so far:
      - 5'-end (adapter-ligation) libraries are required; fragment libraries
        do not mark TSS (section A; Frontiers 2026 data).
      - Biological replicates are worth more than one deep library: they
        let the local test stop sensibly and measure between-library
        variation; with one library, read stacks are called (D, F).
      - A non-enriched control is what removes processed ends: without one,
        no method separates them from TSS in the simulation (D).
      - Depth: the RPM cut-off calls more noise at lower depth; the
        statistical tests lose weak TSS gradually (C, D).
- [x] **One tool with a standard mode and a development mode.**
      Done: `call --mode standard` (default; design from replicates,
      control and library check; confidence tiers validated in
      benchmarks/README.md) and `--mode dev` (all methods and settings).
      - Standard mode: `call --design auto` (or a named design) chooses the
        method and settings from the inputs (replicates, control, library
        check) using the benchmark results, and reports per-TSS confidence
        (q-value, replicate support, consistency, later the promoter score).
      - Development mode: every method and setting exposed, as now, plus the
        benchmark tools, for new designs and for re-tuning on new organisms.
      - Both modes documented with the design table they rely on.

## 9. Open questions

- [ ] *Buchnera*: species/strain and genome accession, host aphid, library
      kit (5'-ligation Cappable-seq or fragment library), read layout,
      replicates, conditions, and whether a non-enriched control exists.
- [ ] Are the Suliman et al. unenriched libraries available?
- [ ] Where should the full benchmark run? This cloud session has about 30 GB
      of disk and 2-hour background jobs; a cluster run of the harness would
      allow all 63 Suliman libraries.
- [ ] Acceptance thresholds in section 6: agree or adjust.
