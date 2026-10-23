# Experimental design for Cappable-seq TSS mapping: recommendations

Ranked by how much each choice changed the TSS that could be found on the
public data benchmarked here (`benchmarks/README.md`): *E. coli* MG1655
(Ettwiller et al. 2016, doi:10.1186/s12864-016-2539-z), *Salmonella* SL1344
(PRJEB58986), *Wolbachia* wAlbB and wMelPop-CLA in mosquito cells (Suliman et
al. 2025, doi:10.1099/mgen.0.001542) and wAlbB fragment libraries (Frontiers
in Microbiology 2026, doi:10.3389/fmicb.2026.1860204). Each recommendation
gives the evidence and what is lost without it. Where a point is an
inference that these data did not test, it says so.

"High-confidence TSS" are those the standard mode (`call`, with `--genome`)
rates high: ≥ 10 reads, a positive genome-trained promoter score, consistent
replicates and not a shoulder. When TSS called from 20% of the reads are
checked against the full data, 99–100% of high-confidence ones are
confirmed.

## Summary

| Rank | Recommendation | Effect on the benchmarks |
|---|---|---|
| 1 | 5'-end libraries (adapter ligation to the RNA 5' end) | Essential: fragment libraries do not mark TSS |
| 2 | Independent replicate libraries (≥ 2, better 3) instead of one deeper library | 1.1–2.2× more high-confidence TSS for the same reads; the statistical test needs them |
| 3 | Sequence deeper (count reads on the bacterial genome) | No saturation up to 16–32 M reads; about +1,500–2,000 high-confidence TSS per doubling |
| 4 | A non-enriched control | Removes processed RNA ends: nearly all of them in the simulation; in *E. coli* about 2% of calls, mostly at rRNA/tRNA |
| 5 | Genome sequence and annotation for the analysis | Promoter score, rRNA/tRNA-aware normalisation, gene classes; sharper confidence tiers |
| 6 | Read processing: right mate, mapping quality, host reads, duplicates | Wrong mate or host reads ruin calls; duplicates (UMIs) untested here |

## 1. Make 5'-end libraries

The 5' end of each read must be the 5' end of the RNA: adapters ligated to
the RNA 5' end after capture, as in the original Cappable-seq protocol
(NEBNext small RNA library prep in Ettwiller et al. 2016 and Suliman et al.
2025).

*Evidence.* In the wAlbB fragment libraries (random-primed, directional
RNA-seq of the captured RNA), only 22–27% of reads start at sharp peaks,
against 77–97% in every 5'-end library, so read starts do not mark TSS.
TSSpredator's coverage steps on these data placed TSS a median 12 nt
downstream of the 5'-end TSS and recovered 14% of them
(`benchmarks/frontiers2026/`). The library check (`qc`, and the first step
of `call`) flags such libraries; the standard mode stops on them.

*Without it:* TSS can only be approximated from coverage rises, at lower
resolution and sensitivity.

## 2. Independent replicate libraries rather than one deeper library

Prefer two, better three, independent libraries (biological replicates:
separate cultures or cell batches) to one library sequenced to the same
total depth.

*Evidence.* Reads of the replicates pooled into one library (same reads, no
replicates) against the replicates as given, standard mode
(`benchmarks/design.py`, `results/design_*.tsv`):

| | Replicates | High-confidence TSS, replicates / one library | All TSS, replicates / one library |
|---|---|---|---|
| wAlbB | 3 | 7,522 / 3,384 | 20,485 / 5,436 |
| *Salmonella* | 2 | 13,831 / 8,160 | 40,173 / 12,036 |
| wMelPop-CLA | 3 | 1,780 / 1,646 | 5,779 / 2,975 |

With replicates the tool can use the test against the local background,
which ranked TSS best on every dataset and stopped where calls stopped
looking like promoters (approach D). With one library that test also calls
isolated read stacks (51,000 calls in *E. coli*), so the tool falls back to
fixed cut-offs, whose meaning changes with depth (*E. coli*: 16,639 TSS at
full depth, 28,374 at 10% of the reads). Replicates also measure the
variation between libraries (dispersion 0.01–0.05 here) and flag sites
carried by one library (`consistency_p`).

Weak sites missing from some replicates were as promoter-like as those
seen in all (approach F), so replicates need not each be deep; their value
is independence, not agreement at every site.

*Biological or technical?* Not tested directly. Resequencing the same
library adds depth but cannot reveal artefacts made in that library, such
as PCR-duplicate stacks; separate library preparations can. (The *E. coli*
"technical replicate" ERR930222 varied from its partner about as much as
separate libraries do, dispersion 0.046, so it is probably a separate
library preparation.)

*Without it:* fewer TSS, no statistical stop, no replicate consistency per
TSS.

## 3. Sequence deep enough, counting reads on the bacterial genome

Plan depth as 5'-end reads on the bacterial genome, not total reads: in
symbiont samples most reads can come from the host.

*Evidence.* Standard mode on thinned data (`benchmarks/design.py`):

| Reads on the genome | wAlbB high / all TSS | | Reads on the genome | *Salmonella* high / all TSS |
|---|---|---|---|---|
| 3.2 M | 2,401 / 7,835 | | 1.6 M | 3,892 / 12,340 |
| 8.0 M | 4,008 / 12,077 | | 4.0 M | 6,628 / 20,484 |
| 15.9 M | 5,497 / 15,993 | | 8.0 M | 9,769 / 28,996 |
| 31.9 M | 7,522 / 20,485 | | 16.0 M | 13,831 / 40,173 |

Neither curve has levelled off: each doubling still adds high-confidence
TSS, mostly weaker ones (internal and antisense TSS). With the
statistical test, fewer reads give fewer but still reliable TSS (calls from
10% of the wAlbB reads were 96% confirmed by the full data); with a fixed
cut-off and one library, fewer reads give more noise (see 2).

*Without it:* weak TSS are missed; strong ones are found even at 1.6–3 M
reads.

## 4. Include a non-enriched control

Sequence a non-enriched library of the same RNA, prepared the same way. Its
role is to remove processed RNA ends (5'-monophosphate RNA from rRNA and
tRNA maturation and from cleavage), which survive the enrichment at low
levels.

*Evidence.* In *E. coli*, the same calls without the control: 17,005 instead
of 16,639 TSS, and the share within 5 nt of rRNA or tRNA genes rises from
0.4% to 2.3%; the promoter rate is unchanged (`results/design_ecoli.tsv`).
In the simulation (wAlbB-like), processed ends planted at 5% of their
control strength were called 96–97% of the time at full depth and 72–78% at
a quarter of the reads without a control (4% at 5% of the reads), and ≤ 1%
with one (`results/sim_walbb_local_bias0.tsv`). Used as
a statistical test of enrichment, the control ranked TSS no better than
read count (approach C), so its value is as a filter.

*Inference, not tested:* the *E. coli* control was 76% rRNA reads; depleting
rRNA from the control would give more informative control reads for the
same sequencing.

*Without it:* processed ends cannot be told from TSS; the effect was small in
*E. coli*, but may be larger where enrichment is less efficient. The tool
notes it.

## 5. Provide the genome and annotation to the analysis

Give `--genome` (FASTA) and `--annotation` (GFF3, with rRNA, tRNA and CDS
features).

*Evidence.* The genome lets the tool learn the genome's own promoter (it
found the σ70 −35 and −10 in all four genomes, without *E. coli* training)
and score every TSS; with it, weak candidates later confirmed by deeper data
were ranked better (AUROC 0.72 → 0.80 in *E. coli*, 0.74 → 0.79 in wAlbB)
and the confidence tiers separate better (approach H, standard mode). The
annotation gives the gene classes, the strand check of the library test and
library sizes without rRNA/tRNA reads (approach B).

## 6. Read processing

* **Mate.** For paired-end data count the mate whose 5' end is the RNA 5' end
  (`--mate`); the library check warns when fewer than 30% of reads in genes
  are on the gene's strand.
* **Host reads.** For symbionts, map to the symbiont and host (or host
  mitochondrion) together and filter on mapping quality (Suliman et al.
  used MAPQ ≥ 20; the *Wolbachia* benchmarks here mapped to the symbiont plus
  the host mitochondrion with `--min-mapq 20`).
* **Adapters.** Trim adapters and drop reads that become very short; on the
  public data, empty reads after trimming slowed mapping greatly.
* **Duplicates.** *Not tested here:* unique molecular identifiers would let
  PCR duplicates be removed, which matters most with a single library. No
  benchmarked dataset had UMIs.

## Designs and what the tool does with them

| Design | Standard mode | Main limitation |
|---|---|---|
| ≥ 2 replicate 5'-end libraries + control | Local-background test, control as filter | – (recommended) |
| ≥ 2 replicate 5'-end libraries, no control | Local-background test | Processed ends not removable |
| 1 library + control | Ettwiller et al. cut-offs | Fixed cut-off depends on depth; no replicate checks |
| 1 library, no control | RPM cut-off only | Both of the above |
| Fragment (random-primed) libraries | Stops | Read starts do not mark TSS |

Still to test: *Buchnera* (own data), and whether Suliman et al. sequenced
non-enriched libraries for the *Wolbachia* samples.
