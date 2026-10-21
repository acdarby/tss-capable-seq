# TSS calling for Cappable-seq

`cappable_tss.py` is a standalone tool that finds transcription start
sites (TSS) in bacterial Cappable-seq data and annotates them in the style
of ANNOgesic (Primary, Secondary, Internal, Antisense, Orphan). ANNOgesic's
own TSS step relies on TSSpredator, which was designed for dRNA-seq;
this tool works directly on Cappable-seq read 5' ends instead. It can also
prepare Cappable-seq data for the original ANNOgesic pipeline.

## Why not use ANNOgesic's TSS step as it is?

ANNOgesic's `tss_ptss` module runs TSSpredator (Dugar et al. 2013). For
dRNA-seq, TSSpredator compares a library treated with terminator
exonuclease (TEX+), which is enriched for primary transcripts, with an
untreated library (TEX-), and looks for sharp rises ("steps") in read
coverage. Cappable-seq instead captures 5'-triphosphorylated (primary)
transcripts directly, so the 5' end of each read in the enriched library
marks the TSS at single-base resolution (Ettwiller et al. 2016). The
natural signal is therefore the **number of read 5' ends at each
position**, compared with a non-enriched control, rather than a coverage
step.

## Requirements

* Python 3 (standard library only).
* `pysam`, only if you give BAM files (`pip install pysam`).

## Commands

### `call`: find and classify TSS

1. Counts read 5' ends per position and strand (BAM, BED or 5'-end wiggle).
2. Normalises each library to reads per million (RPM), using the total
   number of 5' ends counted in that library.
3. Keeps a position in a replicate if the enriched RPM is at least
   `--min-rpm` and `log2((enriched + pc) / (control + pc))` is at least
   `--min-log2fc` (`pc` is `--pseudocount`). A position must pass in at
   least `--min-replicates` replicates. Reported values are means over
   replicates.
4. Merges calls on the same strand that lie no more than
   `--cluster-distance` nt apart and keeps the strongest one.
5. If `--annotation` is given, classifies each TSS against the GFF3
   features chosen with `--feature-types` (default `gene`; `CDS` is often
   better because gene features can include UTRs). The categories follow
   Sharma et al. 2010, as in ANNOgesic and TSSpredator:
   * **Primary**: within `--utr-length` (300) nt upstream of a gene on
     the same strand, and the strongest such TSS for that gene.
   * **Secondary**: within that window but not the strongest.
   * **Internal**: inside a gene, on the same strand.
   * **Antisense**: inside a gene or within `--antisense-distance` (100)
     nt of it, on the opposite strand.
   * **Orphan**: none of the above.

   A TSS can have more than one class (for example `Primary&Internal`).
   Without `--annotation`, the class is reported as `NA`.
6. Writes `<prefix>.gff` and `<prefix>.tsv`, creating the output folder
   if needed, and prints a count per class. The GFF uses ANNOgesic-style
   attributes (`type`, `associated_gene`, `UTR_length`) plus
   `enriched_rpm`, `control_rpm`, `log2_enrichment` and `replicates`.

```bash
python cappable_tss.py call \
    --enriched enr_rep1.bam enr_rep2.bam \
    --control  ctl_rep1.bam ctl_rep2.bam \
    --annotation genome.gff --feature-types CDS \
    --min-replicates 2 \
    --out-prefix results/TSS
```

Default thresholds: `--min-rpm 1`, `--min-log2fc 1`, `--pseudocount 1`,
`--min-replicates 1`, `--cluster-distance 5`. These are starting points,
not values taken from a paper. Tune them for your sequencing depth and
enrichment efficiency; checking a few known promoters in a genome browser
is the quickest test.

Input notes:

* **BAM**: unmapped, secondary, supplementary and QC-fail reads are
  skipped; use `--min-mapq` to drop low-quality alignments. Duplicates
  are not removed. For paired-end data only the mate given by `--mate`
  (default 1) is counted, and its 5' end must be the RNA 5' end, so check
  your library protocol.
* **BED**: needs the strand in column 6 (`bedtools bamtobed` output works).
* **Wiggle**: given as `forward.wig,reverse.wig`, and must contain 5'-end
  counts, not coverage. Negative reverse-strand values (READemption and
  ANNOgesic style) are fine.
* **Controls** are paired with enriched libraries in the order given.
  `--control` can be left out; then only the RPM cut-off applies, which is
  much less specific because processed RNA ends are not filtered out.

### `annogesic-libs`: prepare Cappable-seq data for ANNOgesic

Use this if you want to stay inside the ANNOgesic workflow, for example to
use its TSSpredator parameter optimisation or its downstream modules. It
writes forward and reverse wiggle files for each library and prints the
strings to pass to `annogesic tss_ptss --tex_notex_libs`. The enriched
library is labelled `tex` and the control `notex`, so TSSpredator compares
enriched with control:

```bash
python cappable_tss.py annogesic-libs \
    --enriched enr_rep1.bam enr_rep2.bam \
    --control  ctl_rep1.bam ctl_rep2.bam \
    --outdir ANNOgesic/input/wigs/tex_notex
```

* Input must be BAM or BED, with one control per enriched library.
* By default the wiggle files hold read coverage, which is the signal
  TSSpredator was designed for. `--signal 5prime` writes 5'-end counts
  instead.
* Reverse-strand values are written as negative numbers; use
  `--positive-reverse` to turn this off.
* Replicates are labelled `a`, `b`, `c`… in input order, under condition
  `--condition` (default 1).
* TSSpredator's default parameters are tuned for dRNA-seq, so re-optimise
  them for Cappable-seq (ANNOgesic's `optimize_tss_ptss`, using a set of
  TSS you have checked by hand).
* Check the printed library strings against `annogesic tss_ptss --help`
  for your ANNOgesic version.

## Tests

```bash
python -m unittest discover -s tests -v
```

The BAM test is skipped if `pysam` is not installed.

## Licence

MIT (see `LICENSE`). The tool does not contain or depend on ANNOgesic or
TSSpredator code.

## References

* Yu S-H, Vogel J, Förstner KU (2018). ANNOgesic: a Swiss army knife for
  the RNA-seq based annotation of bacterial/archaeal genomes.
  *GigaScience* 7(9). https://doi.org/10.1093/gigascience/giy096
* Dugar G *et al.* (2013). High-resolution transcriptome maps reveal
  strain-specific regulatory features of multiple *Campylobacter jejuni*
  isolates. *PLoS Genetics* 9(5):e1003495.
  https://doi.org/10.1371/journal.pgen.1003495
* Ettwiller L, Buswell J, Yigit E, Schildkraut I (2016). A novel enrichment
  strategy reveals unprecedented number of novel transcription start sites
  at single base resolution in a model prokaryote and the gut microbiome.
  *BMC Genomics* 17:199. https://doi.org/10.1186/s12864-016-2539-z
* Sharma CM *et al.* (2010). The primary transcriptome of the major human
  pathogen *Helicobacter pylori*. *Nature* 464:250–255.
  https://doi.org/10.1038/nature08756
* Example application of Cappable-seq to an endosymbiont: Suliman Y
  *et al.* (2025). Cappable-Seq reveals the transcriptional landscape of
  stress responses in the bacterial endosymbiont *Wolbachia*.
  *Microbial Genomics* 11(10). https://doi.org/10.1099/mgen.0.001542
