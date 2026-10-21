# ANNOgesic-style TSS annotation for Cappable-seq

`cappable_tss.py` reworks the transcription start site (TSS) step of
ANNOgesic so that it works with Cappable-seq data instead of dRNA-seq.

## Why ANNOgesic needs changing for Cappable-seq

ANNOgesic's `tss_ptss` module runs TSSpredator, which compares a
TEX-treated library against an untreated one (dRNA-seq) and looks for
"steps" in read coverage. Cappable-seq enriches 5'-triphosphorylated
(primary) transcripts directly, so the 5' end of each read in the
enriched library marks the TSS at single-base resolution. The natural
signal is the **count of read 5' ends at each position**, compared with
a non-enriched control, not the coverage step.

## What the script does

### `call`: a native Cappable-seq TSS caller

1. Counts read 5' ends per position and strand (BAM, BED or 5'-end wiggle).
2. Normalises each library to reads per million (RPM).
3. Keeps a position in a replicate if the enriched RPM is at least
   `--min-rpm` and `log2((enriched + pc) / (control + pc))` is at least
   `--min-log2fc`. A position must pass in `--min-replicates` replicates.
4. Merges calls on the same strand within `--cluster-distance` nt and
   keeps the strongest one.
5. Classifies TSS against a GFF3 annotation using the categories ANNOgesic
   and TSSpredator use (Sharma et al. 2010):
   * **Primary**: within `--utr-length` (300) nt upstream of a gene and
     the strongest TSS for that gene.
   * **Secondary**: within that window but not the strongest.
   * **Internal**: inside a gene, on the sense strand.
   * **Antisense**: inside a gene or within `--antisense-distance` (100)
     nt of it, on the opposite strand.
   * **Orphan**: none of the above.

   A TSS can have more than one class (for example `Primary&Internal`).
6. Writes `<prefix>.gff` (ANNOgesic-style attributes: `type`,
   `associated_gene`, `UTR_length`, plus `enriched_rpm`, `control_rpm`,
   `log2_enrichment` and `replicates`) and `<prefix>.tsv`.

```bash
python cappable_tss.py call \
    --enriched enr_rep1.bam enr_rep2.bam \
    --control  ctl_rep1.bam ctl_rep2.bam \
    --annotation genome.gff --feature-types CDS \
    --min-rpm 1.5 --min-log2fc 1 --min-replicates 2 \
    --out-prefix results/TSS
```

The thresholds are defaults, not values taken from a paper. Tune them for
your sequencing depth and enrichment efficiency. Looking at a few known
promoters in a genome browser is the quickest check.

Input notes:

* **BAM** needs `pysam`. Secondary, supplementary, unmapped and QC-fail
  reads are skipped. For paired-end data only `--mate` (default 1) is
  counted; that mate's 5' end must be the RNA 5' end, so check your
  library protocol.
* **BED** needs the strand in column 6 (`bedtools bamtobed` output works).
* **Wiggle** is given as `forward.wig,reverse.wig` and must contain 5'-end
  counts. Negative reverse-strand values (READemption/ANNOgesic style) are
  fine.
* `--control` can be left out. Then only the RPM cut-off applies, which is
  much less specific.

### `annogesic-libs`: feed Cappable-seq into the original ANNOgesic

If you want to stay inside the ANNOgesic workflow (for example to keep
using its TSSpredator parameter optimisation or its downstream modules),
this writes forward and reverse wiggle files and prints the
`--tex_notex_libs` strings. The enriched library is labelled `tex` and the
control `notex`:

```bash
python cappable_tss.py annogesic-libs \
    --enriched enr_rep1.bam enr_rep2.bam \
    --control  ctl_rep1.bam ctl_rep2.bam \
    --outdir ANNOgesic/input/wigs/tex_notex
```

By default the wiggle files hold read coverage, which is the signal
TSSpredator was designed for. Use `--signal 5prime` for 5'-end counts.
Check the printed library strings against `annogesic tss_ptss --help` for
your ANNOgesic version.

## Tests

```bash
python -m unittest discover -s tests -v
```

The BAM test is skipped if `pysam` is not installed.

## References

* Yu S-H, Vogel J, Förstner KU (2018). ANNOgesic: a Swiss army knife for
  the RNA-seq based annotation of bacterial/archaeal genomes.
  *GigaScience* 7(9). https://doi.org/10.1093/gigascience/giy096
* Ettwiller L, Buswell J, Yigit E, Schildkraut I (2016). A novel enrichment
  strategy reveals unprecedented number of novel transcription start sites
  at single base resolution in a model prokaryote and the gut microbiome.
  *BMC Genomics* 17:199. https://doi.org/10.1186/s12864-016-2539-z
* Sharma CM *et al.* (2010). The primary transcriptome of the major human
  pathogen *Helicobacter pylori*. *Nature* 464:250–255.
  https://doi.org/10.1038/nature08756
* Suliman Y *et al.* (2025). Cappable-Seq reveals the transcriptional
  landscape of stress responses in the bacterial endosymbiont *Wolbachia*.
  *Microbial Genomics* 11(10). https://doi.org/10.1099/mgen.0.001542
