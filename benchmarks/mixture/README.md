# Mixed-genome benchmark

Tests `call` with several genomes (one FASTA per organism) on a sample that
contains two related bacteria, and in particular a rare one next to an
abundant one, as with a secondary symbiont next to a primary one.

* Reads: *Salmonella* SL1344 Cappable-seq, two replicates (PRJEB58986), and
  *E. coli* MG1655 Cappable-seq enriched library and its second library
  (ERR930221, ERR930222; Ettwiller et al. 2016, doi:10.1186/s12864-016-2539-z)
  as two replicates.
* Each library is mapped once to both genomes together (bowtie2 local,
  `-L 16`, as for the single-genome data). The *Salmonella* BAMs used as the
  read source held only reads that had mapped to *Salmonella*.
* Mixed samples: each *Salmonella* replicate plus *E. coli* reads at 100%,
  10% or 2% (`samtools view -s`, which keeps the same read names in every
  file). Standard mode with `--genome Ecoli=... Salmonella=...` at
  `--min-mapq 20` and `0`.
* Reference: the same reads of each organism mapped to its own genome only
  and called alone (`--min-mapq 20`).

```bash
PATH=/path/to/bowtie2-and-samtools:$PATH bash run.sh DATA_DIR OUT_DIR
```

Results: `../results/mixture_ecoli_salmonella.tsv` (columns explained in
`compare.py`) and `../results/mixture_mapping.txt`; summary in
`../README.md` ("Several bacteria in one sample").
