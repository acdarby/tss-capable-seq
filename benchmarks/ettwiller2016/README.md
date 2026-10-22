# Benchmark: Ettwiller et al. 2016 (*E. coli* Cappable-seq)

Re-analysis of the *E. coli* K-12 MG1655 Cappable-seq data of Ettwiller et
al. 2016 (*BMC Genomics* 17:199,
[doi:10.1186/s12864-016-2539-z](https://doi.org/10.1186/s12864-016-2539-z)),
scored against the authors' own TSS table.

## Run it

```bash
conda env create -f ../environment.yml && conda activate tss-benchmarks
bash run.sh /path/to/workdir /path/to/12864_2016_2539_MOESM1_ESM.gtf
```

The GTF is the paper's supplementary TSS table (16,359 TSS on U00096.2),
available from the article page. Everything else is downloaded from AWS Open
Data mirrors; only the curated annotation comes from NCBI:

| Input | Source |
|---|---|
| Reads ERR930221 (enriched), ERR930222 (technical replicate), ERR930223 (non-enriched control), ENA PRJEB9717 | `sra-pub-run-odp` S3 bucket |
| Genome U00096.2 (= NC_000913.2, 4,639,675 bp) | AWS iGenomes, *E. coli* K-12 MG1655 NCBI 2001-10-15 |
| Gene annotation | GenBank U00096.2 from NCBI; predicted with Prodigal, barrnap and ARAGORN (`../tools/annotate_genome.sh`) if NCBI is not reachable |

Mapping follows the paper: cutadapt with the NEBNext small RNA adapter and
Bowtie2 `--local -L 16`. Outputs are `comparison.txt`, `offsets.tsv`,
`statistics.txt` and the TSS calls `tss_ettwiller2016.gff/.tsv`.

## Results (October 2026)

Agreement of `cappable_tss.py call --preset ettwiller2016` with Table S1:

| Window | Our calls in Table S1 | Table S1 recovered |
|---|---|---|
| exact | 94.3% | 95.9% |
| ±5 nt | 95.1% | 96.8% |
| ±25 nt | 97.0% | 98.2% |

Most of the difference comes from read counts: newer cutadapt and Bowtie2
versions give counts a few percent different from the paper's, which moves
positions across the 1.5 RPM cut-off. On simulated data the preset gives
exactly the same TSS as Ettwiller's Perl scripts.

`statistics.txt` reproduces the paper's figures on Table S1 within about one
percentage point with the GenBank annotation (for example 80.9% A/G at +1,
63.5% of TSS inside a CDS, and exactly the 82 leaderless TSS reported),
except for C/T at −1 (71.4% against 76% reported). With predicted genes the
figures differ by less than one point.

Other findings from this dataset, which motivate the planned changes to `call`:

* The non-enriched control is a 5'-end library, not standard RNA-seq: 76% of
  its reads start in rRNA genes.
* Normalising by library size therefore places the TSS enrichment mode near
  log2 +3, so the published enrichment cut-off of 0 removes little.
* TSSpredator (via ANNOgesic 1.1.14, default settings) finds 1,134 weak,
  strongly enriched sites below the 1.5 RPM cut-off; 29.5% of them have a
  TATAAT-like −10 box, against 33.6% of TSS found by all methods and 2.3% of
  random positions.

`selftrain.py` is a prototype that learns a promoter sequence model from each
genome's own enriched and depleted positions (no external labels):

```bash
python selftrain.py enriched.bam replicate.bam control.bam genome.fa table_s1.gtf weights.json
```
