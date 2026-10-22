# TSSpredator and ANNOgesic on wAlbB fragment-library Cappable-seq

Data: *Wolbachia* wAlbB in *Aedes albopictus* Aa23 cells, Frontiers in
Microbiology 2026 ([doi:10.3389/fmicb.2026.1860204](https://doi.org/10.3389/fmicb.2026.1860204)),
SRA PRJNA1392518. Whole primary transcripts were captured, then
heat-fragmented and sequenced with NEBNext Ultra II Directional (random-primed,
stranded, read 2 = sense, 2 × 150 bp). The libraries used here are the
Cappable-seq and unenriched libraries of the XRN-1 experiment
(SRR41042961, SRR41042962), mapped with Bowtie2 `--local` to wAlbB
(GCA_004171285.1) plus the *Ae. albopictus* mitochondrion.

These are not 5'-end libraries, so `cappable_tss.py call` does not apply:
read starts are spread over fragments, and coverage rises at a TSS and
continues along the transcript. `run.sh` instead runs TSSpredator (via the
ANNOgesic Docker image) and the ANNOgesic transcript, terminator, UTR and
operon modules, using TSS called from the 5'-end data of Suliman et al. 2025
(same strain and cell line) for the UTR and operon steps.

```bash
bash run.sh WORKDIR cappseq.bam noenrich.bam suliman_tss.gff
```

## Results (October 2026)

**TSSpredator** (enriched as TEX+, unenriched as TEX−, default parameters)
called 250 TSS (73% primary). Compared with the 1,527 TSS called from the
Suliman 5'-end data (MAPQ ≥ 20, ≥ 10 CPM, 10 nt clusters, all three replicates):

| Window | TSSpredator TSS near a 5'-end TSS | 5'-end TSS recovered |
|---|---|---|
| exact | 1.2% | 0.2% |
| ±5 nt | 2.4% | 0.4% |
| ±25 nt | 89.2% | 14.4% |
| ±50 nt | 90.0% | 15.6% |

TSSpredator finds the right promoters but places them a median 12 nt
downstream of the 5'-end TSS: random-primed libraries lose the first bases of
each RNA, so the coverage rise starts after the true TSS. With default
parameters it recovers only 14% of the 5'-end TSS.

The same data give an independent check of 5'-end TSS: 73.1% of the
Suliman-settings TSS show a ≥ 1.5-fold coverage rise (50 nt downstream vs
upstream) in the Cappable-seq library, against 6.7% of positions shifted by
100 nt.

**ANNOgesic modules** (both libraries as fragment libraries):

| Module | Result |
|---|---|
| transcript | 2,085 transcripts (median 204 nt, up to 19.6 kb); 775 of 1,395 CDS lie fully inside one |
| terminator | 1,039 TransTermHP candidates, 597 in expressed regions, 302 with a coverage drop ("best") |
| UTR | 638 5' UTRs (median 28 nt, IQR 18–75) and 450 3' UTRs |
| operon | 608 operons, 228 (37.5%) polycistronic; 36% start at a TSS with `--tss_tolerance 25` (12.5% with the default 5) |

ANNOgesic accepted the TSS GFF written by `cappable_tss.py` in the UTR and
operon modules. Because transcript 5' ends from random-primed coverage lie
about 12 nt downstream of TSS, the operon module needs a wider TSS tolerance
than its dRNA-seq default.
