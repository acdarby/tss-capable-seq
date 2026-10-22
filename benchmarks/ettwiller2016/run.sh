#!/usr/bin/env bash
# Re-analyse the E. coli Cappable-seq data of Ettwiller et al. 2016
# (BMC Genomics 17:199, doi:10.1186/s12864-016-2539-z) and score the calls
# against the authors' TSS table.
#
# usage: run.sh WORKDIR TABLE_S1.gtf
#   TABLE_S1.gtf  the paper's supplementary GTF (12864_2016_2539_MOESM1_ESM.gtf,
#                 16,359 TSS), downloaded from the article page.
#
# needs: bowtie2, samtools, cutadapt, sra-tools, pigz, curl, prodigal, barrnap,
#        aragorn (see environment.yml) and python3 with pysam.
# Only AWS Open Data mirrors are used, so NCBI and ENA access is not needed.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work=$1; table_s1=$2
threads=${THREADS:-4}
mkdir -p "$work" && cd "$work"

# 1. Genome: U00096.2 (= NC_000913.2, 4,639,675 bp) from the AWS iGenomes mirror.
if [ ! -s genome.fa ]; then
  curl -sS https://ngi-igenomes.s3.amazonaws.com/igenomes/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/WholeGenomeFasta/genome.fa \
    | sed '1s/.*/>U00096.2/' > genome.fa
fi
len=$(grep -v '>' genome.fa | tr -d '\n' | wc -c)
[ "$len" -eq 4639675 ] || { echo "unexpected genome length $len"; exit 1; }
# Curated GenBank annotation of U00096.2 (the one the paper used); if NCBI is not
# reachable, fall back to predicted genes.
if [ ! -s annotation.gff ]; then
  curl -sS "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=U00096.2" > annotation.gff || true
  grep -q -P "\tCDS\t" annotation.gff 2> /dev/null \
    || bash "$repo/benchmarks/tools/annotate_genome.sh" genome.fa U00096.2 annotation.gff
fi
[ -e idx.1.bt2 ] || bowtie2-build -q genome.fa idx

# 2. Reads (ENA PRJEB9717): enriched, technical replicate, non-enriched control.
for pair in ERR930221:enriched ERR930222:replicate ERR930223:control; do
  acc=${pair%%:*}; name=${pair##*:}
  [ -s $name.bam ] && continue
  THREADS=$threads bash "$repo/benchmarks/tools/fetch_sra.sh" $acc reads $name
  # As in the paper: cutadapt with the NEBNext small RNA adapter, Bowtie2 local -L 16.
  # -m 1 drops reads that are only adapter; it does not change the mapped reads.
  cutadapt -j $threads -m 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o reads/$name.trim.fastq.gz reads/$name.fastq.gz > $name.cutadapt.log
  bowtie2 --local -L 16 -p $threads -x idx -U reads/$name.trim.fastq.gz 2> $name.bowtie2.log \
    | samtools sort -@ 2 -o $name.bam -
  samtools index $name.bam
  rm -f reads/$name.trim.fastq.gz
done

# 3. TSS calls with the published settings, and comparison with Table S1.
sed 's/^gi|48994873|gb|U00096.2|/U00096.2/' "$table_s1" > table_s1.gtf
"${PYTHON:-python3}" "$repo/cappable_tss.py" call --preset ettwiller2016 \
  --enriched enriched.bam --control control.bam \
  --annotation annotation.gff --feature-types CDS,rRNA,tRNA --out-prefix tss_ettwiller2016
"${PYTHON:-python3}" "$repo/cappable_tss.py" compare --query tss_ettwiller2016.gff --reference table_s1.gtf \
  --windows 0,1,5,25 --out offsets.tsv | tee comparison.txt

# 4. The statistics reported in the paper, for Table S1 and for our calls.
for f in table_s1.gtf tss_ettwiller2016.gff; do
  echo "== $f"; "${PYTHON:-python3}" "$repo/benchmarks/tools/tss_stats.py" $f genome.fa annotation.gff
done | tee statistics.txt
