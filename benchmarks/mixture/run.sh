#!/usr/bin/env bash
# Mixed-genome benchmark: E. coli (Ettwiller et al. 2016) and Salmonella SL1344
# (PRJEB58986) Cappable-seq reads, mapped together to both genomes and mixed
# with E. coli at 100%, 10% and 2% of its reads, against each organism's reads
# mapped to its own genome only. See README.md in this folder.
#
# usage: run.sh DATA_DIR OUT_DIR   (needs bowtie2, samtools, python3 with pysam)
# DATA_DIR is laid out as in ../datasets.toml (data/ecoli/*.bam,
# data/ecoli/genome.fa, data/salmonella/*.bam, genomes/...).
set -euo pipefail
data=$(cd "$1" && pwd); out=$(mkdir -p "$2" && cd "$2" && pwd)
repo=$(cd "$(dirname "$0")/../.." && pwd)
py=${PYTHON:-python3}; threads=${THREADS:-4}
ecoli_fa=$data/data/ecoli/genome.fa
sal_fa=$data/genomes/SL1344.fa
ecoli_gff=$data/genomes/U00096.2.gff3
sal_gff=$data/genomes/GCA_000210855.2/ncbi_dataset/data/GCA_000210855.2/genomic.gff

cd "$out"
cat "$ecoli_fa" "$sal_fa" > combined.fa
[ -f combined.1.bt2 ] || bowtie2-build -q --threads "$threads" combined.fa combined

# 1. Each library mapped once to both genomes (same settings as before: local, -L 16).
map() {  # in.bam out.bam
  [ -f "$2" ] && return
  samtools fastq -F 0x904 "$1" 2>/dev/null \
    | bowtie2 --local -L 16 -p "$threads" -x combined -U - 2> "${2%.bam}.bowtie2.log" \
    | samtools view -b -F 4 - | samtools sort -@ 2 -o "$2" -
  samtools index "$2"
}
map "$data/data/ecoli/ERR930221.bam" ecoli_a.combined.bam
map "$data/data/ecoli/ERR930222.bam" ecoli_b.combined.bam
map "$data/data/salmonella/SL1344_cappable_rep1.bam" sal_1.combined.bam
map "$data/data/salmonella/SL1344_cappable_rep2.bam" sal_2.combined.bam

# 2-4. Mixtures and the single-genome references.
for frac in 1 0.1 0.02; do
  tag=ecoli$frac
  sub() {  # in.bam out.bam: the same read names in every file (hash of the name)
    if [ "$frac" = 1 ]; then samtools view -b -F 4 "$1" -o "$2"
    else samtools view -b -F 4 -s "7${frac#0}" "$1" -o "$2"; fi
    samtools index "$2"
  }
  sub ecoli_a.combined.bam $tag.a.combined.bam
  sub ecoli_b.combined.bam $tag.b.combined.bam
  sub "$data/data/ecoli/ERR930221.bam" $tag.a.alone.bam
  sub "$data/data/ecoli/ERR930222.bam" $tag.b.alone.bam
  samtools merge -f -@ 2 mix_$tag.1.bam sal_1.combined.bam $tag.a.combined.bam
  samtools merge -f -@ 2 mix_$tag.2.bam sal_2.combined.bam $tag.b.combined.bam
  samtools index mix_$tag.1.bam; samtools index mix_$tag.2.bam
  for mapq in 20 0; do
    "$py" "$repo/cappable_tss.py" call --enriched mix_$tag.1.bam mix_$tag.2.bam \
      --genome Ecoli="$ecoli_fa" Salmonella="$sal_fa" \
      --annotation "$ecoli_gff" "$sal_gff" --feature-types CDS \
      --min-mapq $mapq --out-prefix calls/mix_${tag}_q$mapq 2> calls_mix_${tag}_q$mapq.log
  done
  "$py" "$repo/cappable_tss.py" call --enriched $tag.a.alone.bam $tag.b.alone.bam \
    --genome Ecoli="$ecoli_fa" --annotation "$ecoli_gff" --feature-types CDS \
    --min-mapq 20 --out-prefix calls/alone_$tag 2> calls_alone_$tag.log
  rm -f mix_$tag.*.bam* $tag.*.combined.bam*
done
"$py" "$repo/cappable_tss.py" call \
  --enriched "$data/data/salmonella/SL1344_cappable_rep1.bam" "$data/data/salmonella/SL1344_cappable_rep2.bam" \
  --genome Salmonella="$sal_fa" --annotation "$sal_gff" --feature-types CDS \
  --min-mapq 20 --out-prefix calls/alone_salmonella 2> calls_alone_salmonella.log
"$py" "$repo/benchmarks/mixture/compare.py" "$out/calls" > "$out/mixture.tsv"
echo "wrote $out/mixture.tsv"
