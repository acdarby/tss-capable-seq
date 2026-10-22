#!/usr/bin/env bash
# TSSpredator and ANNOgesic modules on the wAlbB Cappable-seq data of
# Frontiers in Microbiology 2026 (doi:10.3389/fmicb.2026.1860204, PRJNA1392518):
# stranded, random-primed RNA-seq of captured primary transcripts (read 2 = sense).
#
# usage: run.sh WORKDIR CAPPSEQ.bam NOENRICH.bam SULIMAN_TSS.gff
#   BAMs: paired reads mapped to wAlbB (CP031221.1); SULIMAN_TSS.gff: TSS called
#   from 5'-end Cappable-seq (e.g. `cappable_tss.py call` on PRJNA1145481).
# needs: docker (silasysh/annogesic:1.1.14 for TSSpredator), ANNOgesic 1.1.14
#        from PyPI, TransTermHP and ViennaRNA (bioconda), curl, python3 + pysam.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd); repo=$(cd "$here/../.." && pwd)
work=$1; enr=$2; ctl=$3; tss=$4
mkdir -p "$work/in" && cd "$work"
acc=GCA_004171285.1
if [ ! -s in/wAlbB.fa ]; then
  curl -sS -o pkg.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA,GENOME_GFF"
  unzip -o -q pkg.zip
  cp ncbi_dataset/data/$acc/*.fna in/wAlbB.fa
  grep -v "^#" ncbi_dataset/data/$acc/genomic.gff | awk -F'\t' 'NF==9' > in/wAlbB.gff
fi
grep "^CP031221.1" "$tss" | sed '1i ##gff-version 3' > in/wAlbB_TSS.gff

# Strand-specific coverage of read 2 (the sense read), chromosome only.
python3 "$repo/cappable_tss.py" annogesic-libs --enriched "$enr" --control "$ctl" \
  --mate 2 --min-mapq 20 --outdir wigs_all > /dev/null
mkdir -p wigs
for f in wigs_all/*.wig; do
  awk '/^variableStep/{keep=($2=="chrom=CP031221.1")} /^track/||keep' "$f" > wigs/$(basename "$f")
done
tf=$(ls wigs/*_tex_forward.wig); tr=$(ls wigs/*_tex_reverse.wig)
nf=$(ls wigs/*_notex_forward.wig); nr=$(ls wigs/*_notex_reverse.wig)
W=$(pwd)

# 1. TSSpredator (enriched as TEX+, unenriched as TEX-), default parameters.
docker run --rm -v "$W:$W" -w "$W" silasysh/annogesic:1.1.14 bash -c "
  annogesic create --project_path $W/proj_tss > /dev/null &&
  annogesic tss_ps --project_path $W/proj_tss --program TSS --fasta_files $W/in/wAlbB.fa \
    --annotation_files $W/in/wAlbB.gff --replicate_tex all_1 --condition_names cappable \
    --tex_notex_libs $W/$tf:tex:1:a:+ $W/$tr:tex:1:a:- $W/$nf:notex:1:a:+ $W/$nr:notex:1:a:-"
python3 "$repo/cappable_tss.py" compare --query proj_tss/output/TSSs/gffs/wAlbB_TSS.gff \
  --reference in/wAlbB_TSS.gff --windows 0,5,25,50 | tee tsspredator_vs_5prime.txt

# 2. ANNOgesic transcript, terminator, UTR and operon modules (both libraries as
#    fragment libraries; TSS from the 5'-end data).
P=$W/proj; annogesic create --project_path $P > /dev/null
frag="$W/$tf:frag:1:a:+ $W/$tr:frag:1:a:- $W/$nf:frag:2:a:+ $W/$nr:frag:2:a:-"
annogesic transcript --project_path $P --annotation_files in/wAlbB.gff --frag_libs $frag --replicate_frag all_1
annogesic terminator --project_path $P --fasta_files in/wAlbB.fa --annotation_files in/wAlbB.gff \
  --transcript_files $P/output/transcripts/gffs/wAlbB_transcript.gff --frag_libs $frag --replicate_frag all_1 \
  --transterm_path "$(command -v transterm)" --rnafold_path "$(command -v RNAfold)" \
  --expterm_path "$(dirname "$(command -v transterm)")/../data/expterm.dat"
annogesic utr --project_path $P --annotation_files in/wAlbB.gff --tss_files in/wAlbB_TSS.gff \
  --transcript_files $P/output/transcripts/gffs/wAlbB_transcript.gff \
  --terminator_files $P/output/terminators/gffs/best_candidates/wAlbB_term.gff
# Transcript 5' ends from random-primed coverage lie ~12 nt downstream of the TSS,
# so the default 5 nt tolerance misses most TSS-operon links.
annogesic operon --project_path $P --annotation_files in/wAlbB.gff --tss_files in/wAlbB_TSS.gff \
  --transcript_files $P/output/transcripts/gffs/wAlbB_transcript.gff \
  --terminator_files $P/output/terminators/gffs/best_candidates/wAlbB_term.gff --tss_tolerance 25
