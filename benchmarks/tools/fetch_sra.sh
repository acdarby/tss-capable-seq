#!/usr/bin/env bash
# Download one SRA run from the AWS Open Data mirror and keep only a compressed
# FASTQ. Paired runs are written interleaved (R1, R2, R1, R2 ...); use
# `bowtie2 --interleaved` or split them later.
#
# usage: fetch_sra.sh ACCESSION OUTDIR [LABEL]
# needs: curl, fasterq-dump (sra-tools), pigz
set -euo pipefail
acc=$1; outdir=$2; label=${3:-$1}
mkdir -p "$outdir"
out="$outdir/$label.fastq.gz"
if [ -s "$out" ]; then echo "exists: $out"; exit 0; fi
curl -sS --retry 4 -o "$outdir/$acc.sra" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$acc/$acc"
# --size-check off: fasterq-dump otherwise refuses when its temp-space estimate
# exceeds free disk; fall back to streaming fastq-dump if it still fails.
fasterq-dump --size-check off -t "$outdir" --split-spot --skip-technical -Z "$outdir/$acc.sra" 2> "$outdir/$label.dump.log" \
    | pigz -p "${THREADS:-4}" > "$out" || true
if [ "$(stat -c %s "$out")" -lt 1000 ]; then
    fastq-dump --split-spot --skip-technical --stdout "$outdir/$acc.sra" 2>> "$outdir/$label.dump.log" \
        | pigz -p "${THREADS:-4}" > "$out"
fi
rm -f "$outdir/$acc.sra"
echo "$label: $(grep -i 'reads written' "$outdir/$label.dump.log")"
