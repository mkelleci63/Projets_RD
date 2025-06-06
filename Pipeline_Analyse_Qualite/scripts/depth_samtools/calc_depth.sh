#!/bin/bash
set -euo pipefail

sample="$1"
region_file="$2"
filtered_bam="$3"
output_dir="$4"
mkdir -p "$output_dir"

while IFS=$'\t' read -r chrom start end name rest; do
  region="${chrom}:${start}-${end}"
  outfile="${output_dir}/${chrom}_${start}_${end}_${name}.txt"
  
  samtools mpileup -B -q 1 -Q 0 -F 1024 -r "${region}" "$filtered_bam" | \
    awk '{
      gsub(/\^./, "", $4)
      gsub(/\$/, "", $4)
      print $1 "\t" $2 "\t" length($4)
    }' > "$outfile"
done < "$region_file"

