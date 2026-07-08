#!/bin/bash

mkdir -p /home/stuart/viral/variants

for f in /home/stuart/viral/mapped/*.sorted.bam; do
  sample=$(basename $f .sorted.bam)
  echo "=== Processing ${sample} ==="

  # Step 1 - Index BAM
  echo "Indexing ${sample}"
  docker run --rm \
    -v /home/stuart/viral:/data \
    quay.io/biocontainers/samtools:1.17--h00cdaf9_0 \
    samtools index /data/mapped/${sample}.sorted.bam

  # Step 2 - Trim primers
  echo "Trimming primers ${sample}"
  docker run --rm \
    -v /home/stuart/viral:/data \
    quay.io/biocontainers/ivar:1.4.4--h077b44d_0 \
    ivar trim \
    -i /data/mapped/${sample}.sorted.bam \
    -b /data/primer_schemes/nCoV-2019/V3/nCoV-2019.bed \
    -p /data/mapped/${sample}.trimmed \
    -m 30 -q 0 
  # Step 3 - Sort trimmed BAM
  echo "Sorting trimmed ${sample}"
  docker run --rm \
    -v /home/stuart/viral:/data \
    quay.io/biocontainers/samtools:1.17--h00cdaf9_0 \
    samtools sort \
    /data/mapped/${sample}.trimmed.bam \
    -o /data/mapped/${sample}.trimmed.sorted.bam

  # Step 4 - Variant calling
  echo "Calling variants ${sample}"
  docker run --rm \
    -v /home/stuart/viral:/data \
    quay.io/biocontainers/ivar:1.4.4--h077b44d_0 \
    bash -c "samtools mpileup -A -d 0 -Q 0 \
    /data/mapped/${sample}.trimmed.sorted.bam | \
    ivar variants \
    -p /data/variants/${sample}.variants \
    -r /data/reference/MN908947.3.fasta"

done
