#!/bin/bash


mkdir -p /home/stuart/viral/qc_results

for f in /home/stuart/viral/ONT_downloads/*.fastq.gz; do
  sample=$(basename $f .fastq.gz)
  echo "Processing $sample"
  docker run --rm \
    -v /home/stuart/viral:/data \
    quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0 \
    NanoStat --fastq /data/ONT_downloads/${sample}.fastq.gz \
    --outdir /data/qc_results \
    --name ${sample}
done


