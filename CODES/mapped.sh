#!/bin/bash

mkdir -p /home/stuart/viral/mapped

for f in  /home/stuart/viral/ONT_downloads/*.fastq.gz; do

sample=$(basename $f .fastq.gz)

echo "Mapping and sam $sample"

docker run --rm \
-v /home/stuart/viral:/data \
quay.io/biocontainers/minimap2:2.24--h7132678_1 \
minimap2 -ax map-ont \
/data/reference/MN908947.3.fasta \
/data/ONT_downloads/${sample}.fastq.gz | \

docker run --rm -i \
-v /home/stuart/viral:/data \
quay.io/biocontainers/samtools:1.17--h00cdaf9_0 \
samtools sort -o /data/mapped/${sample}.sorted.bam

done 
