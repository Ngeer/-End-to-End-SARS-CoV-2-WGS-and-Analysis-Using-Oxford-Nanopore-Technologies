#!/bin/bash
mkdir -p /home/stuart/viral/pangolin_results
mkdir -p /home/stuart/viral/nextclade_results

for f in /home/stuart/viral/consensus/*.fa; do
  sample=$(basename "$f" .fa)
  echo "=== Running Pangolin for $sample ==="
  docker run --rm \
    -v /home/stuart/viral:/data \
    staphb/pangolin:latest \
    pangolin /data/consensus/${sample}.fa \
    --outdir /data/pangolin_results/${sample}

  echo "=== Running Nextclade for $sample ==="
  docker run --rm \
    -v /home/stuart/viral:/data \
    nextstrain/nextclade:latest \
    nextclade run \
    --input-dataset /data/dataset \
    --output-all /data/nextclade_results/${sample} \
    /data/consensus/${sample}.fa
done
