 #!/bin/bash

mkdir -p /home/stuart/viral/consensus

for f in /home/stuart/viral/mapped/*.trimmed.sorted.bam; do
  sample=$(basename $f .trimmed.sorted.bam)
  echo "=== Building consensus for ${sample} ==="

docker run --rm \
  -v /home/stuart/viral:/data \
  quay.io/biocontainers/ivar:1.4.4--h077b44d_0 \
 bash -c "samtools mpileup -A -d 0 -Q 0 \
/data/mapped/${sample}.trimmed.sorted.bam | \
ivar consensus \
-p /data/consensus/${sample}.consensus \
-t 0.5 -m 10 " 

done 
