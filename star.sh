#!/bin/bash
# Run STAR (~10min)
STAR \
--genomeDir star_index \
--readFilesIn filtered/sample_filtered.fq  \
--runThreadN 4 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts
