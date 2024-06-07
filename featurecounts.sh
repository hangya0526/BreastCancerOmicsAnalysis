#!/bin/bash
featureCounts \
-a ../../annotation/* \
-o ../../results/final_counts.txt \
-g 'gene_name' \
-T 4 \
$dirlist
