#!/bin/bash
sortmerna \
--ref $sortmernaREF \
--reads results/sample_trimmed.fq \
--aligned results/sample_aligned.fq \
--other results/sample_filtered.fq \
--fastx \
--log \
-a 4 \
-v
