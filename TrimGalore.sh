#!/bin/bash

# 指定输入目录和输出目录
input_dir="/home/bio/RNAseq/input"
output_dir="/home/bio/RNAseq/results"

# 遍历输入目录下的所有.fastq文件
for file in $input_dir/*.fastq
do
    echo "Running Trim Galore for file: $file"
    trim_galore --quality 20 --length 50 --output_dir $output_dir $file
done
