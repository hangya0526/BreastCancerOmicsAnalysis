#!/bin/bash

# 指定.fastq文件所在的目录
input_dir="/home/bio/RNAseq/input"

# 创建输出目录
output_dir="/home/bio/RNAseq/results"

# 运行FastQC分析每个.fastq文件
for file in $input_dir/*.fastq
do
    echo "Running FastQC for file: $file"
    fastqc $file -o $output_dir
done
