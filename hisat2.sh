#!/bin/bash

# 指定Hisat2的索引文件夹路径
index_dir="/home/bio/RNA-Seqqq/index"

# 指定输入FASTQ文件和输出SAM文件
input_fastq="/home/bio/RNA-Seqqq/input/ERR266411_1000_1.fastq"
output_sam="/home/bio/RNA-Seqqq/output/output.sam"

# 运行Hisat2比对
hisat2 -x $index_dir -U $input_fastq -S $output_sam

echo "Hisat2 alignment finished!"
