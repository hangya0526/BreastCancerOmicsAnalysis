# RNA-seq基本工作流程

#### 作者：袁媛、程愉航

## 摘要：

RNAseq正在成为测量细胞反应的最突出的方法之一。RNAseq不仅能够分析样本之间基因表达的差异，而且可以发现新的亚型并分析SNP变异。本流程将介绍处理和分析差异基因表达数据的基本工作流程，旨在提供设置环境和运行比对工具的通用方法。

## 前言：

##### 项目介绍：

RNA-Seq是一种高通量测序技术，可用于研究生物体内的RNA表达情况。相比传统的基因芯片技术，RNA-Seq具有更高的灵敏度和更广泛的动态范围，因此被广泛应用于分析基因表达、转录本结构、RNA编辑和差异剪接等生物学过程。

##### 工作内容：

1.用FastQC分析序列质量 

2.使用Trim_Galore删除低质量序列 

3.用SortMeRNA去除rRNA序列

4.使用 STAR-aligner 对齐基因组

5.使用 featureCounts 汇总基因计数

6.使用 multiQC 生成分析报告

## 数据集与方法：

#### 1.配置环境

```bash
#配置环境
conda create -n rnaseq python=3.8
conda activate rnaseq
mkdir RNAseq
conda install -c bioconda fastqc
conda install -c bioconda trim-galore
conda install -c bioconda sortmerna=2.1b
conda install -c bioconda star
conda install -c bioconda subread
conda install -c bioconda multiqc

#下载数据
wget -P input/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/001/SRR1374921/SRR1374921.fastq.gz
wget -P input/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/002/SRR1374922/SRR1374922.fastq.gz
wget -P input/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/003/SRR1374923/SRR1374923.fastq.gz
wget -P input/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/004/SRR1374924/SRR1374924.fastq.gz

```

#### 2.用FastQC分析序列质量

工具介绍：

FastQC旨在提供一种简单的方法，对来自高通量测序管道的原始序列数据进行一些质量控制检查。它提供了一组模块化的分析，可以使用它来快速了解您的数据是否存在任何问题。在处理任何样本之前，第一步是分析数据的质量。文件中包含质量信息，这些信息涉及每个基本调用的准确性（置信度百分比），查看样品序列的不同方面，以确定影响结果的任何不规则或特征。

命令：

```shell
bash fastqc.sh
```

输出：

```
sample_fastqc.html
sample_fastqc.zip
```

#### 3.用Trim_Galore去除低质量序列

工具介绍：

Trim_galore用于自动执行质量和接头修剪以及质量控制，并增加了一些功能以消除RRBS序列文件的偏置甲基化位置。结合了 Cutadapt和FastQC，在执行质量分析以查看过滤效果的同时删除低质量序列。

命令：

```shell
bash TrimGalore.sh 
```

输出：

```shell
sample_trimmed.fq                 
sample.fastq.trimming_report.txt 
```

#### 4.用SortMeRNA去除rRNA序列

工具介绍：

SortMeRNA 是一种程序工具，用于在宏转录组学和宏基因组学数据中过滤、定位和 OTU 拾取 NGS 读数。核心算法基于近似种子，可以对核苷酸序列进行快速、灵敏的分析。SortMeRNA的主要应用是从宏转录组数据中过滤核糖体RNA。一旦我们去除了低质量序列并去除了任何接头污染，我们就可以从样品中去除rRNA序列。

生成索引：

```sh
# 建立索引库
wget -P sortmerna_db https://github.com/biocore/sortmerna/archive/2.1b.zip
unzip sortmerna_db/2.1b.zip -d sortmerna_db
mv sortmerna_db/sortmerna-2.1b/rRNA_databases/ sortmerna_db/
rm sortmerna_db/2.1b.zip
rm -r sortmerna_db/sortmerna-2.1b

# 将所有数据库的位置保存到一个文件夹sortmerna_db中
sortmernaREF=sortmerna_db/rRNA_databases/silva-arc-16s-id95.fasta,sortmerna_db/index/silva-arc-16s-id95:\
sortmerna_db/rRNA_databases/silva-arc-23s-id98.fasta,sortmerna_db/index/silva-arc-23s-id98:\
sortmerna_db/rRNA_databases/silva-bac-16s-id90.fasta,sortmerna_db/index/silva-bac-16s-id95:\
sortmerna_db/rRNA_databases/silva-bac-23s-id98.fasta,sortmerna_db/index/silva-bac-23s-id98:\
sortmerna_db/rRNA_databases/silva-euk-18s-id95.fasta,sortmerna_db/index/silva-euk-18s-id95:\
sortmerna_db/rRNA_databases/silva-euk-28s-id98.fasta,sortmerna_db/index/silva-euk-28s-id98

indexdb_rna --ref $sortmernaREF
```

命令：

```bash
bash sortmerna.sh  
```

输出：

```
sample_aligned.fq
sampleLog.final.out                
sampleLog.out
```

#### 5.使用 STAR-aligner 对齐基因组

（本步骤也可以用hisat2、tophat等）

工具介绍：

STAR基于以前未描述的RNA-seq比对算法的剪接转录本比对，该算法在未压缩的后缀数组中使用顺序最大可映射种子搜索，然后进行种子聚类和拼接程序。STAR在映射速度方面比其他对齐器高出>50倍，在适度的12核服务器上每小时对齐人类基因组5.5亿次2×76 bp的双端读长，同时提高了对准灵敏度和精度。STAR 是一种非常快速和高效的剪接对准器工具，用于将 RNAseq 数据与基因组进行比对，具有发现非经典剪接和嵌合转录本的能力，但对于我们的用例，我们将用于将全长 RNA 序列与基因组进行比对。

（1）生成索引

    STAR \
    --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles genome/* \
    --sjdbGTFfile annotation/* \
    --runThreadN 4

（2）运行命令

    bash star.sh

（3）处理输出

```bash
mv -v results/4_aligned_sequences/sampleAligned.sortedByCoord.out.bam results/4_aligned_sequences/aligned_bam/

mv -v results/4_aligned_sequences/${BN}Log.final.out results/4_aligned_sequences/aligned_logs/
mv -v results/4_aligned_sequences/sample*Log.out results/4_aligned_sequences/aligned_logs/
```

输出：

``` bash
sampleAligned.sortedByCoord.out.bam
sampleLog.final.out                 
sampleLog.out                        
```

------------------------------------------------------------------------

#### 6.使用 featureCounts 汇总基因计数

工具介绍：

featureCounts 是一个高效的通用读取摘要程序，可以计算基因组特征（如基因、外显子、启动子、基因体、基因组箱和染色体位置）的映射读取。它可用于计数 RNA-seq 和基因组 DNA-seq 读数。featureCounts 将 SAM/BAM 文件和包含特征染色体坐标的注释文件作为输入。它输出分配给特征的读取次数。它还输出总体汇总结果的统计信息，包括成功分配的读取次数和由于各种原因而未能分配的读取次数。

命令：

    # Store list of files as a variable
    dirlist=$(ls -t ./*.bam | tr '\n' ' ')
    echo $dirlist
    
    bash featurecounts.sh

输出：

``` bash
final_counts.txt                
final_counts.txt.summary       
```

------------------------------------------------------------------------

#### 7.使用 multiQC 生成分析报告

工具介绍：

MultiQC是一种创建单一报告的工具，将多个工具在多个样本中的输出可视化，从而能够快速识别全球趋势和偏见。MultiQC可以绘制来自许多常见生物信息学工具的数据，并且易于扩展和定制。在质量过滤、rRNA 去除、STAR比对和基因总结过程中，创建了多个日志文件，其中包含衡量相应步骤质量的指标。我们可以使用汇总工具MultiQC，而不是迭代许多不同的日志文件，它将搜索所有相关文件并生成丰富的图形，显示来自不同步骤日志文件的数据。在确定序列与基因组的对齐程度并提取每个步骤中丢失的序列数量时，此步骤非常有用。

命令：

    # Run multiqc and output results into final folder
    multiqc results \
    --outdir results

输出 :

``` bash
multiqc_report.html
multiqc_data
```

## 结果：

可以成功生成FastQC的质量分析报告

![image-20240606145015812](C:\Users\糖糖\AppData\Roaming\Typora\typora-user-images\image-20240606145015812.png)可以得到样本删除低质量序列后的文件以及其质量分析报告、Trim_Galore删减的报告、sortmeRNA去除的报告（不过由于文件采用二进制格式，因此无法查看内容）

得到final_counts报告

![image-20240606151043789](C:\Users\糖糖\AppData\Roaming\Typora\typora-user-images\image-20240606151043789.png)

## 讨论：

项目的关键：

（1）完整运行了RNA-seq的基本工作流程

（2）利用了FastQC、Trim_Galore、SortMeRNA等工具对于目标基因进行处理

不足与改进：

（1）对于目标基因样本的处理方式较少，产生的分析报告不够全面

（2）虽然项目中包括了FastQC和Trim_Galore进行质量控制和序列修剪，但可能还需要进一步的质控和数据清洗步骤，例如去除PCR重复序列和低复杂度序列。

（3）由于时间与技术限制，项目中没有包括差异表达分析步骤，这是RNA-Seq数据分析的重要环节之一，后续可用R继续进行差异性分析。

（4）最后生成的分析报告可能需要进一步解释和验证，以确保分析结果的可靠性和可解释性。

## 贡献：

由袁媛和程愉航共同完成

## 参考文献：

1. Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. 
2. Martin, Marcel. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, [S.l.], v. 17, n. 1, p. pp. 10-12, may. 2011. ISSN 2226-6089. 
3. Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611
4. Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21. doi:10.1093/bioinformatics/bts635.
5. Lassmann et al. (2010) "SAMStat: monitoring biases in next generation sequencing data." Bioinformatics doi:10.1093/bioinformatics/btq614 [PMID: 21088025]

## 附录：

所有内容同步发布在github上：[hangya0526/RNA-SeqWorkflow (github.com)](https://github.com/hangya0526/RNA-SeqWorkflow)

项目核心脚本

(1) fastqc.sh

```bash
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
```

(2) TrimGalore.sh

```bash
#!/bin/bash

# 指定输入目录和输出目录
input_dir="/home/bio/RNA
    echo "Running Trim Galore for file: $file"
    trim_galore --quality 20 --length 50 --output_dir $output_dir $file
done
```

(3)sortmerna.sh

```bash
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
```

(4)star.sh

```bash
#!/bin/bash
# Run STAR (~10min)
STAR \
--genomeDir star_index \
--readFilesIn filtered/sample_filtered.fq  \
--runThreadN 4 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts
```

(5)featurecounts.sh

```bash
#!/bin/bash
featureCounts \
-a ../../annotation/* \
-o ../../results/final_counts.txt \
-g 'gene_name' \
-T 4 \
$dirlist
```

