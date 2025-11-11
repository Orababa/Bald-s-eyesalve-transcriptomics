# Balds-eyesalve-transcriptomics
Understanding the mechanism of action of a natural product (Bald's eyesalve) using gene expression (transcriptomics) analysis


# RNA-Seq Differential Expression Analysis Pipeline

This repository contains the workflow and scripts for a comparative transcriptomic analysis of **_Staphylococcus aureus_** and **_Acinetobacter baumannii_** using RNA sequencing (RNA-Seq).

---

## Overview

The pipeline performs the following steps:

1. **Quality control** of raw FASTQ files  
2. **Adapter trimming and quality filtering**  
3. **Read alignment** to reference genomes  
4. **Read counting** at the gene level  
5. **Differential expression analysis**  
6. **Functional enrichment analysis** (COG, GO, KEGG)

---

## Pipeline Steps

### 1. Quality Control

```bash
fastqc *.fastq.gz
```
Raw FASTQ files were inspected using **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** to assess sequence quality, GC content, and potential contaminants.

---

### 2. Adapter Trimming
```bash
trimmomatic PE -phred33 \
    input_R1.fastq.gz input_R2.fastq.gz \
    output_R1_paired.fastq.gz output_R1_unpaired.fastq.gz \
    output_R2_paired.fastq.gz output_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
Illumina adapter sequences were removed using **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)**.

---

### 3. Read Alignment

#### _S. aureus_
```bash
bowtie2 -x NC009641.1 -1 R1_paired.fastq.gz -2 R2_paired.fastq.gz -S saureus.sam
```
Reads aligned to **_S. aureus_ Newman**  
- **Accession**: [NC_009641.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_009641.1)  
- **Assembly**: GCF_000010465.1

#### _A. baumannii_
```bash
bowtie2 -x NZ_CP046654.1 -1 R1_paired.fastq.gz -2 R2_paired.fastq.gz -S abaumannii.sam
```
Reads aligned to **_A. baumannii_ ATCC 19606**  
- **Accession**: [NZ_CP046654.1](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP046654.1)  
- **Assembly**: GCF_009759685.1

Alignment performed using **[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)** [@LangmeadSS].

---

### 4. SAM to BAM Conversion & Sorting
```bash
samtools view -bS input.sam | samtools sort -o sorted.bam -
samtools index sorted.bam
```
Processed using **[SAMtools](https://www.htslib.org/)** [@LiAD].

---

### 5. Gene-Level Read Counting
```bash
featureCounts -p -t CDS -g locus_tag -a annotation.gtf -o counts.txt *.bam
```
Gene counts generated using **[featureCounts](https://subread.sourceforge.net/)** (Subread package) [@LiaoSS].

---

### 6. Differential Expression Analysis
```R
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","treated","control"))
```
Performed in **R** using the **[DESeq2](https://bioconductor.org/packages/DESeq2/)** package [@LoveHA].

---

### 7. Functional Enrichment Analysis
- **COG**, **GO**, and **KEGG** pathway analysis  
- Performed using **[FUNAGE-Pro v3](https://github.com/Cyanolab/FUNAGE-Pro)** [@deJongKJ]

---

## Software & Tools

| Tool           | Version | Reference |
|----------------|---------|---------|
| FastQC         | v0.11.9 | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| Trimmomatic    | v0.39   | http://www.usadellab.org/cms/?page=trimmomatic |
| Bowtie2        | v2.4.5  | @LangmeadSS |
| SAMtools       | v1.15   | @LiAD |
| featureCounts  | v2.0.3  | @LiaoSS |
| DESeq2         | v1.36.0 | @LoveHA |
| FUNAGE-Pro     | v3      | @deJongKJ |

---

## References

```bibtex
@article{LangmeadSS,
  author = {Langmead, Ben and Salzberg, Steven L.},
  title = {Fast gapped-read alignment with Bowtie 2},
  journal = {Nature Methods},
  year = {2012},
  doi = {10.1038/nmeth.1923}
}

@article{LiAD,
  author = {Li, H. and Handsaker, B. and Wysoker, A. and others},
  title = {The Sequence Alignment/Map format and SAMtools},
  journal = {Bioinformatics},
  year = {2009},
  doi = {10.1093/bioinformatics/btp352}
}

@article{LiaoSS,
  author = {Liao, Yang and Smyth, Gordon K. and Shi, Wei},
  title = {featureCounts: an efficient general purpose program for assigning sequence reads to genomic features},
  journal = {Bioinformatics},
  year = {2014},
  doi = {10.1093/bioinformatics/btt656}
}

@article{LoveHA,
  author = {Love, Michael I. and Huber, Wolfgang and Anders, Simon},
  title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
  journal = {Genome Biology},
  year = {2014},
  doi = {10.1186/s13059-014-0550-8}
}

@article{deJongKJ,
  author = {de Jong, K. J. and others},
  title = {FUNAGE-Pro: Functional Annotation of Genomes and Metagenomes},
  journal = {bioRxiv},
  year = {2023},
  doi = {10.1101/2023.xx.xx.xxxxxx}
}
```
---

*Last updated: November 11, 2025, 03:36 PM GMT*
```
```
