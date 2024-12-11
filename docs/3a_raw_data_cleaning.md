---
layout: default
title: Lesson 1 - raw data cleaning
nav_order: 1
parent: 3. Tutorial
description: A comprehensive guide to understanding epigenetics.
---

# Brief description of the lesson


## Table of contents

- [Quality control](# 1. First step check fastq quality using fastqc software )
    - [FastQC](#fasta)
- [Trimming](# 2. Second step: perform trimming of raw data )
    - [TrimGalore](#Trimgalore)

---

# Quality control of raw fastq files
Sequencing raw data are divided in read1 and read2 and store in `fastq` format. Fastq files are compressed in `gzip` format (*.fastq.gz)

We need to check the quality of the raw data in order to be sure that sequencing worked!

We will use FastQC software to verify if raw data quality is appropriate and thus planning the quality trimming.
Fastqc is availbale both as graphical and textual interface (we will use the textual).  

---

## 1. First step check fastq quality using fastqc software 

### Activate the conda environment
```bash
conda activate epigenomics
```

#### Test if fastqc is working
```bash
fastqc --help
```
# add succes 

### Copy the raw data from the folder to our working directory 
```bash
cp ... /data2/student_space/st24_01_folder/...
```

### Run Fastqc 

```bash
fastqc ..R1.fastq.gz ...R2.fastq.gz
```

#### Check the output files obtained 
Refer to the presentation for additional informations.  


----

## 2. Second step: perform trimming of raw data 

Once the quality is evaluated we can procede by removing low quality bases from the fastq files.
Native reads will be subject to quailty and adapter trimming before the alignment. Clipping of additional bases at 5' and/or 3' end may deemed necessary in certain circumstances.

### Peform trimming of raw data
We will use TrimGalore to remove adapter and low quality data from fastq file [TrimGalore short manual][trimgalore short manual] and [TrimGalore on Github][trimgalore_github]



```bash
trim_galore --path_to_cutadapt cutadapt --phred33 --illumina --paired --trim1 [file R1.fastq.gz pathway] [file R2.fastq.gz pathway]
```
```bash
trim_galore --path_to_cutadapt cutadapt --phred33 --illumina --paired --trim1 --clip_R1 20 --clip_R2 6 --three_prime_clip_R1 4 --three_prime_clip_R2 4 [file R1.fastq.gz pathway] [file R2.fastq.gz pathway]
```

The most important options are:

- –path_to_cutadapt cutadapt  trimming software loading
- --phred33  Phred quality scores (DEFAULT)
- --illumina  for Illumina adapters
- --paired  for paired sequences
- --trim1  to elude the software that discards overlapping reads
- --clip_R1 20  cut 20 nt in R1 (5’)
- --clip_R2 6  cut 6 nt in R2 (5’)
- --three_prime_clip_R1 4 #  cut 4 additional nt in 3’ because adapters clipping leads to residual low-quality bases
- --three_prime_clip_R2 4 #  cut 4 additional nt in 3’ because adapters clipping leads to residual low-quality bases

### Perform a second round of quality control on the trimmed data

### Run fastQC on the trimmed files
```bash
fastqc ..R1.fastq.gz ...R2.fastq.gz
```

Open the obtained figures from the output folder in order to evaluate the quality of the processed data. 

---

# 3. Alignment of fastq files 
In order to perform the alignment we will use the Bismark suite [Bismark short manual][bismark short manual] and [TrimGalore on Github][trimgalore_github].

# First we need to 
## #### bismark alignment 

> bismark alignment

`bismark` 
#### perform deduplicate 
```bash
deduplicate_bismark --bam rkatsiteli.leaves.rkatsiteli.leaves.R1_bismark_bt2_pe.bam
```

### extract methylation information 
```bash
bismark_methylation_extractor --genome_folder /projects/novabreed/share/gmagris/collaboration/lezioni/2024/EEA/reference/ -p --bedGraph --cytosine_report --CX_context --multicore 1 --gzip rkatsiteli.leaves.rkatsiteli.leaves.R1_bismark_bt2_pe.deduplicated.bam
```


{: .success-title }
> STDOUT of command
> Finished writing out cytosine report for covered chromosomes (processed 343 chromosomes/scaffolds in total)
>
> Now processing chromosomes that were not covered by any methylation calls in the coverage file...
>
> Writing cytosine report for not covered chromosome h1tg000079l
> Writing cytosine report for not covered chromosome h1tg000179l
> Finished writing out cytosine report (processed 345 chromosomes/scaffolds in total). coverage2cytosine processing complete.
>
> Finished generating genome-wide cytosine report!'



[trimgalore short manual]: https://gabbo89.github.io/EEA2024/docs/2a_trim_galore_manual.html
[trimgalore_github]: https://github.com/FelixKrueger/TrimGalore
[bismark short manual]: https://github.com/FelixKrueger/TrimGalore
[bismark_github]: https://github.com/FelixKrueger/TrimGalore