---
layout: default
title: Lezione 1 - raw data cleaning
nav_order: 1
parent: 3. Tutorial
description: A comprehensive guide to understanding epigenetics.
---


### Check quality of fastq files 
----
# 1. first step check fastq quality using fastqc software 
## Fastqc is availbale both as graphical and textual interface (we will use the textual)

#### upload the module required 
`conda activate epigenomics`

#### test if fastqc is working
`fastqc --help`

#### copy the raw data from the folder to our working directory 
`cp ... /data2/student_space/st24_01_folder/...`

#### run Fastqc 

`fastqc ..R1.fastq.gz ...R2.fastq.gz`

----

### Peform trimming of raw data 
Usie TrimGalore to remove adapter and low quality data from fastq file [TrimGalore manual][def] and [TrimGalore on Github][trimgalore_github]

**dasda**

[def]: https://gabbo89.github.io/EEA2024/docs/3a_trim_galore_manual.html
[def2]: https://github.com/FelixKrueger/TrimGalore