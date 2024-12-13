---
layout: default
title: Lesson 1 - Raw data cleaning and alignment
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

- `--path_to_cutadapt cutadapt` trimming software loading
- `--phred33`  Phred quality scores (DEFAULT)
- `--illumina` for Illumina adapters
- `--paired` for paired sequences
- `--trim1` to elude the software that discards overlapping reads
- `--clip_R1 20` cut 20 nt in R1 (5’)
- `--clip_R2 6` cut 6 nt in R2 (5’)
- `--three_prime_clip_R1 4`  cut 4 additional nt in 3’ because adapters clipping leads to residual low-quality bases
- `--three_prime_clip_R2 4`  cut 4 additional nt in 3’ because adapters clipping leads to residual low-quality bases

### Perform a second round of quality control on the trimmed data

### Run fastQC on the trimmed files
```bash
fastqc ..R1.fastq.gz ...R2.fastq.gz
```

Open the obtained figures from the output folder in order to evaluate the quality of the processed data. 

---

# 3. Alignment of fastq files 
In order to perform the alignment we will use the Bismark suite [Bismark on Github][bismark_github]<sup>[1]</sup>.

<!--
In order to perform the alignment we will use the Bismark suite [^Bismark short manual] [Bismark short manual][bismark short manual] and [^TrimGalore on Github][trimgalore_github].
-->

{: .warning }
> Be sure that the reference genome has the required indexes

### Create the indexes required by Bismark (only once)
```bash
bismark_genome_preparation --path_to_bowtie bowtie2_folder --verbose genome_folder
```
1. bowtie2_folder is the folder containing the bowtie2 software (the command requires the folder name rather the executable).
The folder is:
2. genome_folder is the folder containing the reference fasta file. The preparation command should create and additional folder, inside the `genome_folder`, called `Bisulfite_Genome` > Bisulfite_Genome


### Perform the paired-end mapping 
USAGE: bismark [options] <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}

```bash
bismark --bowtie2 --bam --phred33-quals -N 1 -p 2 genome_folder -1 [file R1.fq.gz pathway] -2 [file R2.fq.gz pathway]
```

#### the options 
- `--bowtie2 bowtie2` is used as the backend (DEFAULT).
- `--bam` alignment is written in bam format (DEFAULT).
- ../genome genome directory (not entered as a parameter per se, but rather directly in the line).
- `--phred33-quals` Quality format: ASCII chars equal to the Phred quality plus 33 (valid for current Illumina data) (DEFAULT).
- `-N 1` Sets the number of mismatches to be allowed in a seed alignment during multiseed alignment (a bowtie2 property that allows for higher sensitivity).
- `-p 2` Number of cores used for bowtie alignment.
- `-1` read1 file
- `-2` read2 file

The ouput of the aligment process is a `bam file` containing mapping results that can be read using `samtools`.
`Samtools` is a suite of commands that can be used for manipulating sam/bam files. 

Also a .txt file is stored with the report of mapping efficiency that can be read with a normal textual reader command:
For example:
```bash
less -S ...
```

{: .note-title }
>The most important values are:
>
>> Sequences analyzed in total
>
>> Mapping efficiency 
>
>> Number of alignments with a unique best hit from the different alignments 
>
>> Sequences did not map uniquely 


It is possible to calculate:
<!--
$$
Total efficiency(%) = (Number of alignments with a unique best hit from the different alignments + sequences did not map uniquely) / Sequences in input
$$

$$ 
\text{Total efficiency(\%)} = \frac{(\text{Number of alignments with a unique best hit from the different alignments} + \text{sequences did not map uniquely})}{\text{Sequences in input}} 
$$

\[
\text{Total efficiency(\%)} = \frac{(\text{Number of alignments with a unique best hit from the different alignments} + \text{sequences did not map uniquely})}{\text{Sequences in input}} \text{ (1)}
\]
-->
> Total efficiency(%) = (Number of alignments with a unique best hit from the different alignments + sequences did not map uniquely) / Sequences in input
<!--
> % methylation (context) = 100 \* methylated Cs (context) / (methylated Cs (context) + unmethylated Cs (context)).
-->
---

# 4. Deduplication and methylome extraction
We need to remove duplicated reads from the alignment file that may have originated from PCR errors.3
- add a comment to why duplicated reads need to be removed

### Perform deduplicate 
```bash
deduplicate_bismark --bam rkatsiteli.leaves.rkatsiteli.leaves.R1_bismark_bt2_pe.bam
```
In order to extract methylation information we will run the command `bismark_methylation_extractor`. The script will operate on Bismark result files and extracts the methylation call for every single C analysed. 

The position of every single C will be written to a new output file, dependending on the context (CG, CHG or CHH), whereby methylated Cs will be labelled as forward read (+) and non-methylated Cs as reverse reads (-). 

### Extract methylation information 
```bash
bismark_methylation_extractor --genome_folder /projects/novabreed/share/gmagris/collaboration/lezioni/2024/EEA/reference/ -p --bedGraph --cytosine_report --CX_context --multicore 1 --gzip rkfatsiteli.leaves.rkatsiteli.leaves.R1_bismark_bt2_pe.deduplicated.bam
```

>#### the options 
> `-p` for processing paired-end data
> `--bedGraph` 
> `--genome-folder`
> `--cytosine_report` 
> `--CX_context` 
> `--multicore` 
> `--gzip` 


{: .success-title }
> STDOUT
>
> Writing genome-wide cytosine report to: rkatsiteli.leaves.rkatsiteli.leaves.R1_bismark_bt2_pe.deduplicated.CX_report.txt.gz
>
> Writing all cytosine context summary file to: rkatsiteli.leaves.rkatsiteli.leaves.R1_bismark_bt2_pe.deduplicated.cytosine_context_summary.txt
>
> Finished writing out cytosine report for covered chromosomes (processed 343 chromosomes/scaffolds in total)
>
> Now processing chromosomes that were not covered by any methylation calls in the coverage file...
>
> Writing cytosine report for not covered chromosome h1tg000079l
> Writing cytosine report for not covered chromosome h1tg000179l
> Finished writing out cytosine report (processed 345 chromosomes/scaffolds in total). coverage2cytosine processing complete.
>
> Finished generating genome-wide cytosine report!


Several files will be produced in this last step. The most important file is the `CX_report.txt` that contain the methylome data.

The files obtained have the following structure.
For example the file `CpG_OT_rkatsiteli.leaves.rkatsiteli.leaves.R1_bismark_bt2_pe.deduplicated.txt.gz` has the following columns:
1. seq-ID
2. methylation state
3. chromosome
4. start position
5. methylation call

Methylated Cs have `+` orientation, while unmethylated Cs have `-`

For example: 

>SEQILMN03:348:CAG91ANXX:8:1101:13304:1989_1:N:0:TGGTGA  +       chr17   13398887        Z
>SEQILMN03:348:CAG91ANXX:8:1101:13304:1989_1:N:0:TGGTGA  +       chr17   13398905        Z
>SEQILMN03:348:CAG91ANXX:8:1101:13304:1989_1:N:0:TGGTGA  +       chr17   13398949        Z
>SEQILMN03:348:CAG91ANXX:8:1101:14588:1992_1:N:0:TGGTGA  +       chr18   34541292        Z
>SEQILMN03:348:CAG91ANXX:8:1101:14588:1992_1:N:0:TGGTGA  +       chr18   34541341        Z
>SEQILMN03:348:CAG91ANXX:8:1101:14588:1992_1:N:0:TGGTGA  +       chr18   34541361        Z
>SEQILMN03:348:CAG91ANXX:8:1101:14588:1992_1:N:0:TGGTGA  +       chr18   34541344        Z
>SEQILMN03:348:CAG91ANXX:8:1101:19965:1990_1:N:0:TGGTGA  -       chr18   35977426        z

It will produce a strand-specific output which will use the following abbreviations in the output file name that indicate the strand the alignment came from:

> OT    - original `TOP` strand
> CTOT  - complementary to original `TOP` strand
> OB    - original `BOTTOM` strand
> CTOB  - complementary to original `BOTTOM` strand

{: .note }
Methylation calls from OT and CTOT will be informative for cytosine methylation positions on the original top strand, calls from OB and CTOB will be informative for cytosine methylation positions on the original bottom strand. Please note that specifying the --directional (the default mode) option in the Bismark alignment step will not report any alignments to the CTOT or CTOB strands.

## BedGraph output
The Bismark methylation extractor can optionally also output a file in bedGraph format which uses 0-based genomic start and 1- based end coordinates. 
The columns are as follows:
1. `chromosome`
2. ``start position``
3. `end position`
4. value (methylation %)

Since the methylation percentage is _per se_ not informative of the read coverage at the specific posiition, a `.cov.gz` file is also created (1-based genomic coordinates) that feature 2 additional columns:
1. `chromosome`
2. `start position`
3. end position
4. value (methylation %)
5. number of methylated Cs
6. number of unmethylated Cs
[link to descriptor](/docs/2a_file_formats.md)
From this file, downstream processing of the file. 

{: .note}
Only performed on CG sites

## M-bias output 
<!--
This allows generating nice graphs by alternative means, e.g. using R or Excel. The plot is also drawn into a .png file which requires the Perl module GD::Graph (more specifically, both modules GD::Graph::lines and GD::Graph::colour are required); if GD::Graph cannot be found on the system, only the table will be printed.
-->
The Bismark methylation extractor can optionally also output a file in M-bias format which uses 0-based genomic start and 1- based end coordinates.?
Methylation bias plot which shows the methylation proportion across each possibile position in the read (cumulatively)[^2]

The output is a tabular file with the following format:
1. `read position`
2. `count methylated`
3. `count unmethylated`
4. `% methylation`
5. `total coverage`



# splitting_report
It represent a summary of the splitting step executed by bismark_methylation_extractor. It the report the % of methylated Cs in the different contexts. 


---

# 5. Manipulating the bam file 
The filtered bam file obtained after deduplicate_bismark, is still unsorted for coordinates.

### sort the bam file 
```bash
samtools sort \
-@ 2 \
-o rkfatsiteli.leaves.bismark_bt2_pe.deduplicated.sort.bam \
rkfatsiteli.leaves.rkatsiteli.leaves.R1_bismark_bt2_pe.deduplicated.bam
```
### index the bam file 
```bash
samtools index -@ 2 rkfatsiteli.leaves.bismark_bt2_pe.deduplicated.sort.bam
```

---

# 6. Conversion rate evaluation 
In order to understand if the conversion rate of the cytosine worked, we need to verify the bisulfite conversion rate. 
We can use the chloroplast genome (or lambda genome)

### Index the chloroplast fasta
```bash
bismark_genome_preparation `
--path_to_aligner /iga/scripts/dev_modules/mambaforge/envs/epigenomics/bin/ \
--bowtie2 \
--parallel 10 \
--verbose \
/projects/novabreed/share/gmagris/collaboration/lezioni/2024/EEA/chloroplast/
```

{: .success-title }
> STDOUT
>
>'========================================='
>
>Parallel genome indexing complete. Enjoy!

Now we are ready to perform the aligment of the reads to the chloroplast genome

### Perform the paired-end mapping 
```bash
bismark \
--bowtie2 \
--bam \
--phred33-quals \
-N 1 \
-p 2 \
genome_folder \
-1 [file R1.fq.gz pathway] \
-2 [file R2.fq.gz pathway]
```

Results are reported in *bismark_bt2_PE_report.txt file!

[trimgalore short manual]: https://gabbo89.github.io/EEA2024/docs/2a_TrimGalore_manual.html
[trimgalore_github]: https://github.com/FelixKrueger/TrimGalore
[^1]: https://gabbo89.github.io/EEA2024/docs/2a_Bismark_manual.html
<sup>[1]</sup> [Bismark short manual](https://gabbo89.github.io/EEA2024/docs/2a_Bismark_manual.html)
[bismark_github](https://felixkrueger.github.io/Bismark/){: .btn }
[^2]: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r83
