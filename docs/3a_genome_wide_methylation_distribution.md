---
layout: default
title: Lesson 4 - genome wide methylation distribution analysis
nav_order: 4
parent: 3. Tutorial
description: A comprehensive guide to understanding epigenetics.
published: true
---

{: .important-title }
> Aim
>
> Obtain a graph with the window distribution of methylation values across a chromosome in three contexts `CG`, `CHG` and `CHH`.
> 
> We will use average values in windows in order to understand how we can represent methylation values across a chromosome as depicted in [Figure 1](#figure-1).

({{ "/assets/images/chr5_CTX_report.png" | relative_url }})
> Figure 1: Methylation values across a chromosome. Each point represents the average methylation value in window
<!--
# Analysis of Methylation

In this analysis, we will refer to the methylation distribution shown in [Figure 1](#figure-1-methylation-distribution).

## Figure 1: Methylation Distribution

![Methylation Distribution](path/to/your/figure.png)
-->

The suffix of the file is `.CX_report.txt` and is as follows:

![CTX_example]({{ "/assets/images/chr5_CTX_report.png" | relative_url }})

The file is tab separated and the columns are in the following order:
1. chromosome
2. coordinate
3. strand
4. number of reads with methylation for the C
5. number of reads without methylation
6. methylation context (CG, CHG, CHH)
7. genomic context

# 1. Filtering of the dataset 
Before reading the file in `R` we need to filter the file in order to remove positions without coverage and by selecting the methylation contexts (`CG`) of interest.

awk command 

# 2. Upload the data in `R`


# 3. Creating genomic windows of a fixed size
The chromosome will be divided in windows of a fixed size (e.g. 100,000 bp). The size is totally arbitrary. The methylation values will be averaged in each window. The windows are not overlapping.

Next we need to assing each C position to a genomic window (for example the C in position 11,450 will be assigned to the window 1-100,000, while the C in position 153,000 will be assigned to the window 100,000-200,000).

<!--We will use the `cut` function in `R` to assign each C position to a genomic window-->
We can perform this operation in a vectorized way using the column with the genomic coordinates (from the CG data.frame) and a new vector of assigned intervals (will become a new column in the CG data.frame).

We need to get the start and end of each chromosome. We can do this by using the `cut` function in `R` with the `breaks`. The `breaks` argument is used to specify the points at which to split the data. The 


--- 

We will now repeat the same analysis as done previously, but in this case the windows and raw data to be imported in `R` will be obtained using `bedtools`.

Bedtools is a software package for the manipulation of genomic datasets. It is designed to be fast, flexible, and relatively easy to use. It is particularly useful for the analysis of large genomic datasets. It is written in C++ and is available for Linux, MacOS and Windows. It is free and open source. It is available at [http://bedtools.readthedocs.io/en/latest/](http://bedtools.readthedocs. 

We will use the same raw table obtained from Bismark:

# 1. Filtering of the dataset 
Before reading the file in `R` we need to filter the file in order to remove positions without coverage and by selecting the methylation contexts (`CG`) of interest.

```bash
awk 'OFS="\t" {if($1=="Chr1" && ($4+$5)>0 && $6=="CG") {meth=100*($4)/($4+$5); print $0,meth}}' file > ..
```

# 2. Create the windows of a fixed size using bedtools 
We will use `bedtools makewindows` to create the windows. It requires the size of the window and the chromosome length. We will use the same size of the window as previously.

#### The chromosome size is obtained using `bedtools getfasta` with the `-fo` option to get the length of the chromosome. It is a tab separated file with the following columns:
1. chromosome
2. length


```bash
bedtools makewindows -g chromosome_size.txt -w 100000 > windows.bed
```

In order to assign the single sites to the windows we will use `bedtools intersect` with the `-wa` option. We need to create a file with the single sites in bed format (or bedgraph format).

```bash
awk 'OFS="\t" {print $1,$2-1,$2,$8}' .. > methylome.bed 
```

Now we can **intersect** the Cs sites with the windows:

```bash
bedtools intersect -a methylome.bed -b windows.bed -wa -wb > methylome_windows.bed
```

Now we have the individuals methylation values assigned to the windows. We can obtain a mean methylation value for each window by using `bedtools groupby`:

```bash
bedtools groupby -i methylome_windows.bed -g 1-3 -c 7 -o mean > methylome_windows_mean.bed
```