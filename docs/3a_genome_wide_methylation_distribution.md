---
layout: default
title: Lesson 3 - genome wide methylation distribution analysis
nav_order: 3
parent: 3. Tutorial
description: A comprehensive guide to understanding epigenetics.
published: true
---

{: .important-title }
> Aim
>
> Obtain a graph with the window distribution of methylation values across a chromosome in three contexts `CG`, `CHG` and `CHH`.
> 
> We will use average values in windows in order to undersatnd the how we can represent methylation values across a chromosome as depicted in Figure 1.

> Figure 1: Methylation values across a chromosome. Each point represents the average methylation value in


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