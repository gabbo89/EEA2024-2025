---
layout: default
title: Lesson 4 - Genome wide methylation distribution analysis
nav_order: 4
parent: 3. Tutorial
description: A comprehensive guide to understanding epigenetics.
has_children: true
published: true
---

{: .important-title }
> Aim
>
> Obtain a graph with the window distribution of methylation values across a chromosome in three contexts `CG`, `CHG` and `CHH`.


<br>
<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
- TOC
{:toc}
</details>
<br>

We will use average values in windows in order to understand how we can represent methylation values across a chromosome as depicted in [Figure 1](#figure-1).


![Figure 1: Methylation values across a chromosome. Each point represents the average methylation value in window]({{ "/assets/images/chr5_CTX_report.png" | relative_url }})

<!--
# Analysis of Methylation

In this analysis, we will refer to the methylation distribution shown in [Figure 1](#figure-1-methylation-distribution).

## Figure 1: Methylation Distribution

![Methylation Distribution](path/to/your/figure.png)
-->

The suffix of the file is `.CX_report.txt` and is as follows:



The file is tab separated and the columns are in the following order:
1. chromosome
2. coordinate
3. strand
4. number of reads with methylation for the C
5. number of reads without methylation
6. methylation context (CG, CHG, CHH)
7. genomic context

<!--
Now we will perform the analysis of the methylation distribution in the three contexts `CG`, `CHG` and `CHH` across a chromosome. 
-->
We will use two different approaches: one using the `R` library and the other using the `bedtools` package.