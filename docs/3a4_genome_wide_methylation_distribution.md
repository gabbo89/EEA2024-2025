---
layout: default
title: Lesson 4 - Genome wide methylation distribution analysis
nav_order: 4
parent: 3. Tutorial
description: A comprehensive guide to understanding epigenetics.
has_children: true
published: true
---

INCOMPLETE
{: .label .label-red }

{: .important-title }
> Aim
>
> Obtain a graph with the window distribution of methylation values across a chromosome in three contexts `CG`, `CHG` and `CHH`.


<br>
<details open markdown="block">
  <summary>
    <strong>Table of contents</strong>
  </summary>
  {: .text-delta }
- TOC
{:toc}
</details>
<br>



We will use the methylation table obtained from Bismark. The file represent the result of wgbs performed in _Arabidopsis thaliana_ sample.

The file is located at the following path:

`/data2/biotecnologie_molecolari_magris/epigenomics/meth_distribution/arabidopsis_wgbs.CX_report.txt`

The suffix of the file is `.CX_report.txt` as already seen in the previous lessons.
The structure is as follows:

![Figure 1: header of the CG data frame]({{ "/assets/images/3a3-0_methylation_distribution_arabidopsis.png" | relative_url }})
<br>

**Figure 1:** First rows of the *arabidopsis_wgbs.CX_report.txt* file.

The file is tab separated and the columns are in the following order:
1. **chromosome**
2. **coordinate**
3. **strand**
4. **number of reads with methylation for the C**
5. **number of reads without methylation**
6. **C-context**
7. **trinucleotide context**

We will use average values in windows in order to understand how we can represent methylation values across a chromosome as depicted in [Figure 1](#figure-1).

Aggiungere FIGURA genome wide 

![Figure 1: chromosome_wide methylation]({{ "/assets/images/chr5_CTX_report.png" | relative_url }})

<!--
# Analysis of Methylation

In this analysis, we will refer to the methylation distribution shown in [Figure 1](#figure-1-methylation-distribution).

## Figure 1: Methylation Distribution

![Methylation Distribution](path/to/your/figure.png)
-->


<!--
Now we will perform the analysis of the methylation distribution in the three contexts `CG`, `CHG` and `CHH` across a chromosome. 
-->
We will use two different approaches: one using the `R` library and the other using the `bedtools` package.