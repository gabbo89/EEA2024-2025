---
layout: default
title: Bismark manual
nav_order: 3
parent: 2. General Guides
description: A brief description of the software used for performing methylation analyses
has_children: true
published: true
---

FINAL VERSION
{: .label .label-green }

![Bismark_figure]({{"/assets/images/bismark.png" | relative_url }})


Bismark is a program to map bisulfite treated sequencing reads to a genome of interest and perform methylation calls in a single step. It's main features are:

- Bisulfite mapping and methylation calling in one single step
- Supports single-end and paired-end read alignments
- Supports ungapped, gapped or spliced alignments
- Alignment seed length, number of mismatches etc. are adjustable
- Output discriminates between cytosine methylation in `CpG`, `CHG` and `CHH` context


Bismark performs alignments of bisulfite treated reads to a reference genome using the short reads aligner Bowtie 2 [link].

## How does it work
Sequence reads are first transformed into fully bisulfite-converted forward (C->T) and reverse read (G->A conversion of the forward strand) versions, before they are aligned to similarly converted versions of the genome (also C->T and G->A converted). Sequence reads that produce a unique best alignment from the four alignment processes against the bisulfite genomes (which are running in parallel) are then compared to the normal genomic sequence and the methylation state of all cytosine positions in the read is inferred. A read is considered to align uniquely if an alignment has a unique best alignment score (as reported by the AS:i field). If a read produces several alignments with the same number of mismatches or with the same alignment score (AS:i field), a read (or a read-pair) is discarded altogether.

![bismark](image-5.png)


## Links

- Bismark Publication
  - http://www.ncbi.nlm.nih.gov/pubmed/21493656
- Our review about primary data analysis in BS-Seq
  - http://www.ncbi.nlm.nih.gov/pubmed/22290186


## Credits

Bismark was written by Felix Krueger, part of the [Babraham Bioinformatics](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) group.

## Licences

Bismark itself is free software, `bismark2report` and `bismark2summary` produce HTML graphs powered by [Plot.ly](https://plot.ly/javascript/) which are also free to use and look at!


[Back to the tutorial](https://gabbo89.github.io/EEA2024/docs/3a1_WGBS_cleaning_and_alignment.html#bismark-meth_extract)