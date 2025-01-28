---
layout: default
title: Using bedtools
parent: Lesson 4 - Genome wide methylation distribution analysis
nav_order: 2
description: A comprehensive guide to understanding epigenetics.
published: true
---
INCOMPLETE
{: .label .label-red }

{: .important-title }
> Aim
>
> Obtain a graph with the window distribution of methylation values across a chromosome in three contexts `CG`, `CHG` and `CHH`, using `bedtools`.

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


We will now repeat the same analysis as done previously (a chromosome distribution of the methylation), but in this case the windows and raw data to be imported in `R` for the final graph will be obtained using `bedtools`.

Bedtools is a software package for the manipulation of genomic datasets. It is designed to be fast, flexible, and relatively easy to use. It is particularly useful for the analysis of large genomic datasets. It is written in C++ and is available for Linux, MacOS and Windows. It is free and open source. It is available at [http://bedtools.readthedocs.io/en/latest/](http://bedtools.readthedocs.io/en/latest/).

The file has been already used in the previous tutorial. (add link)

The file is located at the following path:

`/data2/biotecnologie_molecolari_magris/epigenomics/meth_distribution/arabidopsis_wgbs.CX_report.txt`

It should be already available in your directory:
`/data2/student_space/st24_16_folder/epigenomics/methylation_distribution/`

# 1. Filter the dataset 
We need to filter the file in order to remove positions without coverage and by selecting the methylation contexts (`CG`) for the chromosome of interest.

```bash
# Move the working directori
cd /data2/student_space/st24_16_folder/epigenomics/

# Create a new directory for this tutorial
mkdir -p genome_wide_meth/

# Filter the input file in order to keep only the methylation context of interest (CG) and to keep sites located on Chr1 with a coverage greater than 0
awk -v "OFS=\t" '{if($1=="Chr1" && ($4+$5)>0 && $6=="CG") {meth=100*($4/($4+$5)); print $0,meth}}' methylation_distribution/arabidopsis_wgbs.CX_report.txt > genome_wide_meth/arabidopsis_chr1_CG_meth.txt
```

![bedtools_table]({{"/assets/images/image-15.png" | relative_url }})

# 2. Create windows of fixed size
We will use `bedtools makewindows` to create the windows. It requires the size of the **window** and the **chromosome length**. We will use the same size of the window as previously.

The **chromosome length** is obtained using `samtools faidx`. The output is a tab separated file with the following columns:
1. **chromosome name**
2. **sequence length**
3. **offset** _# byte offset of the chromosome in the FASTA file_{: .label .label-green-100 }
4. **line bases**
5. **line width** _# number of bytes in each line_ {: .label .label-green-100 }


Thus using the first two columns, we are able to get the chromsome size. But we will need the fasta file of the reference genome!

{: .highlight-title}
> Question
>
> How do we get the reference genome?
>

<details>
    <summary>Show answer</summary>
We can look for the fasta sequencing by performing a search on google for example.<br>
<br>
Try to type in google: Arabidopsis thaliana genome fasta<br>
<br>

da rimuovere eventualmente!!!!
![google search Arabidopsis]({{"/assets/images/image-16.png" | relative_url }})

<img src="{{ '/assets/images/image-16.png' | relative_url }}" alt="google search Arabidopsis">
<br>
You will find different options. Try to navigate and look for the fasta file on ncbi database. <br>

</details>

<br>
<!-- Hidden link -->
<div id="hidden-link" style="display:none;">
  <a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.3/" target="_blank">https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.3/</a>
</div>

<!-- Button to toggle the visibility of the link -->
<button onclick="document.getElementById('hidden-link').style.display='block'; this.style.display='none';">Show Fasta link</button>

You can find an example PowerPoint presentation on OneDrive [here](https://onedrive.live.com/?cid=YOUR_CID&resid=YOUR_RESID&authkey=YOUR_AUTHKEY&action=embedview).

We can use the following command to download the fasta file:

```bash
# Create a directory to store the reference genome 
mkdir -p genome_wide_meth/reference

# Download the reference genome using wget 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz -P genome_wide_meth/reference/

# Unzip the fasta file 
gunzip genome_wide_meth/reference/GCF_000001735.4_TAIR10.1_genomic.fna.gz

# Create the index file
samtools faidx genome_wide_meth/reference/GCF_000001735.4_TAIR10.1_genomic.fna
```

The reference genome is also available at [http://plants.ensembl.org/Arabidopsis_thaliana/Info/Index](http://plants.ensembl.org/Arabidopsis_thaliana/Info/Index)



We will set the **windows size** to 100,000bp.

In order to generate the windows, we will use the following command:

```bash
bedtools makewindows \
-g <(cut -f 1,2 genome_wide_meth/reference/GCF_000001735.4_TAIR10.1_genomic.fna | grep "Chr1") -w 100000 > genome_wide_meth/100k_windows.bed
```
# 3. Intersect the two tables 
In order to assign the single sites to the windows we will use `bedtools intersect` with the `-wa` option. We need to create a file with the single sites in bed format (or bedgraph format). We need to format the columns of the input file by creating a bed-like structure:

```bash
awk 'OFS="\t" {print $1,$2-1,$2,$8}' genome_wide_meth/arabidopsis_chr1_CG_meth.txt > genome_wide_meth/methylome.bed 
```

Now we can **intersect** the Cs sites with the previously created windows:

```bash
bedtools intersect `
-a genome_wide_meth/methylome.bed \
-b genome_wide_meth/100k_windows.bed \
-wa -wb > genome_wide_meth/methylome_windows.bed
```

Now we have the individuals methylation values assigned to the windows. 
We can obtain a mean methylation value for each window by using `bedtools groupby`:

```bash
bedtools groupby \
-i genome_wide_meth/methylome_windows.bed \
-g 1-3 \
-c 7 \
-o mean > genome_wide_meth/methylome_windows_mean.bed
```
