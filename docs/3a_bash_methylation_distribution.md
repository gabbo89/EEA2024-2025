---
layout: default
title: Lesson 3 - methylation distribution analysis
nav_order: 3
parent: 3. Tutorial
description: A comprehensive guide to understanding epigenetics.
published: false
---

{: .important-title }
> Aim
>
> Obtain a graph with the distribution of methylation values in three contexts `CG`, `CHG` and `CHH`.

We will us the methylation table obtained froom Bismark.

The file is located at the following path

`/data2/biotecnologie_molecolari_gmagris/bismark_methylome/...CX_report.txt`

The suffix of the file is `.CX_report.txt` and is as follows:

> insert image with the methylome of the extracted region 

The file is tab separated and the columns are in the following order:
1 chromosome
2 coordinate
3 strand
4 number of reads with methylation for the C
5 number of reads without methylation 

![CTX_example]({{ "/assets/images/chr5_CTX_report.png" | relative_url }})

---

The first step is to filter the data in order to select:
- positions covered by at least 1 read
- positions that belong to a specific context (CG, CHG, CHH)
- calculate the methylation level for each position (and append to the table)

We can opt to use `R`, but since the file is quite large, we will use `linux` commands (in particular `awk`).

The input file is located at `/data2/biotecnologie_molecolari_gmagris/bismark_methyl`

### Copy the table of interest to your working directory
```bash
cp /data2/biotecnologie_molecolari_gmagris/bismark_methylome/... /data2/student_space/st24_01_folder/epigenomics/methylation distribution
```
If I want to select only C in CG context and with a coverage greater than 0, i can use `awk` in order to filter the table and select data of interest:

```bash
awk '{ if (($4+$5)>0 && $6=="CG") {meth = $4/($4+$5); print $0"\t"meth}}' file_bismark > Arabidopisis_metiloma_CG.txt
```

The same `awk` command may be used to extract values for CHG and CHH contexts:
```bash
awk '{ if (($4+$5)>0 && $6=="CHG") {meth = $4/($4+$5); print $0"\t"meth}}' file_bismark > Arabidopisis_metiloma_CHG.txt
```

and 

```bash
awk '{ if (($4+$5)>0 && $6=="CHH") {meth = $4/($4+$5); print $0"\t"meth}}' file_bismark > Arabidopisis_metiloma_CHH.txt
```

Now that we obtained a filtered dataset (and lighter ), we can proceed with the analysis of the methylation distribution. We can use `R` for this purpose.
### Load the libraries and the data 
R 

upload the tydiverse library 
```r
library(tidiverse)
```

Now we can load the data and check the structure of the table.
We will store the data in a data.frame called CG, using the `read.table` function. We will also specify the path to the file and the separator used in the file (tab).
```r
CG=read.table(“Arabidopisis_metiloma_CG.txt ", stringsAsFactors=F, header=F, sep="\t")
```

