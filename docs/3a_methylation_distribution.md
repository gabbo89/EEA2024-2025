---
layout: default
title: Lesson 3 - methylation distribution analysis
nav_order: 3
parent: 3. Tutorial
description: A comprehensive guide to understanding epigenetics.
published: true
---

{: .important-title }
> Aim
>
> Obtain a graph with the distribution of methylation values in three contexts `CG`, `CHG` and `CHH`.

We will us the methylation table obtained froom Bismark.

The file is located at the following path

`/data2/biotecnologie_molecolari_gmagris/bismark_methylome/...CX_report.txt`

The suffix of the file is `.CX_report.txt` and the structure is as follows:

![CTX_example]({{ "/assets/images/chr5_CTX_report.png" | relative_url }})


The file is tab separated and the columns are in the following order:
1. chromosome
2. coordinate
3. strand
4. number of reads with methylation for the C
5. number of reads without methylation 

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

### Filter the data and calculate the methylation level
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
```R
library(tidiverse)
```

Now we can load the data and check the structure of the table.
We will store the data in a data.frame called CG, using the `read.table` function. We will also specify the path to the file and the separator used in the file (tab).

```r
CG=read.table(“Arabidopisis_metiloma_CG.txt ", stringsAsFactors=F, header=F, sep="\t")
```

# rename the columns 
```r
names(CG)=c('chr', 'pos', 'strand', 'c', 't', 'context', 'real_context', 'methylation')
```

# add a new named `coverage` which include the total coverage

```r
CG$coverage = CG$c + CG$t
```

Add now a new column named `methR` which represent the methylation level calculated in a different way in the R environment. The value is calculated as % value, rounded to 0 decimal places:

```r
CG$methR = round(100*CG$c / CG$coverage, 0)
```
Now we can filter the table by removing the rows where the coverage is lower than a certain threshold (e.g. 10). We haven't done it previously with AWK in order to test the different coverage thresholds in R. Removing the non covered positions (done with [awk](#filter-the-data-and-calculate-the-methylation-level)) can be done at the beginning because they are not informative. 

We will use the `dplyr` library (from the `tidyverse` package) to filter the data.


```r
CG_coverage_filtered = CG %>% filter(coverage > 10)
```

If, as commonly happens, the number of Cs with methylation values = 0 is extremely high, the graph may appear compressed and hard to understand on the __X__ axis. Thus it might be useful to remove the rows where the methylation is 0. This can be done with the following command:

```r
CG_coverage_filtered = CG %>% filter(coverage > 10 & methR > 0)
```

#### Repeat now the same analysis for CHG and CHH contexts.
## CHG 


## CHH


# Draw the plot in R
We need to upload the necessary libraries . We will use `ggplot2` for the plot and `dplyr` for the filtering. We will also use `tidyverse` for the data manipulation. We will load the libraries with the following command:

```r
library(ggplot2)
library(dplyr)
```

Draw the graph as histogram:

```r
ggplot(CG_coverage_filtered,aes(x=methR)) + \
geom_histogram(colour=4,fill="white",binwidth=1)
```

Draw the graph as density plot:

```r
ggplot(CG_coverage_filtered,aes(x=methR)) + \
geom_density(alpha=.2,fill="#FF6666")
```

#### Repeat now the same for CHG and CHH.

We can save the graph in a pdf file with the following command:

```r
pdf("CG_density.pdf",paper="A4")
ggplot(CG_coverage_filtered,aes(x=methR)) + \
geom_density(alpha=.2,fill="#FF6666")
dev.off()
```

Or we can save all the plots in a single file:

```r
pdf("CG_density.pdf",paper="A4")
ggplot(CG_coverage_filtered,aes(x=methR)) + \
geom_density(alpha=.2,fill="#FF6666")

ggplot(CHG_coverage_filtered,aes(x=methR)) + \
geom_density(alpha=.2,fill="#FF6666")

ggplot(CHH_coverage_filtered,aes(x=methR)) + \
geom_density(alpha=.2,fill="#FF6666")

dev.off()
```

We need to pair `pdf()` and `dev.off()` in order to open (draw the graph) and close the pdf file. We can also use `png()` or `jpeg()` to save the plots in a png or jpeg file. With `png()` or `jpeg()` we need an extra step to arrange the figures in multiple panels.