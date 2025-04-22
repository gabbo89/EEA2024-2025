---
layout: default
title: Lesson 3 - Methylation distribution analysis
nav_order: 3
parent: 3. Tutorial
description: A comprehensive guide to understanding epigenetics.
published: true
---

Final version
{: .label .label-green }

{: .important-title }
> Aim
>
> Obtain a graph with the distribution of methylation values in three contexts `CG`, `CHG` and `CHH`.


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

---

# 1. Filter the Bismark output file.
The first step is to filter the data in order to select:
- positions covered by at least 1 read
- positions that belong to a specific context (CG, CHG, CHH)
- calculate the methylation level for each position (and append to the table)

We can opt to use `R`, but since the file is quite large, we will use `linux` commands (in particular `awk`).

### Activate the conda environment
{: .no_toc }

```bash
conda activate epigenomics
```

### Set the working directory
{: .no_toc }

```bash
# move to the working directory
cd /data2/student_space/st24_16_folder

# create the folder structure
mkdir -p epigenomics/methylation_distribution

# move to the new working directory
cd epigenomics/methylation_distribution
```

### Copy the table of interest to your working directory
{: .no_toc }

```bash
cp /data2/biotecnologie_molecolari_magris/epigenomics/meth_distribution/arabidopsis_wgbs.CX_report.txt .

# in order to check the file format, read the first rows of the tabular file
head arabidopsis_wgbs.CX_report.txt
```


```
Chr3    101     +       0       0       CHH     CCC
Chr3    102     +       0       0       CHH     CCT
Chr3    103     +       0       0       CHH     CTA
Chr3    108     +       0       0       CHH     CCC
Chr3    109     +       0       0       CHH     CCT
Chr3    110     +       0       0       CHH     CTA
Chr3    115     +       0       0       CHH     CCC
Chr3    116     +       0       0       CHH     CCT
Chr3    117     +       0       0       CHH     CTA
Chr3    122     +       0       0       CHH     CCC
```


### Filter the data and calculate the methylation level
{: .no_toc }
If we want to select only C in CG context and with a coverage greater than 0, we can use `awk` in order to filter the table and select data of interest:

```bash
awk '{ if (($4+$5)>0 && $6=="CG") {meth = $4/($4+$5); print $0"\t"meth}}' \
arabidopsis_wgbs.CX_report.txt > arabidopsis_metilome_CG.txt
```

The same `awk` command may be used to extract values for CHG context:

```bash
awk '{ if (($4+$5)>0 && $6=="CHG") {meth = $4/($4+$5); print $0"\t"meth}}' \
arabidopsis_wgbs.CX_report.txt > arabidopsis_metilome_CHG.txt
```

and for CHH context:

```bash
awk '{ if (($4+$5)>0 && $6=="CHH") {meth = $4/($4+$5); print $0"\t"meth}}' \
arabidopsis_wgbs.CX_report.txt > arabidopsis_metilome_CHH.txt
```

# 2. Upload the methylation table in `R`.
Now that we obtained a **filtered** dataset (that impact less memory), we can proceed with the analysis of the methylation distribution. We will use `R` for this purpose.

### Load the R environment
{: .no_toc }

```bash
R 
```

```r
# Load the dplyr library 
library(dplyr)
```

Now we can load the data and check the structure of the table from the command line.
We will store the table in a data.frame called CG, using the `read.table` function. We will also specify the path to the file and the separator used in the file (tab).


```r
# read the input file, which is missing the header
CG=read.table("arabidopsis_metilome_CG.txt", stringsAsFactors=F, header=F,sep="\t")

# we can now check the data.frame using for example the head() function
head(CG)
```
![Figure 2: header of the CG data frame]({{ "/assets/images/3a3-1_methylation_distribution_arabidopsis.png" | relative_url }})
<br>
**Figure 2:** First rows of the CG data frame.


### Rename the columns 
{: .no_toc }

```r
names(CG)=c('chr', 'pos', 'strand', 'c', 't', 'context', 'genome_context', 'methylation')
```

### Add a new column named `coverage` which include the total coverage
{: .no_toc }

```r
# total coverage is calculated by summing the columns c and t
CG$coverage = CG$c + CG$t
```

Add now a new column named `methR` which represent the methylation level calculated in a different way. The value is calculated as % value, rounded to 0 decimal places:

```r
CG$methR = round(100*CG$c / CG$coverage, 0)

# check now the new dataframe
head(CG)
```

![Figure 3: header of the modified CG data frame]({{ "/assets/images/3a3-2_methylation_distribution_arabidopsis.png" | relative_url }})
<br>
**Figure 3:** First rows of the modified CG data frame.


Now we can filter the table by removing the rows where the coverage is lower than a certain threshold (e.g. 10). We haven't done it previously with `awk` in order to test the different coverage thresholds in `R`. Removing the non covered positions (done previously with [awk](#filter-the-data-and-calculate-the-methylation-level)) can be done at the beginning because they are not informative. 

We will use now the `dplyr` library to filter the data.


```r
# select only the rows where the coverage is higher than 10
CG_coverage_filtered = CG %>% filter(coverage > 10)
```

If, as commonly happens, the number of Cs with methylation values = 0 is extremely high, the graph may appear compressed and hard to understand on the **_x_** axis. Thus it might be useful to remove the rows where the methylation is 0. This can be done with the following command:

```r
# select only the rows where the methylation is higher than 0 and the coverage is higher than 10
CG_coverage_filtered = CG %>% filter(coverage > 10 & methR > 0)
```

## Repeat now the same analysis for CHG.
{: .no_toc }

```r
# read the input file, which is missing the header
CHG=read.table("arabidopsis_metilome_CHG.txt", stringsAsFactors=F, header=F,sep="\t")

# rename the columns
names(CHG)=c('chr', 'pos', 'strand', 'c', 't', 'context', 'genome_context', 'methylation')

# Add a new column named coverage which include the total coverage
# total coverage is calculated by summing the columns c and t
CHG$coverage = CHG$c + CHG$t

# Add now a new column named `methR` which represent the methylation level calculated in a different way. The value is calculated as % value, rounded to 0 decimal places:
CHG$methR = round(100*CHG$c / CHG$coverage, 0)

# check now the new dataframe
head(CHG)

# The output should look like
   chr pos strand  c  t context genome_context methylation coverage methR
1 Chr3 381      + 15  5     CHG            CTG    0.750000       20    75
2 Chr3 383      - 29  3     CHG            CAG    0.906250       32    91
3 Chr3 430      + 13  8     CHG            CTG    0.619048       21    62
4 Chr3 432      - 14 10     CHG            CAG    0.583333       24    58
5 Chr3 470      +  9  4     CHG            CTG    0.692308       13    69
6 Chr3 472      - 24  2     CHG            CAG    0.923077       26    92

```


Now we can filter the table by removing the rows where the coverage is lower than a certain threshold (e.g. 10). We haven't done it previously with `awk` in order to test the different coverage thresholds in `R`. Removing the non covered positions (done previously with [awk](#filter-the-data-and-calculate-the-methylation-level)) can be done at the beginning because they are not informative. 

We will use now the `dplyr` library to filter the data.


```r
# select only the rows where the coverage is higher than 5
CHG_coverage_filtered = CHG %>% filter(coverage > 5)
```

If, as commonly happens, the number of Cs with methylation values = 0 is extremely high, the graph may appear compressed and hard to understand on the **_x_** axis. Thus it might be useful to remove the rows where the methylation is 0. This can be done with the following command:

```r
# select only the rows where the methylation is higher than 0 and the coverage is higher than 5
CHG_coverage_filtered = CHG %>% filter(coverage > 5 & methR > 0)
```






## Repeat now the same analysis for CHH.
{: .no_toc }

```r
# read the input file, which is missing the header
CHH=read.table("arabidopsis_metilome_CHH.txt", stringsAsFactors=F, header=F,sep="\t")

# rename the columns
names(CHH)=c('chr', 'pos', 'strand', 'c', 't', 'context', 'genome_context', 'methylation')

# Add a new column named coverage which include the total coverage
# total coverage is calculated by summing the columns c and t
CHH$coverage = CHH$c + CHH$t

# Add now a new column named `methR` which represent the methylation level calculated in a different way. The value is calculated as % value, rounded to 0 decimal places:
CHH$methR = round(100*CHH$c / CHH$coverage, 0)

# check now the new dataframe
head(CHH)

# The output should look like
   chr pos strand c t context genome_context methylation coverage methR
1 Chr3 166      + 0 1     CHH            CCC           0        1     0
2 Chr3 167      + 0 1     CHH            CCT           0        1     0
3 Chr3 168      + 1 0     CHH            CTA           1        1   100
4 Chr3 174      + 0 1     CHH            CCA           0        1     0
5 Chr3 175      + 0 1     CHH            CAT           0        1     0
6 Chr3 181      + 0 1     CHH            CCT           0        1     0

```



Now we can filter the table by removing the rows where the coverage is lower than a certain threshold (e.g. 10). We haven't done it previously with `awk` in order to test the different coverage thresholds in `R`. Removing the non covered positions (done previously with [awk](#filter-the-data-and-calculate-the-methylation-level)) can be done at the beginning because they are not informative. 

We will use now the `dplyr` library to filter the data.


```r
# select only the rows where the coverage is higher than 5
CHH_coverage_filtered = CHH %>% filter(coverage > 5)
```

If, as commonly happens, the number of Cs with methylation values = 0 is extremely high, the graph may appear compressed and hard to understand on the **_x_** axis. Thus it might be useful to remove the rows where the methylation is 0. This can be done with the following command:


```r
# select only the rows where the methylation is higher than 0 and the coverage is higher than 5
CHH_coverage_filtered = CHH %>% filter(coverage > 5 & methR > 0)
```

# 3. Draw the methylation distribution in `R`.

We will use `ggplot2` in order to draw the plot.

```r
# load the required packages 
library(ggplot2)
```

### Draw the graph as histogram:
{: .no_toc }
ggplot use the filtered dataset to draw the histogram. The `geom_histogram` function is used to draw the histogram. The `fill` argument is used to specify the color of the bars. 

```r
ggplot(CG_coverage_filtered,aes(x=methR)) +
geom_histogram(colour=4,fill="white",binwidth=1)
```

### Draw the graph as density plot:
{: .no_toc }

```r
ggplot(CG_coverage_filtered,aes(x=methR)) +
geom_density(alpha=.2,fill="#FF6666")
```

The obtained graphs should look like: 

| **histogram** | **density** |
|:--------:|:---:|
|  ![histo]({{"/assets/images/3a3-3_methylation_distribution_arabidopsis.png" | relative_url }})  | ![density]({{"/assets/images/3a3-4_methylation_distribution_arabidopsis.png" | relative_url }})  |


## Repeat now the same for CHG
{: .no_toc }

### Draw the graph as histogram:
{: .no_toc }

```r
ggplot(CHG_coverage_filtered,aes(x=methR)) +
geom_histogram(colour=4,fill="white",binwidth=1)
```

### Draw the graph as density plot:
{: .no_toc }

```r
ggplot(CHG_coverage_filtered,aes(x=methR)) +
geom_density(alpha=.2,fill="#FF6666")
```


The obtained graphs for CHG should look like: 

| **histogram** | **density** |
|:--------:|:---:|
|  ![histo]({{"/assets/images/3a3-5_methylation_distribution_arabidopsis.png" | relative_url }})  | ![density]({{"/assets/images/3a3-6_methylation_distribution_arabidopsis.png" | relative_url }})  |


## Repeat now the same for CHH
{: .no_toc }


### Draw the graph as histogram:
{: .no_toc }

```r
ggplot(CHH_coverage_filtered,aes(x=methR)) +
geom_histogram(colour=4,fill="white",binwidth=1)
```

### Draw the graph as density plot:
{: .no_toc }

```r
ggplot(CHH_coverage_filtered,aes(x=methR)) +
geom_density(alpha=.2,fill="#FF6666")
```

The obtained graphs for CHH should look like: 

| **histogram** | **density** |
|:--------:|:---:|
|  ![histo]({{"/assets/images/3a3-7_methylation_distribution_arabidopsis.png" | relative_url }})  | ![density]({{"/assets/images/3a3-8_methylation_distribution_arabidopsis.png" | relative_url }})  |


{: .note}
Until now we have only seen the graphs in an interactive way, but what if we want to save them?

We can save the graph in a pdf file by using the `pdf()` function. 

```r
# open the pdf device (already present in r-base)
pdf("CG_density.pdf",paper="A4")

# draw the plot with ggplot
ggplot(CG_coverage_filtered,aes(x=methR)) +
geom_density(alpha=.2,fill="#FF6666")

# close the device
dev.off()
```

In this way we will create a pdf document with a single page (A4) where the plot of the CG density is reported.
We can save all the plots (`CG`, `CHG` and `CHH` contexts) in a single file (as separate pages) by combining multiple plot operations:

```r
# open the pdf device 
pdf("arabidopsis_methylation_density.pdf",paper="A4")
# add the plot of CG context
ggplot(CG_coverage_filtered,aes(x=methR)) +
geom_density(alpha=.2,fill="#FF6666")

# add the plot of CHG context
ggplot(CHG_coverage_filtered,aes(x=methR)) +
geom_density(alpha=.2,fill="#FF6666")

# add the plot of CHH context
ggplot(CHH_coverage_filtered,aes(x=methR)) +
geom_density(alpha=.2,fill="#FF6666")

# close the device
dev.off()
```

{: .note}
>We need to pair `pdf()` and `dev.off()` in order to open (draw the graph) and close the pdf file. We can also use `png()` or `jpeg()` to save plots in a png or jpeg format. With `png()` or `jpeg()` we need an extra step to arrange the figures in multiple panels.

