---
layout: default
title: Lesson 7 - Chip-seq analysis
parent: 3. Tutorial
nav_order: 8
description: A comprehensive guide to understanding epigenetics.
published: true
---

INCOMPLETE
{: .label .label-red }

{: .important-title }
> Aim
>
> Perform a Chipseq analysis

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

By combining chromatin immunoprecipitation (ChIP) assays with sequencing, ChIP sequencing (ChIP-Seq) is a powerful method for identifying genome-wide DNA binding sites for transcription factors and other proteins.

In this technique, we first cross-link chromatin complexes, isolate them from the cell nuclei and fragment them. We can then purify chromatin fragments containing our protein of interest by immunoprecipitation. After this, the DNA fragments are purified and sequenced. We can use the sequencing results to determine the DNA regions our protein of interest interacts with.

After isolation and fragmentation of chromatin, the protein–DNA complexes are captured using antibodies specific to the histone or transcription factor of interest. After reversal of crosslinks, the ChIP DNA is then purified and analyzed  by high-throughput sequencing (ChIP-seq).

![alt text](image-63.png)


Our epigenome is defined as methylated DNA and modified histone proteins (around which both methylated and unmethylated DNA are wrapped). DNA methylation and histone modifications undergo global changes during transitions in developmental states and in diseases such as cancer and therefore are major contributors to the dynamic nature of chromatin.

![alt text](image-61.png)
**Figure 1:** Schematic nucleosome structure. A nucleosome consists of two copies of each core histone (H2A, H2B, H3 and H4) and ∼150 bp DNA. The N-terminal tail of each histone is extruded from the nucleosome.


Specific histone modifications may still serve as good epigenetic indicators of chromatin state, even if they are not directly involved in the regulation of gene expression. Trimethylated H3K4 (H3K4me3) is a good marker associated with actively transcribed genes. In active genes, H3K4me3 is enriched around transcription start sites (TSS) ([Figure 2](#figure-2)), whereas H3K4me2 peaks just downstream, followed by monomethylated H3K4 (H3K4me1) further downstream towards the gene body. To summarize, TSSs of actively transcribed genes are marked by H3K4me3 and H3K27ac, and active enhancers can be identified by enrichments of both H3K4me1 and H3K27ac.

![Figure2](image-62.png)

**Figure 2: **Distributions of six modifications with respect to genes are schematically illustrated. TSS, transcription start site; TES, transcription end site. H3K4me3 is enriched around TSSs. H3K4me1 is enriched around enhancers and more downstream. H3K27ac is enriched around active enhancers and TSSs. In undifferentiated stem cells, both H3K4me3 and H3K27me3 (active and inactive marks, respectively) are enriched around TSSs on many genes. H3K27me3 is enriched around inactive TSS in somatic cells. H3K9me3 is broadly distributed on inactive regions. H3K27me3 and H3K9me3 are usually not colocalized. TSSs are generally devoid of nucleosomes.


We need to add a control sample to our ChIP-seq experiment to account for non-specific binding of the antibody. Two commonly controls are used:

1. DNA isolated from the same cells but without immunoprecipitation (input DNA). Cells are cross-linked and fragmented, but no immunoprecipitation is performed. This DNA is then sequenced and used as a control to account for non-specific binding of the antibody [Control for subtracting the effect of chromatin accessibility]

2. Performing  a ChIP-seq experiment with an antibody that does not bind to the protein of interest (not known to be involved in DNA bindin or chromatin modifications, such as IgG). This is called an (Mock IP). [Control for antibody specificity] 


#

### Activate the conda environment
{: .no_toc }

```bash
conda activate chipseq_data
```


```bash
# Create a new directory for this tutorial
mkdir -p /data2/student_space/st24_16_folder/epigenomics/chip_seq/

# Move the working directory
mkdir -p /data2/student_space/st24_16_folder/epigenomics/chip_seq/

# Create the output folders 
mkdir -p fastq trimmed aligned reference bigwig plots
```

We need to create the indexes for bowtie2 which is the aligner we will use to map the reads to the genome.

```bash
# Copy the fasta file 
cp \
/data2/biotecnologie_molecolari_magris/epigenomics/wgbs/reference/vitis_vinifera.fasta reference/

# Create the indexes 
bowtie2-build reference/vitis_vinifera.fasta  reference/vitis_vinifera
```


# 1. H3K4me1 data analysis

```bash
# Define a set of variables
dataset=h3k4me1
R1=fastq/rkatsiteli.${dataset}.R1.fastq.gz
R2=fastq/rkatsiteli.${dataset}.R2.fastq.gz
threads=2

# path to bowtie2 indexes
ref_index=reference/vitis_vinifera
```

### Trimming with fastp
{: .no_toc }

```bash
fastp \
  -i $R1 -I $R2 \
  -w 3 \
  -o trimmed/${dataset}_1.trimmed.fastq.gz \
  -O trimmed/${dataset}_2.trimmed.fastq.gz \
  -h logs/${dataset}_fastp.html \
  -j logs/${dataset}_fastp.json

```

### Alignment with bowtie2
{: .no_toc }

```bash
bowtie2 -x $ref_index \
  -1 trimmed/${dataset}_1.trimmed.fastq.gz \
  -2 trimmed/${dataset}_2.trimmed.fastq.gz \
  -S aligned/${dataset}.sam \
  --very-sensitive \
  -p $threads 2> logs/${dataset}_bowtie2.log
```


### Convert SAM to BAM and remove duplicated reads
{: .no_toc }

```bash
samtools view  -@ ${threads} -b -o aligned/${dataset}.bam aligned/${dataset}.sam 

# remove the sam file (since we have the bam file now)
rm aligned/${dataset}.sam 

# Before sorting the bam file, we will fixmate the reads 
# to ensure that the mate information is correct, 
# the output will be redirect to the sort command and then
# duplicated reads will be marked
# -u options of samtools sort specify an uncompressed output
samtools fixmate \
-@ ${threads} \
-m \
-u \
aligned/${dataset}.bam - | \
samtools sort \
-@ ${threads} \
-u \
- | \
samtools markdup \
-@ ${threads} \
--write-index \
- \
aligned/${dataset}.dedup.bam
```


### Peak calling
{: .no_toc }

```bash
macs3 callpeak \
  -t aligned/${dataset}.dedup.bam \
  -c aligned/rkatsiteli_input.chr5.bam \
  -f AUTO \
  -g hs \
  -n ${SAMPLE}_chr5 \
  --outdir macs2_output \
  --nomodel \
  --shift -100 \
  --extsize 200 \
  -q 0.01

echo "Step 7: Generate bigWig for chr5 with deepTools"
bamCoverage \
  -b aligned/${dataset}.dedup.bam \
  -o bigwig/${dataset}.dedup.bw \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 26902106 \
  --binSize 25 \
  --extendReads 200 \
  -p $threads



# Compute matrix around TSS
computeMatrix reference-point \
  -S bigwig/${SAMPLE}.chr5.bw \
  -R chr5_tss.bed \
  --referencePoint TSS \
  -b 3000 -a 3000 \
  --skipZeros \
  -bs 50 \
  -p $THREADS \
  -out plots/${SAMPLE}_TSS_matrix.gz \
  --outFileSortedRegions plots/${SAMPLE}_TSS_genes.bed

# Plot heatmap
plotHeatmap \
  -m plots/${SAMPLE}_TSS_matrix.gz \
  -out plots/${SAMPLE}_TSS_heatmap.png \
  --colorMap RdBu \
  --heatmapHeight 10 \
  --heatmapWidth 6 \
  --dpi 300
