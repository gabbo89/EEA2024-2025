---
layout: default
title: File formats
parent: 2. General Guides
nav_order: 1
published: true
---

FINAL VERSION
{: .label .label-green }

A brief introduction to various file formats used in bioinformatics.

<br>
<details open markdown="block">
  <summary>
    <strong>Table of contents</strong>
  </summary>
  {: .text-delta }
- TOC
{:toc}
</details>

<!--
## Table of contents

- [Sequence formats](#sequence-formats)
    - [FASTA](#fasta)
	- [Example FASTA file](#example-fasta-file)
	- [FASTA filename extensions](#fasta-filename-extensions)
    - [FASTQ](#fastq)
	- [Example FASTQ file](#example-fastq-file)
- [Alignment formats](#alignment-formats)
    - [SAM](#sam)
    - [BAM](#bam)
    - [CRAM](#cram)
    - [Stockholm format](#stockholm-format)
	- [Example Stockholm file](#example-stockholm-file)
    - [VCF](#vcf)
    - [BCF](#bcf)
- [Generic Feature Formats](#generic-feature-formats)
    - [GFF general structure](#gff-general-structure)
    - [GTF](#gtf)
	- [Example GTF file](#example-gtf-file)
    - [GFF3](#gff3)
	- [Example GFF3 file](#example-gff3-file)
- [References](#references)
-->


# Sequence formats

These are file formats for storing nucleotide sequences and/or amino acid (protein) sequences.

## FASTA

FASTA is a ubiquitous text-based format for representing nucleotide sequences or amino acid sequences. A FASTA file can contain one sequence or multiple sequences. If a FASTA file contains multiple sequences, it may sometimes be referred to as a "multi-FASTA" file.

Each FASTA entry begins with a `>` (greater-than) symbol, followed by a comment on the same line describing the sequence that will follow. The actual sequence begins on the line after this comment. Another `>` (greater-than) symbol denotes the beginning of another FASTA entry (comment describing the sequence + the sequence itself).

##### Example FASTA file
{: .no_toc }
This example FASTA file contains two linear nucleotide sequences.

```
>gi|1817694395|ref|NZ_JAAGMU010000151.1| Streptomyces sp. SID7958 contig-52000002, whole genome shotgun sequence
CCGGCTGGCGCGGCTGGCGCTGGCGGTGGGGCTGCGGCTGCTGGAGCTGGGGGTGGCGCTGGAGGCGCAC
GGCCAGAACCTGCTGGTGGTGCTGTCGCCGTCCGGGGAGCCGCGGCGGCTGGTCTACCGCGATCTGGCGG
ACATCCGGGTCTCCCCCGCGCGGCTGGCCCGGCACGGTATCCGGGTTCCGGACCTGCCGGCG

>gi|1643051563|gb|SZWM01000399.1| Citrobacter sp. TBCS-14 contig3128, whole genome shotgun sequence
GCACAGTGAGATCAGCATTCCGTTGGATCTACTGGTCAATCAAAACCTGACGCTGGGTACTGAATGGAAC
CAGCAGCGCATGAAGGACATGCTGTCTAACTCGCAGACCTTTATGGGCGGTAATATTCCAGGCTACAGCA
GCACCGATCGCAGCCCATATTCGAAAGCCGAGATCTTCTCTTTGTTTGCCGAAAACAACATG
```

##### FASTA filename extensions
{: .no_toc }

FASTA files usually end with the extension `.fasta`. This extension is arbitrary, as the content of the file determines its format, not its extension. More descriptive filename extensions can be used instead of `.fasta`, which are useful as they describe the type of sequence(s) in the file at a glance. 

Here are some examples...
- **`.fna`** can be used for **F**ASTA **n**ucleic **a**cids
- **`.faa`** can be used for **F**ASTA **a**mino **a**cids
- **`.frn`** can be used for **F**ASTA non-coding **RN**A

## FASTQ

The FASTQ format is an extension of [FASTA](#fasta) that stores both biological sequences (usually nucleotide sequences) and their corresponding quality scores. Both the sequence letter and quality score are encoded with a single character for brevity.

A FASTQ file normally uses four lines per sequence:
1. A line beginning with `@` followed by a sequence identifier and optional description (like the comment line at in a FASTA file)
2. The raw sequence letters
3. A line beginning with `+`, sometimes followed by the same comment as the first line
4. A line encoding the quality values for the sequence in line 2, with the same numbers of symbols as letters in the sequence

##### Example FASTQ file from sequencing center
{: .no_toc }

![fastq_example]({{"/assets/images/image-2.png" | relative_url }})

##### Example FASTQ file downloaded from an online database
{: .no_toc}

```
@SRR8933535.1 1 length=75
NAGGAAACAAAGGCTTACCCGTTATCATTTCCGCAAGAATGCACCCACACGACCATATATCAATGGATGTGGAGT
+SRR8933535.1 1 length=75
#AAAAEEEEEEEEEEEEEEEEEEEEEAEEEEEAEEEEEEEAEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAE
@SRR8933535.2 2 length=75
NGAGGAGTGGTGGTAGTGTTGCTTGGTGGCAAAGATGTAGTTGGTGGGAAAGCTGAAGTGGTACCGTTGGTTGGA
+SRR8933535.2 2 length=75
#AAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAEEEE
```
The quality values are encoded using the [Sanger Phred quality scores](https://en.wikipedia.org/wiki/Phred_quality_score), which are a way of encoding the probability of a base call being incorrect. The probability of a base call being incorrect is calculated as `10^(-Q/10)`, where `Q` is the quality score. For example, a quality score of `30` corresponds to a probability of 0.1% of the base call being incorrect.

![Phred quality]({{ "/assets/images/1a4-1_phred_quality.png" | relative_url }})

The quality scores are usually encoded as ASCII characters, with higher quality scores corresponding to higher ASCII values. 

![ASCII conversion]({{ "/assets/images/1a4-2_ascii_conversion.png" | relative_url }})


Quality values are of paramount importance in bioinformatics, as they allow us to assess the reliability of the sequencing data. This is particularly important when working with high-throughput sequencing data, where the quality of the data can vary significantly across the sequence. Low quality values may indicate errors during the sequencing process, while high values may indicate accurate base calls.

Fastq files are widely used to store the sequence information obtained from sequencers (as for example Illumina) and are extensively used in bioinformatics in the pipeline for trimming (removing of low quality bases and/or sequencing adapters), quality control, and alignment of the sequences to a reference genome.

# Alignment formats

These are formats for storing alignments of nucleotide or amino acid sequences.
The alignment formats discussed here can be manipulated using [SAMtools](samtools.md).

## SAM

Sequence Alignment/Map (SAM) is a text-based alignment format that supports single- and paired-end reads produced by different sequencing platforms.
It can support short and long reads (up to 128Mbp).
The format has been extended to include unmapped sequences, and it may contain other data such as base-call & alignment qualities.

The SAM format consists of a header and an alignment section.
Headings begin with the `@` symbol, which distinguishes them from the alignment section.
All lines in a SAM file are tab-delimited.
Alignment sections contain 11 mandatory fields, with other fields being optional.
Although the mandatory fields must be present, their value can be a `*` or a `0` depending on the field.
Optional fields are presented as key-value pairs in the format `TAG:TYPE:VALUE`.

The mandatory fields in a SAM file are...

1. **QNAME** Query template name. Reads/segments with the same QNAME are from the same template. A read may occupy multiple alignment lines when its alignment is chimeric or when multiple mappings are given.
2. **FLAG** Bitwise flag (pairing, strand, mate strand, etc.).
3. **RNAME** Reference sequence name.
4. **POS** 1-based leftmost mapping position. The first base in a reference sequence has coordinate 1. `POS` is set to 0 for an unmapped read.
5. **MAPQ** Mapping quality. If equal to 255, the mapping quality is not available.
6. **CIGAR** Extended Consise Idiosyncratic Gapped Alignment Report string. `M` for match/mismatch, `I` for insertion and `D` for deletion compared with the reference, `N` for skipped bases on the reference, `S` for soft clipping, `H` for hard clipping, and `P` for padding.
7. **RNEXT** Reference name of the mate/next read in the template. For the last read, the next read is the first read in the template.
8. **PNEXT** Position of the primary alignment of the mate/next read.
9. **TLEN** Observed template length.
10. **SEQ** Segment sequence.
11. **QUAL** Phred-based base quality (same as the quailty string in [FASTQ](#fastq)) plus 33.

SAM files can be manipulated using [SAMtools](samtools.md).

## BAM

Binary Alignment Map (BAM) is a binary representation of [SAM](#sam), containing the same information in binary format for improved performance.
A position-sorted BAM file can be indexed to allow random access.

## CRAM

CRAM is a sequencing read file format that is highly space efficient by using reference-based compression of sequence data and offers both lossless and lossy modes of compression.
CRAM files are typically 30 to 60% smaller than their [BAM](#bam) equivalents.
CRAM has the following major objectives:

1. Significantly better lossless compression than BAM
2. Full compatibility with BAM
3. Effortless transition to CRAM from using BAM files
4. Support for controlled loss of BAM data


## PAF

PAF is a text format describing the approximate mapping positions between two
set of sequences. PAF is TAB-delimited with each line consisting of the
following predefined fields:


1. **Query sequence name**
2. **Query sequence length**
3. **Query start coordinate (0-based; BED-like-closed)**
4. **Query end coordinate (0-based; BED-like-open)**
5. **Relative strand** (+ on same strand/- if opposite)
6. **Target sequence name**
7. **Target sequence length**
8. **Target start on original strand (0-based)**
9. **Target end on original strand (0-based)**
10. **Number of matching bases**
11. **Number bases, including gaps, in the mapping**
12. **Mapping quality (0-255; 255 for missing)**

Column 11 gives the total number of sequence matches, and gaps in the alignment. Column 10 divided by column 11 gives the BLAST-like alignment identity. PAF may optionally have additional fields in the SAM-like typed key-value format. 

[Back to the ONT tutorial](https://gabbo89.github.io/EEA2024-2025/docs/3a2_ONT_cleaning_and_alignment.html#paf-format)

<!--
### Stockholm format

The Stockholm format is a system for marking up features in a multiple alignment, used by [HMMER](http://hmmer.org/), [Pfam](https://pfam.xfam.org/), and [Rfam](https://rfam.xfam.org/).
The first line of a Stockholm file (`.sto`, or `.stk`) states the format and version identifier, currently `# STOCKHOLM 1.0`.
The header is followed by mark-up lines beginning with `#`.
These mark-up lines can annotate features of the alignment file (`#=GF`, generic per-file annotation), or features of the aligned sequences (`#=GS`, generic per-sequence annotation).
The sequence alignment itself is a series of lines with sequence names (typically in the form `name/start-end`) followed by a space and the aligned sequence.
A line with two forward slashes (`//`) indicates the end of the alignment.

#### Example Stockholm file

This is a Stockholm format alignment for Alpha-haemoglobin stabilising protein [(AHSP) from Pfam](https://pfam.xfam.org/family/PF09236):

```
# STOCKHOLM 1.0
#=GS G1TJ87_RABIT/6-91   AC G1TJ87.1
#=GS H0XCX0_OTOGA/5-91   AC H0XCX0.1
#=GS F6RTV5_HORSE/5-91   AC F6RTV5.2
#=GS AHSP_BOVIN/5-91     AC Q865F8.1
#=GS AHSP_MOUSE/5-91     AC Q9CY02.1
#=GS G3IJS0_CRIGR/5-91   AC G3IJS0.1
#=GS L9KTP8_TUPCH/5-91   AC L9KTP8.1
#=GS G5BYB6_HETGA/5-104  AC G5BYB6.1
#=GS G3T391_LOXAF/5-91   AC G3T391.1
#=GS F7EDV6_ORNAN/6-92   AC F7EDV6.1
#=GS G3WLQ8_SARHA/5-86   AC G3WLQ8.1
#=GS F7DP52_MONDO/4-90   AC F7DP52.1
#=GS H2NS04_PONAB/5-91   AC H2NS04.1
G1TJ87_RABIT/6-91              .TNKDLISMGLKEF....NVLLNQ.........QVFSDPL.LSQEAMQTVLDDWVNLYVNYYRQQMTGEQQELDKALEELRLELNGLAKPFLNKYSVFLKS
H0XCX0_OTOGA/5-91              QANEDLISAGVKEF....NILLNQ.........QVFNEPF.VSEEAMETVVNDWVNFYMNYYKKQMTGEQGEQEKALQELKQKLNSLANPFLAKYRAFLKS
F6RTV5_HORSE/5-91              QANRDLISTAIKEF....NVLLNQ.........QVFSDPP.VSEEAMVTVVNDWVNFYINYYRRQVVGEQQEKDRALQELRQELNILSAPFLAKYRAFLKS
AHSP_BOVIN/5-91                QTNKDLISKGIKEF....NILLNQ.........QVFSDPA.ISEEAMVTVVNDWVSFYINYYKKQLSGEQDEQDKALQEFRQELNTLSASFLDKYRNFLKS
AHSP_MOUSE/5-91                QSNKDLISTGIKEF....NVLLDQ.........QVFDDPL.ISEEDMVIVVHDWVNLYTNYYKKLVHGEQEEQDRAMTEFQQELSTLGSQFLAKYRTFLKS
G3IJS0_CRIGR/5-91              QTNKELISEGIKQF....NVLLGQ.........QVFDDPL.IPEENMVTVVNDWVNLYINYYKPLVFGKQQEQDKALQELQQELNTLGSQFLTKYRTILKS
L9KTP8_TUPCH/5-91              QVNKDIIATGMKKF....SVLLDQ.........QVFSEPP.ISEEAMVVVVNDWVNFYVNYYGQQVTGEQQEQDRALNELRQELTTMASPFLAKYRAFLKS
G5BYB6_HETGA/5-104             QANKDLIALGMKEFPADYSDMLESHSLSPASHPQVFNYPL.ITEEDMVVVVDDWVNIYINYYRKRLTGEKQDQDRALQELRQELKTLASPFLAKYRACLES
G3T391_LOXAF/5-91              QANKDLISTGMKEF....SILLNQ.........QDMRDNP.IPEEAMVIVVNDWMSFYINYYRQKMTGEQQEQDRALQELQQGLNTLANPFLTKYRDFLKT
F7EDV6_ORNAN/6-92              .SNQDVINSAMAAF....QALLNQ.........QVFSPQIPIPMEAMKIIVRDWIEFYISYFAPKLRGDRQERERAQEDLWETLQAIARPFLDKYRDFLNA
G3WLQ8_SARHA/5-86              QSNQDVISSAMQEF....SKLLDQ.........QEFTKPA.FSETDMVTIVDDWIKFYLSYYSKKMTGNEQEQERAMQKLQEELRTSASPFLDKSQ.....
F7DP52_MONDO/4-90              QSNQDVISSAMQEF....NKLLNQ.........QDFTYAV.ISEKDMVTIVDDWMNYYLSFFSQKMSGDQQEQERAMQKLQEELRSSANPFLDKYRAFLKS
H2NS04_PONAB/5-91              KANKDLISAGLKEF....SVLLNQ.........QVFNDPL.ISEEDMVTVVEDWMNFYINYYRQQVTGEPQERDKALQELRQELNTLANPFLAKYRDFLKS
#=GC seq_cons                  QuNKDLISsGhKEF....slLLNQ.........QVFs-Ph.ISEEsMVTVVsDWVNFYlNYY+pploGEQQEQDRALQELpQELsTLAsPFLsKYRsFLKS
//
```
-->

# Variation

## VCF

Variant Call Format (VCF) is a format for storing variations between a reference genome and sequences aligned to it, based on SAM/BAM alignments.
VCF files begin with a header section: lines in the header section begin with `##`.
The last line in the header section begins with `#`; this line gives the headers of the columns used in the VCF file:

1. **CHROM** The name of the sequence (typically a chromosome) on which the variation is being called. This sequence is usually known as 'the reference sequence', i.e. the sequence against which the given sample varies.
2. **POS** The 1-based position of the variation on the given sequence.
3. **ID** The identifier of the variation, e.g. a dbSNP rs identifier, or if unknown a ".". Multiple identifiers should be separated by semi-colons without white-space.
4. **REF** The reference base (or bases in the case of an indel) at the given position on the given reference sequence.
5. **ALT** The list of alternative alleles at this position.
6. **QUAL** A quality score associated with the inference of the given alleles.
7. **FILTER** A flag indicating which of a given set of filters the variation has passed.
8. **INFO** An extensible list of key-value pairs (fields) describing the variation. Multiple fields are separated by semicolons with optional values in the format: <key>=<data>[,data].
9. **FORMAT** An (optional) extensible list of fields for describing the samples.
10. **SAMPLEs** For each (optional) sample described in the file, values are given for the fields listed in FORMAT. If multiple samples have been aligned to the reference sequence, each sample will have its own column.

## BCF

Binary Call Format (BCF) is a binary representation of [VCF](#vcf), containing the same information in binary format for improved performance.

# Annotations

The Generic Feature Formats (`.gff`) are tab-delimited text file formats used for describing genes and other features of DNA, RNA and protein sequences.
GFF files are used to annotate genomes, as they describe functional regions of genomes.

There are two widely used versions of the GFF file format:

1. [Gene Transfer Format](#gtf), a variation of GFF version 2.
2. [Generic Feature Format version 3 (GFF3)](#gff3).

## Generic Feature Format (GFF) general structure

All GFF formats (GFF2, GFF3 and GTF) are tab delimited with 9 fields per line. They all share the same structure for the first 7 fields, while differing in the content and format of the ninth field. The general structure is as follows: 

1. **sequence** The name of the sequence where the feature is located.
2. **source** Keyword identifying the source of the feature, like a program (e.g. Augustus) or an organization (e.g. [SGD](https://www.yeastgenome.org/)).
3. **feature** The feature type name, like `gene` or `exon`. In a well-structured GFF file, all the children features always follow their parents in a single block (so all exons of a transcript are put after their parent `transcript` feature line and before any other parent transcript line).
4. **start** Genomic start of the feature, with a 1-base offset.
5. **end** Genomic end of the feature, with a 1-base offset.
6. **score** Numeric value that generally indicates the confidence of the source in the annotated feature. A value of `.` (a dot) is used to define a null value.
7. **strand** Single character that indicates the strand of the feature; it can assume the values of `+` (positive, or `5'->3'`), `-`, (negative, or `3'->5'`), `.` (undetermined).
8. **phase** Phase of coding sequence (CDS) features, indicating where the feature starts in relation to the reading frame. It can be either one of `0`, `1`, `2` (for CDS features) or `.` (for everything else).
9. **attributes** All the other information pertaining to this feature. The format, structure and content of this field is the one which varies the most between GFF formats.

## GTF

GTF files (`.gtf`) have the following tab-delimited fields (see also [GFF general structure](#gff-general-structure)):

1. **seqname** The name of the sequence. Commonly, this is the chromosome ID or contig ID.
2. **source**
3. **feature** The following feature types are required: `CDS`, `start_codon`, `stop_codon`. The features `5UTR`, `3UTR`, `inter`, `inter_CNS`, `intron_CNS` and `exon` are optional.
4. **start**
5. **end**
6. **score**
7. **strand**
8. **frame** Similar to `phase` in the [GFF general structure](#gff-general-structure).
9. **attributes** All nine features have the same two mandatory attributes at the end of the record: `gene_id "value"` and `transcript_id "value"`. These are globally unique identifiers for the genomic locus of the transcript and the predicted transcript, respectively. If empty, no gene/transcript is associated with this feature. Attributes are separated by `;` (semi-colon).

GTF files also support comments, beginning with `#` and running until the end of the line.
Nothing beyond a hash will be parsed.
These may occur anywhere in the file, including at the end of a feature line. 

##### Example GTF file
{: .no_toc}

Here are the first 12 lines of [example_genome_annotation.gtf](../data/example_genome_annotation.gtf):

```
#gtf-version 2.2
#!genome-build R64
#!genome-build-accession NCBI_Assembly:GCF_000146045.2
#!annotation-source SGD R64-2-1
NC_001133.9	RefSeq	gene	1807	2169	.	-	.	gene_id "YAL068C"; db_xref "GeneID:851229"; gbkey "Gene"; gene "PAU8"; gene_biotype "protein_coding"; locus_tag "YAL068C"; partial "true"; 
NC_001133.9	RefSeq	exon	1807	2169	.	-	.	gene_id "YAL068C"; transcript_id "NM_001180043.1"; db_xref "GeneID:851229"; gbkey "mRNA"; gene "PAU8"; locus_tag "YAL068C"; partial "true"; product "seripauperin PAU8"; exon_number "1"; 
NC_001133.9	RefSeq	CDS	1810	2169	.	-	0	gene_id "YAL068C"; transcript_id "NM_001180043.1"; db_xref "SGD:S000002142"; db_xref "GeneID:851229"; experiment "EXISTENCE:mutant phenotype:GO:0030437 ascospore formation [PMID:12586695]"; gbkey "CDS"; gene "PAU8"; locus_tag "YAL068C"; note "hypothetical protein; member of the seripauperin multigene family encoded mainly in subtelomeric regions"; product "seripauperin PAU8"; protein_id "NP_009332.1"; exon_number "1"; 
NC_001133.9	RefSeq	start_codon	2167	2169	.	-	0	gene_id "YAL068C"; transcript_id "NM_001180043.1"; db_xref "SGD:S000002142"; db_xref "GeneID:851229"; experiment "EXISTENCE:mutant phenotype:GO:0030437 ascospore formation [PMID:12586695]"; gbkey "CDS"; gene "PAU8"; locus_tag "YAL068C"; note "hypothetical protein; member of the seripauperin multigene family encoded mainly in subtelomeric regions"; product "seripauperin PAU8"; protein_id "NP_009332.1"; exon_number "1"; 
NC_001133.9	RefSeq	stop_codon	1807	1809	.	-	0	gene_id "YAL068C"; transcript_id "NM_001180043.1"; db_xref "SGD:S000002142"; db_xref "GeneID:851229"; experiment "EXISTENCE:mutant phenotype:GO:0030437 ascospore formation [PMID:12586695]"; gbkey "CDS"; gene "PAU8"; locus_tag "YAL068C"; note "hypothetical protein; member of the seripauperin multigene family encoded mainly in subtelomeric regions"; product "seripauperin PAU8"; protein_id "NP_009332.1"; exon_number "1"; 
NC_001133.9	RefSeq	gene	2480	2707	.	+	.	gene_id "YAL067W-A"; db_xref "GeneID:1466426"; gbkey "Gene"; gene_biotype "protein_coding"; locus_tag "YAL067W-A"; partial "true"; 
NC_001133.9	RefSeq	exon	2480	2707	.	+	.	gene_id "YAL067W-A"; transcript_id "NM_001184582.1"; db_xref "GeneID:1466426"; gbkey "mRNA"; locus_tag "YAL067W-A"; partial "true"; product "uncharacterized protein"; exon_number "1"; 
NC_001133.9	RefSeq	CDS	2480	2704	.	+	0	gene_id "YAL067W-A"; transcript_id "NM_001184582.1"; db_xref "SGD:S000028593"; db_xref "GeneID:1466426"; gbkey "CDS"; locus_tag "YAL067W-A"; note "hypothetical protein; identified by gene-trapping, microarray-based expression analysis, and genome-wide homology searching"; product "uncharacterized protein"; protein_id "NP_878038.1"; exon_number "1"; 
```

## GFF3

GFF3 files (`.gff3` or `.gff`) have the same tab-delimited fields as [GTF](#gtf), with the following differences:

1. Column/field 3 is referred to as `type` instead of `feature`.
2. In column/field 9, feature attributes are in the format `tag=value`, as opposed to `tag "value"`. Multiple `tag=value` pairs are still separated by semicolons.
3. All attributes that begin with an uppercase letter are reserved for later use. Attributes that begin with a lowercase letter can be used freely by applications.
4. The attribute `Parent` indicates the parent of a feature. A parent ID can be used to group exons into transcripts, transcripts into genes, an so forth. A feature may have multiple parents. `Parent` can only be used to indicate a "part of" relationship (i.e. that a feature is a smaller *part of* a larger feature).

##### Example GFF3 file
{: .no_toc}
From the [specification page for GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md):

```
##gff-version 3.2.1
##sequence-region ctg123 1 1497228
ctg123 . gene            1000  9000  .  +  .  ID=gene00001;Name=EDEN
ctg123 . TF_binding_site 1000  1012  .  +  .  ID=tfbs00001;Parent=gene00001
ctg123 . mRNA            1050  9000  .  +  .  ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123 . mRNA            1050  9000  .  +  .  ID=mRNA00002;Parent=gene00001;Name=EDEN.2
ctg123 . mRNA            1300  9000  .  +  .  ID=mRNA00003;Parent=gene00001;Name=EDEN.3
ctg123 . exon            1300  1500  .  +  .  ID=exon00001;Parent=mRNA00003
ctg123 . exon            1050  1500  .  +  .  ID=exon00002;Parent=mRNA00001,mRNA00002
ctg123 . exon            3000  3902  .  +  .  ID=exon00003;Parent=mRNA00001,mRNA00003
ctg123 . exon            5000  5500  .  +  .  ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123 . exon            7000  9000  .  +  .  ID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123 . CDS             1201  1500  .  +  0  ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
```

## Browse Extensible Data (BED)
BED files (`.bed`) are used to store genomic features in a tab-delimited format. Each line represents a feature with 3-12 columns of data. The first three fields are required:

1. **chromosome**: The name of the chromosome or scaffold where the feature is located.
2. **start**: The starting position of the feature in the chromosome (0-based).
3. **end**: The ending position of the feature in the chromosome.


```
chr1  213941196  213942363
chr1  213942363  213943530
chr1  213943530  213944697
chr2  158364697  158365864
chr2  158365864  158367031
chr3  127477031  127478198
chr3  127478198  127479365
chr3  127479365  127480532
chr3  127480532  127481699
```

{: .note-title}
>Interval types
>
> Different types of intervals exists, which may be classified as follows:
> 
> **FULLY-OPEN** 
> 
>![fully-open]({{ "/assets/images/1a4-3_fully-open.png" | relative_url }})
>
> **FULLY-CLOSED**
>
>![fully-open]({{ "/assets/images/1a4-4_fully-closed.png" | relative_url }})
>
> **HALF-OPEN**
>
>![half-open]({{ "/assets/images/1a4-5_half-open.png" | relative_url }})


{: .note-title}
>Starting types
>
>Different types of start exists, which may be classified as follows:
>
> **0-based**
>
> **1-based**

BED format is **0-start**{: .label .label-blue }, **half-open**{: .label .label-blue }, which means that the start position is inclusive and the end position is exclusive.
`Closed-start` A is included and `Open-end` B is excluded.

![bed-format]({{ "/assets/images/1a4-6_bed_format.png" | relative_url }})

**1-based, fully-closed** system is more intuitive, but it is not always the most efficient method for performing calculations in bioinformatic systems, because an additional step is required to calculate the size of the bp range.


Nine additional fields are optional. Columns cannot be empty, but can be filled with a dot (`.`) if the value is not applicable. The fields are:

4. **name**: A name for the feature.
5. **score**: A numerical value representing the score of the feature.
6. **strand**: The strand where the feature is located, either `+` or `-`.
7. **thickStart**: The starting position at which the feature is drawn thickly.
8. **thickEnd**: The ending position at which the feature is drawn thickly.
9. **itemRgb**: A comma-separated list of RGB values for the feature.
10. **blockCount**: The number of blocks (exons) that make up the feature.
11. **blockSizes**: A comma-separated list of the sizes of the blocks.
12. **blockStarts**: A comma-separated list of the starting positions of the blocks relative to the start of the feature.

##### Example BED file
{: .no_toc}
From the [UCSC Genome Browser](https://genome.ucsc.edu/FAQ/FAQformat.html#format1):


```
track name=pairedReads description="Clone Paired Reads" useScore=1
chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
```

or

```
chr7  127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
chr7  127472363  127473530  Pos2  0  +  127472363  127473530  255,0,0
chr7  127473530  127474697  Pos3  0  +  127473530  127474697  255,0,0
chr7  127474697  127475864  Pos4  0  +  127474697  127475864  255,0,0
chr7  127475864  127477031  Neg1  0  -  127475864  127477031  0,0,255
chr7  127477031  127478198  Neg2  0  -  127477031  127478198  0,0,255
chr7  127478198  127479365  Neg3  0  -  127478198  127479365  0,0,255
chr7  127479365  127480532  Pos5  0  +  127479365  127480532  255,0,0
chr7  127480532  127481699  Neg4  0  -  127480532  127481699  0,0,255
```

## BedGraph 

BedGraph format is a suitable format for storing the density of features along a chromosome. It is a tab-delimited format, based on the BED format with four columns where the score or density is placed in column 4. 
Data are usually preceded by a track definition line which adds a number of options for setting the default display of the track.

```
track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20
chr19 49302000 49302300 -1.0
chr19 49302300 49302600 -0.75
chr19 49302600 49302900 -0.50
chr19 49302900 49303200 -0.25
```



# References
{: .no_toc }
- [SAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
- [CRAM format specification (version 3.0)](https://github.com/samtools/hts-specs/blob/5a5d05fa157c679f34db8920ce3acab1d9f3dfd1/CRAMv3.pdf)
- [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352)
- [GTF2.2: A Gene Annotation Format (Revised Ensembl GTF)](http://mblab.wustl.edu/GTF22.html)
- [GFF3 Specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
- [Stockholm format](https://sonnhammer.sbc.su.se/Stockholm.html)
- [Specifications of SAM/BAM and related high-throughput sequencing file formats](https://github.com/samtools/hts-specs)
- [VCF poster](http://vcftools.sourceforge.net/VCF-poster.pdf)