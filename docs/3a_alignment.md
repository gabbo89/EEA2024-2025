---
layout: default
title: Lezione 1 - raw data cleaning
nav_order: 1
parent: 3. Tutorial
description: A comprehensive guide to understanding epigenetics.
---

#### perform deduplicate 
```deduplicate_bismark --bam rkatsiteli.leaves.rkatsiteli.leaves.R1_bismark_bt2_pe.bam ```


### extract methylation information 
```bismark_methylation_extractor --genome_folder /projects/novabreed/share/gmagris/collaboration/lezioni/2024/EEA/reference/ -p --bedGraph --cytosine_report --CX_context --multicore 1 --gzip rkatsiteli.leaves.rkatsiteli.leaves.R1_bismark_bt2_pe.deduplicated.bam```

<p class="callout success">
Finished writing out cytosine report for covered chromosomes (processed 343 chromosomes/scaffolds in total)

Now processing chromosomes that were not covered by any methylation calls in the coverage file...
Writing cytosine report for not covered chromosome h1tg000079l
Writing cytosine report for not covered chromosome h1tg000179l
Finished writing out cytosine report (processed 345 chromosomes/scaffolds in total). coverage2cytosine processing complete.



Finished generating genome-wide cytosine report
</p>

