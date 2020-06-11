---
title: "R Notebook"
output: html_notebook
---

# Authors: Martin Ostrowski, Claudia 
## Contact: martin.ostrowski@uts.edu.au

### Description

This notebook outlines The process is required for assembling the managing the metagenomes produced as part of the Australian Microbiome initiative. In total we expect approximately 200 new metagenomes to be added to a collection of 570 existing metagenomes. Each metagenome consists of approximately 90 million reads of illlumina Novaseq fastq sequences.

The general steps are as follows

1. Concatenation of files from the same Sample (into R1 and R2 paired ends)
2. Trimmomatic quality filtering and tag removal
  a. generate a merged fasta file
3. Digital normalisation with bbnorm and error correction (This is to remove the duplicated sequences and make their files smaller in order to assemble them with available resources.)
4. Metagenomic assembly using spades or megahit
5. Determine the assembly quality and collate statistics using meta quast, assembly_stats,
6. Hand-off to annotation or metagenomic binning


# -------------------------------------------------------
Source data: downloaded from the [Australian Microbiome data portal](https://data.bioplatforms.com/organization/australian-microbiome)


Working Directories:


Repositories:

UTS HPC: /shared/c3/bio_db/BPA/metaG
        /shared/c3/bio_db/BPA/amplicon
        /shared/c3/bio_db/BPA/assembly
        

Mainly a Business Dropbox - but will begin moving data to a Cloudstor repository



Inputs:

