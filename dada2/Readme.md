# Sequence variants, denoising and the large scale Microbiome Initiative

Single base pair accuracy is important for addressing fundamental questions in ecology

Large scale sequencing projects seek to describe the biodiversity of (micro)-organisms across a broad range of environmental and host associated systems. Continual improvements in the quality, length and quantity of environmental sequences produced with contemporary methods now offer the opportunity to integrate data obtained from different studies and sequencing runs to provide comprehensive, standardised data sets that cover much larger spatial and temporal scales. The development of controls and several 'denoising' methods to correct and potentially eliminate errors introduced during the amplification and sequencing processes also provide a significant boost to the taxonomic resolution of single marker gene surveys. such that 'Actual Sequence Variants' of true biological origin can be inferred from 

Often involve Dna extracted from a wide range of different samples, using different procedures and different operators, and sequencing data combined across multiple sequencing runs, potentially produced by different sequencing centres, using different primers for marker gene surveys and different post-sequencing analyses pipelines - all of these factors can impact upon the reproducibility and of the results and the interpretation of the data. This study we examine the impact of different filtering and denoising methods on the reported Sequence Variant Tables by comparing a version of the [Unoise3 workflow](https://data.bioplatforms.com/dataset/c532203a-bd1f-4565-bbb7-df3cad5f53c5/resource/92f6a989-bb6f-4005-9b51-0620b431af9f/download/sequence_analysis.pdf)*, with an 'unfiltered' dataset that has not been denoised or screened for chimeras, and the output of a DADA2 workflow.

This repository contains code for a DADA2 reanalysis of the AMMBI and Marine Microbes amplicon datasets. The code is largely based on [DADA2 tutorials](https://benjjneb.github.io/dada2/tutorial_1_8.html) and the work of Dr Anna Bramucci to implement the DADA2 workflow on the UTS HPCC.

The key advantages of Dada2 over other approaches (e.g. Accuracy, Interoperability, Scaling and Open Source) are highlighted on the [DADA2 Github page](https://benjjneb.github.io/dada2/index.html), along with many useful and up-to-date ressources for the statistical analyses of sequencing data post processing.

By comparing these outputs of these three analyses workflows we should be able to determine:

* how well amplicon datasets represent the true diversity within a sample
* potentially identify any systematic artefacts
* obtain a better estimate of the microbial sequence diversity across the Australian datasets and begin to compare directly with global efforts
* develop open source resources to assist with the statistical analyses of the data, or suybsets corresponding to lineages of interest to non-microbiologist researchers


1. 18S Coastal Primary Analysis pipeline
1. [16S Primary Analysis pipeline](./16S/do-dada2f.1.r)
1. [18S Pelagic Primary Analysis pipeline](./18s/do-dada2f.1.r)
1. [A16S Pelagic Primary Analysis pipeline](./a16s/do-dada2f.1.r)
