# Sequence variants, denoising and the large scale Microbiome Initiatives: How important is a single base pair for addressing fundamental questions in ecology?

Large scale sequencing projects seek to describe the biodiversity of (micro)-organisms across a  broad range of environmental and host associated systems. Often involve Dna extracted from a wide range of different samples, using different procedures and different operators, and sequencing data combined across multiple sequencing runs, potentially produced by different sequencing centres, using different primers for marker gene surveys and different post-sequencing analyses pipelines - all of these factors can impact upon the reproducibility and of the results and the interpretation of the data. This study we examine the impact of different filtering and denoising methods on the reported Sequence Variant Tables by comparing a version of the [Unoise3 workflow](https://data.bioplatforms.com/dataset/c532203a-bd1f-4565-bbb7-df3cad5f53c5/resource/92f6a989-bb6f-4005-9b51-0620b431af9f/download/sequence_analysis.pdf)*, an 'unfiltered' dataset that has no been denoised or screened for chimeras, and the output of a DADA2 workflow. N.B  **This seems to be an old version of the BASE workflow** implemented by the AMI data team because there is no mention of primer removal or Unoise.

As we integrate new molecular dimensions into our understanding of complex living system is edging us towards a better grasp of the universal rules that underpin (microbial) ecology. The interpretation of fundamental questions, such as ‘How many species are there?’ and ‘what processes contribute to the structure microbial communities.

This repository contains code for a DADA2 reanalysis of the AMMBI and Marine Microbes amplicon datasets. The code is largely based on [DADA2 tutorials](https://benjjneb.github.io/dada2/tutorial_1_8.html) and the work of Dr Anna Bramucci to implement the DADA2 workflow on the UTS HPCC.

The key advantages of Dada2 over other approaches (e.g. Accuracy, Interoperability, Scaling and Open Source) are highlighted on the [DADA2 Github page](https://benjjneb.github.io/dada2/index.html), along with many useful and up-to-date ressources for the statistical analyses of sequencing data post processing.

By comparing these outputs of these three analyses workflows we should be able to determine:

* how well amplicon datasets represent the true diversity within a sample
* potentially identify any systematic artefacts
* obtain a better estimate of the microbial sequence diversity across the Australian datasets and begin to compare directly with global efforts
* develop open source resources to assist with the statistical analyses of the data, or suybsets corresponding to lineages of interest to non-microbiologist researchers


1. 18S Coastal Primary Analysis pipeline



Expected impacts on Beta Diversity

The relative abundance of particular sequence variants can be an important indicator of the combined status, status defined as the collective physical, chemiscal, biogeochemical and biological drivers, and the history of an ecosystem. 

WHat are all of the sources of sequencing errors? what assumptions are made in combating them
what are the impacts on betadiversity?
what does that prevent us from doing?
When sequence error correction Unsupervised 
