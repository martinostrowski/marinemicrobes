---
title: "R Notebook"
output: html_notebook
---

# Authors: Martin Ostrowski, Claudia 
## Contact: martin.ostrowski@uts.edu.au

### Description

This notebook outlines The process is required for assembling the managing marine metagenomes produced as part of the Australian Microbiome initiative. There are ~800 metagenomes. Each metagenome consists of approximately 90 million illlumina paired end sequences.

The general steps are as follows

1. Concatenation of files from the same Sample (into R1 and R2 paired ends)
2. Trimmomatic quality filtering and tag removal
  a. generate a merged fasta file
  b. remove redundant reads prior to assembly
3. Digital normalisation with bbnorm and error correction (This is to remove the duplicated sequences and make their files smaller in order to assemble them with available resources.)
4. Metagenomic assembly using spades or megahit
5. Determine the assembly quality and collate statistics using meta quast, assembly_stats,
6. Hand-off to annotation or metagenomic binning


 -------------------------------------------------------

Source data: downloaded from the [Australian Microbiome data portal](https://data.bioplatforms.com/organization/australian-microbiome)

Tools: [Trimmomatic version ](https://github.com/timflutre/trimmomatic)
[bbnorm](bbtools)
[SPAdes  ](http://cab.spbu.ru/software/spades/)


Working Directories:


Repositories:

UTS HPC: /shared/c3/bio_db/BPA/metaG
        /shared/c3/bio_db/BPA/amplicon
        /shared/c3/bio_db/BPA/assembly
        

Mainly a Business Dropbox - but will begin moving data to a Cloudstor repository



Inputs:


processing code

for i in *fasta; do  name=$(basename $i .norm.contigs.fasta); sed -i "s/NODE/"${name}"/" $i; done


### Trimmomatic



### BBnorm 

```r
#!/bin/bash

#PBS -N bbnorm
#PBS -l ncpus=10
#PBS -l mem=128GB
#PBS -l walltime=5:00:00


cd /shared/c3/projects/ami_assembly_2020.martin/
echo "Job ID is ${PBS_JOB_ID}"
echo "Job Array ID is ${PBS_ARRAY_INDEX}"
echo "Timestamp is $(date +%F_%T)"
echo "Directory is $(pwd)"
echo "Running on host $(hostname)"
echo "Working directory is ${PBS_O_WORKDIR}"
echo "Job has the following nodes/cores:"
cat ${PBS_NODEFILE}

#split the input fasta into many parts then distribute widely
CODE=$(awk -v line=${PBS_ARRAY_INDEX} '{if (NR == line) { print $0; };}' /shared/c3/projects/ami_assembly_2020.martin/bbnorm.sort.conf)

date +%F_%T
#echo "file working on is $CODE"
/shared/c3/projects/ami_assembly_2020.martin/bbmap/bbnorm.sh target=30 mindepth=2 in=/shared/c3/bio_db/BPA/metaG/fastq_trimmed/$CODE.R1.fastq.gz in2=/shared/c3/bio_db/BPA/metaG/fastq_trimmed/$CODE.R2.fastq.gz out=/shared/c3/projects/ami_assembly_2020.martin/$CODE.bbnorm.R1.fastq.gz out2=/shared/c3/projects/ami_assembly_2020.martin/$CODE.bbnorm.R2.fastq.gz thread$
date +%F_%T
```
### Spades Assemblies

```r
#!/bin/bash

#PBS -N spades
#PBS -l ncpus=20
#PBS -l mem=240GB
#PBS -l walltime=12:00:00


module load bio/SPAdes-3.14.1

cd /shared/c3/projects/ami_assembly_2020.martin/
echo "Job ID is ${PBS_JOB_ID}"
echo "Job Array ID is ${PBS_ARRAY_INDEX}"
echo "Timestamp is $(date +%F_%T)"
echo "Directory is $(pwd)"
echo "Running on host $(hostname)"
echo "Working directory is ${PBS_O_WORKDIR}"
echo "Job has the following nodes/cores:"
cat ${PBS_NODEFILE}

#split the input fasta into many parts then distribute widely
CODE=$(awk -v line=${PBS_ARRAY_INDEX} '{if (NR == line) { print $0; };}' /shared/c3/projects/ami_assembly_2020.martin/runList3)

#CODE=138100

date +%F_%T
#echo "file working on is $PARAMETERS"

spades.py --tmp-dir /scratch/work/$CODE.spadestmp --meta --only-assembler -k 33,55,77 -m 240 -t 20 -1 /shared/c3/projects/ami_assembly_2020.martin/$CODE.bbnorm.R1.fastq.gz -2 /shared/c3/projects/ami_assembly_2020.martin/$CODE.bbnorm.R2.fastq.gz -o /shared/c3/projects/ami_assembly_2020.martin/$CODE.norm

rm -r /scratch/work/$CODE.spadestmp

date +%F_%T
```
