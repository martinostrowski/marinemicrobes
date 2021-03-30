#!/bin/bash
 
#PBS -N emp1PHB
#PBS -l ncpus=18
#PBS -l mem=128Gb
#PBS -l walltime=12:00:00

echo "Job ID is ${PBS_JOB_ID}"
echo "Job Array ID is ${PBS_ARRAY_INDEX}"
echo "Timestamp is $(date +%F_%T)"
echo "Directory is $(pwd)"
echo "Running on host $(hostname)"
echo "Working directory is ${PBS_O_WORKDIR}"
echo "Job has the following nodes/cores:"
cat ${PBS_NODEFILE}

cd /shared/c3/projects/functions.nrs/

module load devel/c8/miniconda3

module load bio/c8/diamond-2.0.2

module load bio/c8/eggnog-mapper-2.1.0

if [[ ! -e "/scratch/u119966_${PBS_ARRAY_INDEX}/eggnog_proteins.dmnd" ]];     
  then
       mkdir  /scratch/u119966_${PBS_ARRAY_INDEX}

       cp  /shared/c3/bio_db/eggnog_mapper-2.1.0/eggnog.db /shared/c3/bio_db/eggnog_mapper-2.1.0/eggnog_proteins.dmnd /scratch/u119966_${PBS_ARRAY_INDEX}

fi

date +%F_%T


chunk=$(awk -v line=${PBS_ARRAY_INDEX} '{if (NR == line) { print $0; };}' /shared/c3/projects/functions.nrs/MAI_chunks)

emapper.py  --temp_dir /scratch/u119966_${PBS_ARRAY_INDEX} --data_dir /scratch/u119966_${PBS_ARRAY_INDEX}  --scratch_dir /scratch/u119966_${PBS_ARRAY_INDEX} -m diamond --no_annot --no_file_comments --data_dir /scratch/u119966_eggnogdb/  --cpu 18 -i /shared/c3/projects/functions.nrs/${chunk} -o /scratch/u119966_${PBS_ARRAY_INDEX}/${chunk}_EN2.1.0 --override

rm /scratch/u119966_${PBS_ARRAY_INDEX}/eggnog_proteins.dmnd

cp  /scratch/u119966_${PBS_ARRAY_INDEX}/*   /shared/c3/projects/functions.nrs/

rm -r /scratch/u119966_${PBS_ARRAY_INDEX}

date +%F_%T
