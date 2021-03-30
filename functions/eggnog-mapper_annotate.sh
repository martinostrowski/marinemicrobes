
#!/bin/bash
 
#PBS -N emp2PHB
#PBS -l ncpus=20
#PBS -l mem=32Gb
#PBS -l walltime=12:00:00
        



echo "Job ID is ${PBS_JOB_ID}"
echo "Job Array ID is ${PBS_ARRAY_INDEX}"
echo "Timestamp is $(date +%F_%T)"
echo "Directory is $(pwd)"
echo "Running on host $(hostname)"
echo "Working directory is ${PBS_O_WORKDIR}"
echo "Job has the following nodes/cores:"
cat ${PBS_NODEFILE}

module load devel/c8/miniconda3

module load bio/c8/diamond-2.0.2
module load bio/c8/eggnog-mapper-2.1.0

cd /shared/c3/projects/in2016v04.ami/;

if [[ ! -e "/scratch/u119966_${PBS_ARRAY_INDEX}/eggnog.db" ]];     
  then
       mkdir  /scratch/u119966_${PBS_ARRAY_INDEX}

       cp  /shared/c3/bio_db/eggnog_mapper-2.1.0/eggnog.db  /scratch/u119966_${PBS_ARRAY_INDEX}

fi


line=$(awk -v line=${PBS_ARRAY_INDEX} '{if (NR == line) { print $0; };}' /shared/c3/projects/functions.nrs/eggnog_search.chunkx)

mkdir  /scratch/u119966_p2${PBS_ARRAY_INDEX}




emapper.py  --scratch_dir /scratch/u119966_p2${PBS_ARRAY_INDEX} --data_dir /scratch/u119966_${PBS_ARRAY_INDEX} --annotate_hits_table /shared/c3/projects/functions.nrs/${line} --no_file_comments -o /shared/c3/projects/functions.nrs/${line} --decorate_gff yes --cpu 20 --override

rm /scratch/u119966_${PBS_ARRAY_INDEX}/eggnog.db

cp /scratch/u119966_p2${PBS_ARRAY_INDEX}/*  /shared/c3/projects/functions.nrs

#rm /dev/shm/eggnog.db



rm -rf /scratch/u119966_p2${PBS_ARRAY_INDEX}


date +%F_%T

