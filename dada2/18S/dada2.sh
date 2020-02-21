
	
#!/bin/bash
 
#PBS -N DADA2.p2
#PBS -l ncpus=20
#PBS -l mem=64GB
#PBS -l walltime=12:00:00
#PBS -e /shared/c3/bio_db/BPA/20.err
#PBS -o /shared/c3/bio_db/BPA/20.out
#PBS -q c3highmem


module load devel/R-current


cd /shared/c3/bio_db/BPA/18S
echo "Job ID is ${PBS_JOBID}"
echo "Job Array ID is ${PBS_ARRAY_INDEX}"
echo "Timestamp is $(date +%F_%T)"
echo "Directory is $(pwd)"
echo "Running on host $(hostname)"
echo "Working directory is ${PBS_O_WORKDIR}"
echo "Job has the following nodes/cores:"
cat ${PBS_NODEFILE}


#PARAMETERS=$(awk -v line=${PBS_ARRAY_INDEX} '{if (NR == line) { print $0; };}' the.conf)

date +%F_%T


echo "$PBS_ARRAY_INDEX"

Rscript --verbose do-dada2f.r ${PBS_ARRAY_INDEX} > fd.${PBS_ARRAY_INDEX}.out;
