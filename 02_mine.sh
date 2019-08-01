
#!/bin/bash


# qsub script to run MINE on a matrix of 10000 columns 

#PBS -N MINE.java
#PBS -l nodes=1:ppn=1,mem=10GB
#PBS -l walltime=00:60:00
#PBS -j oe
#PBS -t 1-10000
cd ${PBS_O_WORKDIR}

java -jar /home/martino/MINEv2.jar AMI.table.for.modelling.201907.csv ${PBS_ARRAYID} cv=0.1  -equitability  id=AMI.spread.${PBS_ARRAYID}

