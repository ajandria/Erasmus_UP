#!/bin/sh

#PBS -N fastqc_fastp
#PBS -q xeon
#PBS -t 1-29
#PBS -o arrayTest.log
#PBS -e arrayTest.err
#PBS -l nodes=1:ppn=4,mem=4gb,walltime=04:00:00

source /etc/profile.d/modules.sh
module load bioinfo/fastqc/0.11.5

echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"
echo "PBS job array id PBS_ARRAYID value is ${PBS_ARRAYID}"

read1=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/fastp_sample_R1.in)
read2=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/fastp_sample_R2.in)

echo "Read 1 is ${read1}" 
echo "Read 2 is ${read2}"

echo "Try running fastqc with:"
fastqc -t 2 ${read1} ${read2} -o /home/ajan/shared_folder/fastqc_fastp

