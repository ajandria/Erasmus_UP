#!/bin/sh

#PBS -N quant
#PBS -q xeon
#PBS -t 1-29
#PBS -o quant.log
#PBS -e quant.err
#PBS -l nodes=1:ppn=8,mem=16gb,walltime=08:00:00

source /etc/profile.d/modules.sh
module load mbon/miniconda/3

source /cluster/bioinfo/mbon/miniconda/3/etc/profile.d/conda.sh
conda activate salmon

echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"
echo "PBS job array id PBS_ARRAYID value is ${PBS_ARRAYID}"

read1=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/fastp_sample_R1.in)
read2=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/fastp_sample_R2.in)

read_base=$(basename ${read1})
base_name=${read_base%???????????}

echo "Read 1 is ${read1}" 
echo "Read 2 is ${read2}"

echo "Try running fastqc with:"
salmon quant -i /home/ajan/shared_folder/salmon_index/index -l ISR -1 ${read1} -2 ${read2} -p 8 -o /home/ajan/shared_folder/salmon_quant/${base_name} -g /home/ajan/shared_folder/salmon_index/genes.gtf

