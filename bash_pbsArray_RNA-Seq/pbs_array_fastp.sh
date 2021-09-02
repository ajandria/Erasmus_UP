#!/bin/sh

#PBS -N fastp
#PBS -q xeon
#PBS -t 1-29
#PBS -o arrayFASTP.log
#PBS -e arrayFASTP.err
#PBS -l nodes=1:ppn=8,mem=16gb,walltime=04:00:00

source /etc/profile.d/modules.sh
module load mbon/fastp/0.21.0

echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"
echo "PBS job: array id PBS_ARRAYID value is ${PBS_ARRAYID}"

read1=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/samples_R1.in)
read2=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/samples_R2.in)

echo "Read 1 is ${read1}" 
echo "Read 2 is ${read2}"

read_base=$(basename ${read1})
sample_name=${read_base%???????????}
read1_name=$(basename ${read1})
read1_name=${read1_name%?????????}
read2_name=$(basename ${read2})
read2_name=${read2_name%?????????}


echo "Running fastp with:"

echo "fastp -i ${read1} -I ${read2} -o /home/ajan/shared_folder/fastp_fastq/${read1_name}_fastp.fastq.gz -O /home/ajan/shared_folder/fastp_fastq/${read2_name}_fastp.fastq.gz -q 30 -e 25 -n 3 -l 40 -c -x 20 -w 8 -j /home/ajan/shared_folder/fastp_reports/${sample_name}_fastp.json -h /home/ajan/shared_folder/fastp_reports/${sample_name}_fastp.html"

fastp -i ${read1} -I ${read2} -o /home/ajan/shared_folder/fastp_fastq/${read1_name}_fastp.fastq.gz -O /home/ajan/shared_folder/fastp_fastq/${read2_name}_fastp.fastq.gz -q 30 -e 25 -n 3 -l 40 -c -x 20 -w 8 -j /home/ajan/shared_folder/fastp_reports/${sample_name}_fastp.json -h /home/ajan/shared_folder/fastp_reports/${sample_name}_fastp.html
