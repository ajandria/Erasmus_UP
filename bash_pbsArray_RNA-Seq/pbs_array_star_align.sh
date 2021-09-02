#!/bin/sh

#PBS -N star_align
#PBS -q xeon
#PBS -t 1-29
#PBS -o star_align.log
#PBS -e star_align.err
#PBS -l nodes=1:ppn=16,mem=8gb,walltime=04:00:00

source /etc/profile.d/modules.sh
module load mbon/miniconda/3

source /cluster/bioinfo/mbon/miniconda/3/etc/profile.d/conda.sh
conda activate star

echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"
echo "PBS job array id PBS_ARRAYID value is ${PBS_ARRAYID}"

read1=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/fastp_sample_R1.in)
read2=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/fastp_sample_R2.in)

read_base=$(basename ${read1})
base_name=${read_base%?????????????????}

echo "Read 1 is ${read1}" 
echo "Read 2 is ${read2}"

echo "Try running fastqc with:"
cd /home/ajan/shared_folder/star_hg19_run/alignments
STAR --runThreadN 16 --runMode alignReads --genomeDir /home/ajan/shared_folder/star_hg19_run/index --readFilesIn ${read1} ${read2} --readFilesCommand zcat --readStrand Reverse --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /home/ajan/shared_folder/star_hg19_run/chr_primary_assembly/Homo_sapiens.GRCh37.87.chr.gtf --sjdbOverhang 50 --outFileNamePrefix ${base_name}_

