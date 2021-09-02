#!/bin/sh

#PBS -N htseq
#PBS -q xeon
#PBS -t 2
#PBS -o htseq.log
#PBS -e htseq.err
#PBS -l nodes=1:ppn=1,mem=8gb,walltime=48:00:00

source /etc/profile.d/modules.sh
module load mbon/miniconda/3

source /cluster/bioinfo/mbon/miniconda/3/etc/profile.d/conda.sh
conda activate htseq

echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"
echo "PBS job array id PBS_ARRAYID value is ${PBS_ARRAYID}"

star_dir_n=$(sed -n "${PBS_ARRAYID}p" /home/ajan/star_dirs.txt)
sample_bam=$(basename ${star_dir_n})

echo "Analised .bam file: ${sample_bam}" 

echo "Running htseq with:"

echo htseq-count -s reverse -f bam -m intersection-nonempty ${sample_bam} /home/ajan/shared_folder/GCF_000001405.25_GRCh37.p13_genomic.gtf

htseq-count -s reverse -f bam -m intersection-nonempty -r pos -t gene ${star_dir_n} /home/ajan/shared_folder/star_hg19_run/merged_1toY/clean_sorted_Homo_sapiens.GRCh37.87.gtf > /home/ajan/shared_folder/htseq/${sample_bam}.txt
