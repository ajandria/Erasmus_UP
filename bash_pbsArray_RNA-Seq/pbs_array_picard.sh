
#!/bin/sh

#PBS -N picard
#PBS -q xeon
#PBS -t 1-29
#PBS -o $PBS_JOBNAME.out
#PBS -e $PBS_JOBNAME.err
#PBS -l nodes=1:ppn=4,mem=8gb,walltime=24:00:00

source /etc/profile.d/modules.sh
module load mbon/miniconda/3

source /cluster/bioinfo/mbon/miniconda/3/etc/profile.d/conda.sh
conda activate picard

echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"
echo "PBS job array id PBS_ARRAYID value is ${PBS_ARRAYID}"

star_dir_n=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/star_dirs.txt)
sample_bam_temp=$(basename ${star_dir_n})
sample_bam=${sample_bam_temp%??????????????????????????????}

echo "Analised .bam file: ${sample_bam}" 

picard CollectRnaSeqMetrics -I ${star_dir_n} -O /home/ajan/shared_folder/picard/${sample_bam}.txt --REF_FLAT /home/ajan/shared_folder/star_hg19_run/chr_primary_assembly/noChr_refFlat.txt --STRAND SECOND_READ_TRANSCRIPTION_STRAND
