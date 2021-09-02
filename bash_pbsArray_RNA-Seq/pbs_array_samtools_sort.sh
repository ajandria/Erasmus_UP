
#!/bin/sh

#PBS -N sort
#PBS -q xeon
#PBS -t 1-29
#PBS -o $PBS_JOBNAME.out
#PBS -e $PBS_JOBNAME.err
#PBS -l nodes=1:ppn=8,mem=16gb,walltime=24:00:00

source /etc/profile.d/modules.sh
module load mbon/samtools/1.9

echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"
echo "PBS job array id PBS_ARRAYID value is ${PBS_ARRAYID}"

star_dir_n=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/star_dirs.txt)
sample_bam_temp=$(basename ${star_dir_n})
sample_bam=${sample_bam_temp%??????????????????????????????}

samtools sort -n -@ 8 -O BAM -o /home/ajan/shared_folder/star_hg19_run/name_sorted_alignments/${sample_bam}.Name.Sorted.bam ${star_dir_n}
