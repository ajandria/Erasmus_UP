
#!/bin/sh

#PBS -N samtools
#PBS -q xeon
#PBS -t 1-29
#PBS -o $PBS_JOBNAME.out
#PBS -e $PBS_JOBNAME.err
#PBS -l nodes=1:ppn=4,mem=8gb,walltime=24:00:00

source /etc/profile.d/modules.sh
module load mbon/samtools/1.9 

echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"
echo "PBS job array id PBS_ARRAYID value is ${PBS_ARRAYID}"

star_dir_n=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/star_dirs.txt)
sample_bam_temp=$(basename ${star_dir_n})
sample_bam=${sample_bam_temp%??????????????????????????????}

echo "Analised .bam file: ${sample_bam}" 

echo "Running samtools with:"

echo "samtools index ${star_dir_n}"

echo "samtools flagstat ${star_dir_n} -@ 4 > /home/ajan/shared_folder/samtools/${sample_bam}.txt"

samtools index ${star_dir_n}

samtools flagstat ${star_dir_n} -@ 4 > /home/ajan/shared_folder/samtools/${sample_bam}.txt

