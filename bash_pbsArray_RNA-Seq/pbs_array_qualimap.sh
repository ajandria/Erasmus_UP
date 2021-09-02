
#!/bin/sh

#PBS -N qualimap
#PBS -q xeon
#PBS -t 1-29
#PBS -o $PBS_JOBNAME.out
#PBS -e $PBS_JOBNAME.err
#PBS -l nodes=1:ppn=4,mem=16gb,walltime=24:00:00

source /etc/profile.d/modules.sh
module load mbon/miniconda/3

source /cluster/bioinfo/mbon/miniconda/3/etc/profile.d/conda.sh
conda activate qualimap

echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"
echo "PBS job array id PBS_ARRAYID value is ${PBS_ARRAYID}"

star_dir_n=$(sed -n "${PBS_ARRAYID}p" /home/ajan/shared_folder/bams_name_sorted.txt)
sample_bam=$(basename ${star_dir_n})

echo "Analised .bam file: ${sample_bam}" 
cd /home/ajan/shared_folder/qualimap
mkdir ${sample_bam}
qualimap rnaseq -bam ${star_dir_n} -gtf /home/ajan/shared_folder/star_hg19_run/chr_primary_assembly/Homo_sapiens.GRCh37.87.chr.gtf -pe -p strand-specific-reverse -s -outformat PDF:HTML -outdir /home/ajan/shared_folder/qualimap/${sample_bam} --java-mem-size=16G

