
#!/bin/sh

#PBS -N star_index
#PBS -q xeon
#PBS -o $PBS_JOBNAME.out
#PBS -e $PBS_JOBNAME.err
#PBS -l nodes=1:ppn=16,mem=48gb,walltime=10:00:00

source /etc/profile.d/modules.sh
module load mbon/STAR/2.5.3a

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /home/ajan/shared_folder/star_hg19_run/index --genomeFastaFiles /home/ajan/shared_folder/star_hg19_run/chr_primary_assembly/cleaned_Homo_sapiens.GRCh37.dna.primary_assembly.fa
