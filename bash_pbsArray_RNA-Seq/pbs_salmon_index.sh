
#!/bin/sh

#PBS -N index
#PBS -q xeon
#PBS -o $PBS_JOBNAME.out
#PBS -e $PBS_JOBNAME.err
#PBS -l nodes=1:ppn=16,mem=16gb,walltime=01:00:00

source /etc/profile.d/modules.sh
module load mbon/miniconda/3

source /cluster/bioinfo/mbon/miniconda/3/etc/profile.d/conda.sh
conda activate salmon

salmon index -t /home/ajan/shared_folder/salmon_index/Homo_sapiens.GRCh37.75.cdna.all.fa -i /home/ajan/shared_folder/salmon_index/index -p 16
