
#!/bin/sh

#PBS -N featureCount
#PBS -q xeon
#PBS -o $PBS_JOBNAME.out
#PBS -e $PBS_JOBNAME.err
#PBS -l nodes=1:ppn=8,mem=16gb,walltime=24:00:00

source /etc/profile.d/modules.sh
module load mbon/miniconda/3

source /cluster/bioinfo/mbon/miniconda/3/etc/profile.d/conda.sh
conda activate subread

featureCounts -t exon -s 2 -T 16 --verbose -p -a /home/ajan/shared_folder/star_hg19_run/chr_primary_assembly/Homo_sapiens.GRCh37.87.chr.gtf  -o /home/ajan/shared_folder/featureCounts/exon_featureCounts_matrix.txt /home/ajan/shared_folder/star_hg19_run/alignments/*.bam
