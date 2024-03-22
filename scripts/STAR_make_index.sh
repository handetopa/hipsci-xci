#!/bin/bash -l
#SBATCH --job-name=STAR_make_index
#SBATCH --output=STAR_make_index.stdout
#SBATCH --error=STAR_make_index.stderr
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128g
#SBATCH --time=24:00:00
#SBATCH --account=project_2004435

# calculate indexes. You don't need to recalculte the indexes if they already exist.

#module load biokit
#mkdir ../genome-index38
#STAR --runMode genomeGenerate --genomeDir ../genome-index38 --genomeFastaFiles ../genome/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta  --sjdbGTFfile ../genome/gencode.v26.annotation.gtf --runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 75

export PROJAPPL=/projappl/project_2004435
module load bioconda
source activate hipsci_star

mkdir ../genome-index38_star_2.6.0c
STAR --runMode genomeGenerate --genomeDir ../genome-index38_star_2.6.0c --genomeFastaFiles ../genome/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta  --sjdbGTFfile ../genome/gencode.v26.annotation.gtf --runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 75

conda deactivate

