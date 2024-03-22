#!/bin/bash -l
#SBATCH --job-name=RNASEQC
#SBATCH --output=%j.stdout 
#SBATCH --error=%j.stderr
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50g
#SBATCH --time=01:00:00
#SBATCH --account=project_2004435

rootDir=/scratch/project_2004435/hipsci/rnaseq-pipeline
outputDir=${rootDir}/rnaseqc
ftp_file=${rootDir}/INFO/ftp_list.txt

annotation_file=/scratch/project_2004435/hipsci/genome/gencode.v26.GRCh38.genes.gtf

outDir=${rootDir}/finished

#fastQ1=/scratch/project_2004435/hipsci/fastq/dene/ERR914278_1.fastq.gz
#fastQ2=/scratch/project_2004435/hipsci/fastq/dene/ERR914278_2.fastq.gz
outName=ERR914278

echo "Starting RNA-SeQC..."

export PROJAPPL=/projappl/project_2004435
module load bioconda
source activate hipsci_rnaseqc

mkdir -p ${outputDir}

rnaseqc ${annotation_file} ${outDir}/${outName}_Aligned.sortedByCoord.out.bam ${outputDir} --coverage --verbose --sample=${outName}

echo "Finished RNA-SeQC."

conda deactivate

done

