#!/bin/bash -l                                                                                                                                                                                  
#SBATCH --job-name=STAR_alignment                                                                                                                                                                   
#SBATCH --output=%j.stdout                                                                                                                                                                      
#SBATCH --error=%j.stderr                                                                                                                                                                       
#SBATCH --partition=small                                                                                                                                                                       
#SBATCH --ntasks=1                                                                                                                                                                              
#SBATCH --nodes=1                                                                                                                                                                               
#SBATCH --cpus-per-task=4                                                                                                                                                                       
#SBATCH --mem=50g                                                                                                                                                                               
#SBATCH --time=05:00:00                                                                                                                                                                         
#SBATCH --account=project_2004435                                                                                                                                                               
                                                                                                                                                                                      
i=$1

module load samtools/1.12
module load trimmomatic/0.39
#module load biokit                                                                                                                                                                            \                                                                                                                                                                                              
module load r-env-singularity/4.0.5

rootDir=/scratch/project_2004435/hipsci/rnaseq-pipeline
fastqDir=${rootDir}/fastq
outputDir=${rootDir}/rnaseqc
ftp_file=${rootDir}/INFO/ftp_list.txt

genomeIndexDir=${rootDir}/genome-index38_star_2.6.0c
annotation_file=/scratch/project_2004435/hipsci/genome/gencode.v26.GRCh38.genes.gtf
adapters=/scratch/project_2004435/hipsci/illumina_adapters.fa

tempDir=${rootDir}/tmp
outDir=${rootDir}/finished
QC=${rootDir}/metrics

mkdir -p ${fastqDir}
mkdir -p ${tempDir}
mkdir -p ${outDir}
mkdir -p ${QC}

#for i in `seq ${line1} ${line2}`;                                                                                                                                                              
#do                                                                                                                                                                                             
    echo $i
    outName=$(eval sed -n '${i}p' ${ftp_file} | cut -d' ' -f2)
    ftp1=$(eval sed -n '${i}p' ${ftp_file} | cut -d' ' -f3)
    ftp2=$(eval sed -n '${i}p' ${ftp_file} | cut -d' ' -f4)
    fastQ1=$(eval sed -n '${i}p' ${ftp_file} | cut -d' ' -f3 |  sed 's/.*\///')
    fastQ2=$(eval sed -n '${i}p' ${ftp_file} | cut -d' ' -f4 |  sed 's/.*\///')                                                                                                                                                 
                                                                                                                                                                              
    echo ${outName}
    #cd ${fastqDir}   
    #wget ${ftp1}                                                                                                                                                                               
    #wget ${ftp2}                                                                                                                                                                               
    #cd ${rootDir}                                                                                                                                                                              

    echo "running trimmomatic" # minlength changed from 20 to 15 to 36 (recommended), SLIDING WINDOW 5:20 -- less stringent     
 
    trimmomatic PE \
        -threads 2 -trimlog ${tempDir}/${outName}_trimlogfile \
        ${fastqDir}/${fastQ1} ${fastqDir}/${fastQ2} \
        ${tempDir}/${outName}_unaligned_trimmed-1.fastq.gz ${tempDir}/${outName}_unaligned_discard-1.fastq.gz \
        ${tempDir}/${outName}_unaligned_trimmed-2.fastq.gz ${tempDir}/${outName}_unaligned_discard-2.fastq.gz \
        ILLUMINACLIP:${adapters}:2:30:10 \
        SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36

    echo "Finished trimming."

#    rm -rf ${fastqDir}/${fastQ1}                                                                                                                                                               
#    rm -rf ${fastqDir}/${fastQ2}                                                                                                                                                               

    export PROJAPPL=/projappl/project_2004435
    module load bioconda
    source activate hipsci_star

    echo "running STAR version 2.6.0c"

    STAR \
        --runThreadN $SLURM_CPUS_PER_TASK \
        --genomeDir ${genomeIndexDir} \
        --readFilesIn ${tempDir}/${outName}_unaligned_trimmed-1.fastq.gz ${tempDir}/${outName}_unaligned_trimmed-2.fastq.gz \
        --readFilesCommand zcat \
        --outTmpDir ${tempDir}/tmp_${outName} \
        --outFileNamePrefix ${outDir}/${outName}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within                                                                                                                                             

    echo "Finished STAR alignment."

    conda deactivate

#    rm -rf ${tempDir}/${outName}_unaligned_trimmed-1.fastq.gz                                                                                                                                  
#    rm -rf ${tempDir}/${outName}_unaligned_trimmed-2.fastq.gz                                                                                                                                  
#    rm -rf ${tempDir}/${outName}_unaligned_discard-1.fastq.gz                                                                                                                                  
#    rm -rf ${tempDir}/${outName}_unaligned_discard-2.fastq.gz                                                                                                                                  

#    echo "Files deleted."                                                                                                                                                                      

#done                                                                                                                                                                                           