#!/bin/bash -l
#SBATCH --job-name=ase2
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --partition=small
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50g
#SBATCH --time=24:00:00
#SBATCH --account=project_2004435

i=$1

u=updated_

module load gatk/4.2.6.6     
#module load picard/2.27.4
#module load biokit
#module load samtools/1.16.1

rootDir=/scratch/project_2004435/hipsci_ase
commonFiles=${rootDir}/common_vcf_bam.txt
refFile=${rootDir}/genome/hs37d5.fa

vcfDir=${rootDir}/vcf
bamDir=${rootDir}/bam
aseDir=${rootDir}/ase

sample=$(eval sed -n '${i}p' ${commonFiles} | cut -d' ' -f1)
bamFile=${bamDir}/$(eval sed -n '${i}p' ${commonFiles} | cut -d' ' -f3)
vcfFile=${vcfDir}/$(eval sed -n '${i}p' ${commonFiles} | cut -d' ' -f4)
updated_vcfFile=${vcfDir}/${u}$(eval sed -n '${i}p' ${commonFiles} | cut -d' ' -f4)

aseFile=${aseDir}/${sample}.table
mkdir -p ${aseDir}

echo ${i}
echo ${sample}
echo ${bamFile}
echo ${vcfFile}
echo ${aseFile}

echo "SelectVariants started."
gatk SelectVariants \
     -R ${refFile} \
     -V ${vcfFile} \
     --restrict-alleles-to BIALLELIC \
     --select-type-to-include SNP \
     --exclude-non-variants \
     --exclude-filtered \
     -O ${filtered_vcfFile}
echo "SelectVariants finished."


echo "ASEReadCounter started."

gatk ASEReadCounter \
 --min-mapping-quality 10 \
 --min-base-quality 2 \
 -R ${refFile} \
 -I ${bamFile} \
 -V ${filtered_vcfFile} \
 -O ${aseFile}

echo "ASEReadCounter finished."
