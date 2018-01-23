#!/bin/bash
#SBATCH --job-name=MapSamVar
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=4
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=12-0:0

# Variáveis ambientais
REF=/scratch/global/wcdaraujo/exome/ref/GRCh38.86/Homo_sapiens_GRCh38.86.fa 
SAMPLE1=/scratch/global/wcdaraujo/exome/1_sample/
SAMPLE3=/scratch/global/wcdaraujo/exome/3_sample/
SAMPLE4=/scratch/global/wcdaraujo/exome/4_sample/
SAMPLE5=/scratch/global/wcdaraujo/exome/5_sample/
SAM=/scratch/global/wcdaraujo/exome/sam/
BAM=/scratch/global/wcdaraujo/exome/bam/
MPILEUP=/scratch/global/wcdaraujo/exome/mpileup/
VCF=/scratch/global/wcdaraujo/exome/vcf/
VARIANT=/scratch/global/wcdaraujo/exome/variantCalling/

# 1. Uso de bcftools para BCF

samtools mpileup -u -tDP -f ${REF} ${BAM}1_sample.mem.sorted.rmdup.bam | bcftools view -vcg - > variant_sample1.raw.bcf

#2. Trasnformar em VCF
bcftools view variant_sample1.raw.bcf > 1_sample.variant.vcf

