#!/bin/bash
#SBATCH --job-name=Samtools
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=4
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=7-0:0

set -euo pipefail

## Samtools: Chamada de Variantes, usando o RMDUP e sem usar opção Index de Samtools  #4
## Semelhante ao tutorial do curso do Biome

# Variáveis ambientais
REF=/scratch/global/wcdaraujo/exome/bwasw/ref/Homo_sapiens_GRCh38_p7_correct.fa 
SAM=/scratch/global/wcdaraujo/exome/bwasw/sam/
BAM=/scratch/global/wcdaraujo/exome/bwasw/bam/
MPILEUP=/scratch/global/wcdaraujo/exome/bwasw/mpileup/

## Chamada de Variantes com Samtools

# 1. Samtools faidx
samtools faidx ${REF} 

# 2. Samtools view para transformar SAM a BAM
samtools view -b -S ${SAM}1_sample.bwasw.sam -o ${BAM}1_sample.bwasw.bam
samtools view -b -S ${SAM}3_sample.bwasw.sam -o ${BAM}3_sample.bwasw.bam
samtools view -b -S ${SAM}4_sample.bwasw.sam -o ${BAM}4_sample.bwasw.bam
samtools view -b -S ${SAM}5_sample.bwasw.sam -o ${BAM}5_sample.bwasw.bam 
# -b --> output no formato BAM
# -S --> input no formato SAM
# -o --> nome do arquivo de saída no formato BAM

# 3. Samtools sort: ordenar os alinhamentos/mapeamentos nos arquivos .BAM
samtools sort ${BAM}1_sample.bwasw.bam -o ${BAM}1_sample.bwasw.sorted.bam
samtools sort ${BAM}3_sample.bwasw.bam -o ${BAM}3_sample.bwasw.sorted.bam
samtools sort ${BAM}4_sample.bwasw.bam -o ${BAM}4_sample.bwasw.sorted.bam
samtools sort ${BAM}5_sample.bwasw.bam -o ${BAM}5_sample.bwasw.sorted.bam

# 4. Samtools index
# samtools index ${BAM}1_sample.bwasw.sorted.bam ${BAM}1_sample.bwasw.sorted.bam.bai
# samtools index ${BAM}3_sample.bwasw.sorted.bam ${BAM}3_sample.bwasw.sorted.bam.bai
# samtools index ${BAM}4_sample.bwasw.sorted.bam ${BAM}4_sample.bwasw.sorted.bam.bai
# samtools index ${BAM}5_sample.bwasw.sorted.bam ${BAM}5_sample.bwasw.sorted.bam.bai

# 4. Samtools RMDUP: remoção de duplicatas
samtools rmdup ${BAM}1_sample.bwasw.sorted.bam ${BAM}1_sample.bwasw.sorted.rmdup.bam
samtools rmdup ${BAM}3_sample.bwasw.sorted.bam ${BAM}3_sample.bwasw.sorted.rmdup.bam
samtools rmdup ${BAM}4_sample.bwasw.sorted.bam ${BAM}4_sample.bwasw.sorted.rmdup.bam
samtools rmdup ${BAM}5_sample.bwasw.sorted.bam ${BAM}5_sample.bwasw.sorted.rmdup.bam

# 5. Samtools mpileup
samtools mpileup -f ${REF} ${BAM}1_sample.bwasw.sorted.rmdup.bam > ${MPILEUP}1_sample.mpileup
samtools mpileup -f ${REF} ${BAM}3_sample.bwasw.sorted.rmdup.bam > ${MPILEUP}3_sample.mpileup
samtools mpileup -f ${REF} ${BAM}4_sample.bwasw.sorted.rmdup.bam > ${MPILEUP}4_sample.mpileup
samtools mpileup -f ${REF} ${BAM}5_sample.bwasw.sorted.rmdup.bam > ${MPILEUP}5_sample.mpileup

