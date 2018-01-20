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

## Samtools: Chamada de Variantes, usando o RMDUP (#4)
# É preciso ter feito o mapeamento com BWA mem

# Variáveis ambientais
REF=/scratch/global/wcdaraujo/exome/ref/Homo_sapiens_GRCh38.86.fa 
SAM=/scratch/global/wcdaraujo/exome/sam/
BAM=/scratch/global/wcdaraujo/exome/bam/
MPILEUP=/scratch/global/wcdaraujo/exome/mpileup/

## Chamada de Variantes com Samtools

# 1. Samtools faidx
samtools faidx ${REF} 

# 2. Samtools view para converter SAM a BAM
samtools view -h -b -S ${SAM}1_sample.mem.sam -o ${BAM}1_sample.mem.bam
samtools view -h -b -S ${SAM}3_sample.mem.sam -o ${BAM}3_sample.mem.bam
samtools view -h -b -S ${SAM}4_sample.mem.sam -o ${BAM}4_sample.mem.bam
samtools view -h -b -S ${SAM}5_sample.mem.sam -o ${BAM}5_sample.mem.bam 
# -h --> incluir header (cabeçalho) na saída
# -b --> output no formato BAM
# -S --> input no formato SAM
# -o --> nome do arquivo de saída no formato BAM

# 3. Samtools sort: ordenar os alinhamentos/mapeamentos por posição no genoma
samtools sort ${BAM}1_sample.mem.bam -o ${BAM}1_sample.mem.sorted.bam
samtools sort ${BAM}3_sample.mem.bam -o ${BAM}3_sample.mem.sorted.bam
samtools sort ${BAM}4_sample.mem.bam -o ${BAM}4_sample.mem.sorted.bam
samtools sort ${BAM}5_sample.mem.bam -o ${BAM}5_sample.mem.sorted.bam

# 4. Samtools index
samtools index ${BAM}1_sample.mem.sorted.bam
samtools index ${BAM}3_sample.mem.sorted.bam
samtools index ${BAM}4_sample.mem.sorted.bam
samtools index ${BAM}5_sample.mem.sorted.bam

## 4.1. Visualizar:
# samtools tview 1_sample.mem.sorted.bam.bai ${REF}
# samtools tview 3_sample.mem.sorted.bam.bai ${REF}
# samtools tview 4_sample.mem.sorted.bam.bai ${REF}
# samtools tview 5_sample.mem.sorted.bam.bai ${REF}

# 4. Samtools RMDUP: remoção de duplicatas
samtools rmdup ${BAM}1_sample.mem.sorted.bam.bai ${BAM}1_sample.mem.sorted.rmdup.bam
samtools rmdup ${BAM}3_sample.mem.sorted.bam.bai ${BAM}3_sample.mem.sorted.rmdup.bam
samtools rmdup ${BAM}4_sample.mem.sorted.bam.bai ${BAM}4_sample.mem.sorted.rmdup.bam
samtools rmdup ${BAM}5_sample.mem.sorted.bam.bai ${BAM}5_sample.mem.sorted.rmdup.bam
## Opção: não usar rmdup

# 5. Samtools mpileup: Chamada de Variantes
samtools mpileup -f ${REF} ${BAM}1_sample.mem.sorted.rmdup.bam > ${MPILEUP}1_sample.mpileup
samtools mpileup -f ${REF} ${BAM}3_sample.mem.sorted.rmdup.bam > ${MPILEUP}3_sample.mpileup
samtools mpileup -f ${REF} ${BAM}4_sample.mem.sorted.rmdup.bam > ${MPILEUP}4_sample.mpileup
samtools mpileup -f ${REF} ${BAM}5_sample.mem.sorted.rmdup.bam > ${MPILEUP}5_sample.mpileup

