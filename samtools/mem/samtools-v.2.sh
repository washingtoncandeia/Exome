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
VCF=/scratch/global/wcdaraujo/exome/vcf/

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
samtools sort ${BAM}1_sample.mem.bam -o ${BAM}1_sample.mem.sorted
samtools sort ${BAM}3_sample.mem.bam -o ${BAM}3_sample.mem.sorted
samtools sort ${BAM}4_sample.mem.bam -o ${BAM}4_sample.mem.sorted
samtools sort ${BAM}5_sample.mem.bam -o ${BAM}5_sample.mem.sorted

# 4. Samtools index
samtools index ${BAM}1_sample.mem.sorted.bam 
samtools index ${BAM}3_sample.mem.sorted.bam 
samtools index ${BAM}4_sample.mem.sorted.bam 
samtools index ${BAM}5_sample.mem.sorted.bam 

## 4.1. Visualizar:
# samtools tview 1_sample.mem.sorted.bam ${REF}
# samtools tview 3_sample.mem.sorted.bam ${REF}
# samtools tview 4_sample.mem.sorted.bam ${REF}
# samtools tview 5_sample.mem.sorted.bam ${REF}

# 5. Samtools RMDUP: remoção de duplicatas
samtools rmdup ${BAM}1_sample.mem.sorted.bam ${BAM}1_sample.mem.sorted.rmdup.bam
samtools rmdup ${BAM}3_sample.mem.sorted.bam ${BAM}3_sample.mem.sorted.rmdup.bam
samtools rmdup ${BAM}4_sample.mem.sorted.bam ${BAM}4_sample.mem.sorted.rmdup.bam
samtools rmdup ${BAM}5_sample.mem.sorted.bam ${BAM}5_sample.mem.sorted.rmdup.bam
## Opção: não usar rmdup

# 6. Samtools mpileup: Chamada de Variantes
samtools mpileup -uf ${REF} ${BAM}1_sample.mem.sorted.rmdup.bam | bcftools view -bvcg - > ${MPILEUP}1_variants.raw.bcf
samtools mpileup -uf ${REF} ${BAM}3_sample.mem.sorted.rmdup.bam | bcftools view -bvcg - > ${MPILEUP}3_variants.raw.bcf
samtools mpileup -uf ${REF} ${BAM}4_sample.mem.sorted.rmdup.bam | bcftools view -bvcg - > ${MPILEUP}4_variants.raw.bcf
samtools mpileup -uf ${REF} ${BAM}5_sample.mem.sorted.rmdup.bam | bcftools view -bvcg - > ${MPILEUP}5_variants.raw.bcf

# 7. Bcftools view: converter .bcf a .vcf
bcftools view ${MPILEUP}1_variant.raw.bcf > ${VCF}1_variants.vcf
bcftools view ${MPILEUP}3_variant.raw.bcf > ${VCF}3_variants.vcf
bcftools view ${MPILEUP}4_variant.raw.bcf > ${VCF}4_variants.vcf
bcftools view ${MPILEUP}5_variant.raw.bcf > ${VCF}5_variants.vcf

