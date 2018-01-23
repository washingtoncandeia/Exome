#!/bin/bash
#SBATCH --job-name=VarCall3
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=2
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=12-0:0

##-------------------------------------------------
# Indexação do genoma de referência e alinhamentos
# BWA mem --> samtools  (#3)
# Autor: Washington Candeia de Araujo
# Instituto de Medicina Tropical - IMT - UFRN
# Data: 23 Jan 2018
##--------------------------------------------------
set -euo pipefail

# Variáveis 

# Variaveis ambientais
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

# 1. Indexar genoma de referência no formato FASTA
bwa index -a bwtsw ${REF}

# Mapear com MEM
bwa mem -t 32 ${REF} ${SAMPLE1}1_R1.fastq.gz ${SAMPLE1}1_R2.fastq.gz > ${SAM}1_sample.mem.sam
bwa mem -t 32 ${REF} ${SAMPLE3}3_R1.fastq.gz ${SAMPLE3}3_R2.fastq.gz > ${SAM}3_sample.mem.sam
bwa mem -t 32 ${REF} ${SAMPLE4}4_R1.fastq.gz ${SAMPLE4}4_R2.fastq.gz > ${SAM}4_sample.mem.sam
bwa mem -t 32 ${REF} ${SAMPLE5}5_R1.fastq.gz ${SAMPLE5}5_R2.fastq.gz > ${SAM}5_sample.mem.sam
## Opção
# bwa mem -t 32 ${REF} ${SAMPLE1}1_R1.fastq.gz ${SAMPLE1}1_R2.fastq.gz -f ${SAM}1_sample.mem.sam

# bwa mem -t 32 ${REF} ${SAMPLE1}1_R1.fastq.gz ${SAMPLE1}1_R2.fastq.gz | samtools view -h -b -S -F4 - > ${BAM}1_sample.mem.bam

## Samtools: Chamada de Variantes, usando o RMDUP (#4)
# É preciso ter feito o mapeamento com BWA mem

## Chamada de Variantes com Samtools

# 1. Samtools faidx
samtools faidx ${REF} 

# 2. Samtools view para converter SAM a BAM
samtools view -h -b -S -F4 ${SAM}1_sample.mem.sam -o ${BAM}1_sample.mem.bam
samtools view -h -b -S ${SAM}3_sample.mem.sam -o ${BAM}3_sample.mem.bam
samtools view -h -b -S ${SAM}4_sample.mem.sam -o ${BAM}4_sample.mem.bam
samtools view -h -b -S ${SAM}5_sample.mem.sam -o ${BAM}5_sample.mem.bam 
# -h --> incluir header (cabeçalho) na saída
# -b --> output no formato BAM
# -S --> input no formato SAM
# -o --> nome do arquivo de saída no formato BAM
# -F 4 --> retirar o que não mapeou

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
samtools mpileup -f ${REF} \ 
                    ${BAM}1_sample.mem.sorted.rmdup.bam \ 
                    ${BAM}3_sample.mem.sorted.rmdup.bam \ 
                    ${BAM}4_sample.mem.sorted.rmdup.bam \ 
                    ${BAM}5_sample.mem.sorted.rmdup.bam \ 
                    > ${MPILEUP}allSamples.mpileup 

##-----------------------------------------------------------------------------------------------------
# Etapas 6 e 7, exemplo de uso para economizar espaço:
#samtools mpileup -f ${REF} \ 
#                    ${BAM}1_sample.mem.sorted.rmdup.bam \ 
#                    ${BAM}3_sample.mem.sorted.rmdup.bam \ 
#                    ${BAM}4_sample.mem.sorted.rmdup.bam \ 
#                    ${BAM}5_sample.mem.sorted.rmdup.bam | java -jar VarScan.v2.2.jar mpileup2snp
##-----------------------------------------------------------------------------------------------------

# 7. Procedendo ao uso de VarScan para todas as amostras unidas
java -jar VarScan.v2.4.3.jar mpileup2snp ${MPILEUP}allSamples.mpileup --output-vcf 1 --strand-filter 0
# java -jar VarScan.v2.4.3.jar mpileup2snp ${MPILEUP}allSamples.mpileup --output-vcf 1 --strand-filter 0 > ${VARIANT}allSamples.vcf

