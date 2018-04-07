#!/bin/bash
#SBATCH --job-name=mpileup
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=7-0:0

##-------------------------------------------------
# Indexação do genoma de referência e alinhamentos
# BWA mem --> samtools --> .BAM --> .mpileup
# Autor: Washington Candeia de Araujo
# Instituto de Medicina Tropical - IMT - UFRN
# Data: 06 ABRIL 2018
##--------------------------------------------------
set -euo pipefail

# Variaveis 

REF=/home/wcdaraujo/ref/Homo_sapiens.GRCh38.dna.fa

SAMPLE1_R1=/home/wcdaraujo/samples/1sample/1_AL201_20171106000_S1_L006_R1_001.fastq.gz
SAMPLE1_R2=/home/wcdaraujo/samples/1sample/1_AL201_20171106000_S1_L006_R2_001.fastq.gz

SAMPLE3_R1=/home/wcdaraujo/samples/3sample/3_AL208_20171106000_S2_L006_R1_001.fastq.gz
SAMPLE3_R2=/home/wcdaraujo/samples/3sample/3_AL208_20171106000_S2_L006_R2_001.fastq.gz

SAMPLE4_R1=/home/wcdaraujo/samples/4sample/4_SA002_20171106000_S3_L006_R1_001.fastq.gz
SAMPLE4_R2=/home/wcdaraujo/samples/4sample/4_SA002_20171106000_S3_L006_R2_001.fastq.gz

SAMPLE5_R1=/home/wcdaraujo/samples/5sample/5_SA006_20171106000_S4_L006_R1_001.fastq.gz
SAMPLE5_R2=/home/wcdaraujo/samples/5sample/5_SA006_20171106000_S4_L006_R2_001.fastq.gz

BAM=/home/wcdaraujo/bam/
MPILEUP=/home/wcdaraujo/mpileup/
VCF=/home/wcdaraujo/vcf/

# 1. Indexar genoma de referência no formato FASTA
bwa index -a bwtsw ${REF}

# 2. Alinhar com BWA mem
bwa mem -t 32 ${REF} ${SAMPLE1_R1} ${SAMPLE1_R2} | samtools view -h -b -S -F4 - > ${BAM}1_sample.mem.bam
bwa mem -t 32 ${REF} ${SAMPLE3_R1} ${SAMPLE3_R2} | samtools view -h -b -S -F4 - > ${BAM}3_sample.mem.bam
bwa mem -t 32 ${REF} ${SAMPLE4_R1} ${SAMPLE4_R2} | samtools view -h -b -S -F4 - > ${BAM}4_sample.mem.bam
bwa mem -t 32 ${REF} ${SAMPLE5_R1} ${SAMPLE5_R2} | samtools view -h -b -S -F4 - > ${BAM}5_sample.mem.bam
## Samtools
# -h --> incluir header na saída
# -b --> output no formato BAM
# -S --> input no formato SAM
# -o --> nome do arquivo de saída no formato BAM
# -F4 --> retirar o que não foi mapeado

## Chamada de Variantes com Samtools ##

# 1. Samtools faidx
samtools faidx ${REF} 


#2. Picard tools: AddOrReplaceReadGroups
java -jar picard.jar AddOrReplaceReadGroups \
      I=${BAM}1_sample.mem.bam \
      O=${BAM}1_sample.mem.picard.bam \
      RGID="Sample1" \
      RGLB="lib1" \
      RGPL="illumina" \
      RGPU="unit1" \
      RGSM="sample1"
      

java -jar picard.jar AddOrReplaceReadGroups \
      I=${BAM}3_sample.mem.bam \
      O=${BAM}3_sample.mem.picard.bam \
      RGID="Sample3" \
      RGLB="lib3" \
      RGPL="illumina" \
      RGPU="unit3" \
      RGSM="sample3"
      

java -jar picard.jar AddOrReplaceReadGroups \
      I=${BAM}4_sample.mem.bam \
      O=${BAM}4_sample.mem.picard.bam \
      RGID="Sample4" \
      RGLB="lib4" \
      RGPL="illumina" \
      RGPU="unit4" \
      RGSM="sample4"
      

java -jar picard.jar AddOrReplaceReadGroups \
      I=${BAM}5_sample.mem.bam \
      O=${BAM}5_sample.mem.picard.bam \
      RGID="Sample5" \
      RGLB="lib5" \
      RGPL="illumina" \
      RGPU="unit5" \
      RGSM="sample5"
      

# 3. Samtools sort: ordenar os alinhamentos por posição no genoma
samtools sort ${BAM}1_sample.mem.picard.bam -o ${BAM}1_sample.mem.sorted.bam
samtools sort ${BAM}3_sample.mem.picard.bam -o ${BAM}3_sample.mem.sorted.bam
samtools sort ${BAM}4_sample.mem.picard.bam -o ${BAM}4_sample.mem.sorted.bam
samtools sort ${BAM}5_sample.mem.picard.bam -o ${BAM}5_sample.mem.sorted.bam

# 4. Samtools index
samtools index ${BAM}1_sample.mem.sorted.bam 
samtools index ${BAM}3_sample.mem.sorted.bam 
samtools index ${BAM}4_sample.mem.sorted.bam 
samtools index ${BAM}5_sample.mem.sorted.bam 

# 5. Samtools mpileup com todas as amostras
samtools mpileup -f ${REF} \
                    ${BAM}1_sample.mem.sorted.bam \
                    ${BAM}3_sample.mem.sorted.bam \
                    ${BAM}4_sample.mem.sorted.bam \
                    ${BAM}5_sample.mem.sorted.bam \
                    > allSamples.mpileup
                   
# 6. Samtools mpileup: amostras individuais
samtools mpileup -f ${REF} ${BAM}1_sample.mem.sorted.bam > ${MPILEUP}1_sample.mpileup
samtools mpileup -f ${REF} ${BAM}3_sample.mem.sorted.bam > ${MPILEUP}3_sample.mpileup
samtools mpileup -f ${REF} ${BAM}4_sample.mem.sorted.bam > ${MPILEUP}4_sample.mpileup
samtools mpileup -f ${REF} ${BAM}5_sample.mem.sorted.bam > ${MPILEUP}5_sample.mpileup

# 7. Samtools mpileup e VarScan: Chamada de Variantes com todas as amostras
samtools mpileup -f ${REF} \
                    ${BAM}1_sample.mem.sorted.bam \
                    ${BAM}3_sample.mem.sorted.bam \
                    ${BAM}4_sample.mem.sorted.bam \
                    ${BAM}5_sample.mem.sorted.bam \
                    | java -jar VarScan.v2.4.3.jar mpileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}allSamples.vcf
                   

