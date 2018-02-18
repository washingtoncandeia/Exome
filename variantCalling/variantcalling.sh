#!/bin/bash
#SBATCH --job-name=VarCall
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=2
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=12-0:0

##-------------------------------------------------
# Indexação do genoma de referência e alinhamentos
# BWA mem --> samtools --> picard tools
# Economia de espaço em disco.
# Autor: Washington Candeia de Araujo
# Instituto de Medicina Tropical - IMT - UFRN
# Data: 01 FEV 2018
##--------------------------------------------------
set -euo pipefail

# Variaveis 

# Variaveis ambientais
REF=/scratch/global/wcdaraujo/exome/ref/GRCh38.86/Homo_sapiens_GRCh38.86.fa 

SAMPLE1=/scratch/global/wcdaraujo/exome/1_sample/
SAMPLE3=/scratch/global/wcdaraujo/exome/3_sample/
SAMPLE4=/scratch/global/wcdaraujo/exome/4_sample/
SAMPLE5=/scratch/global/wcdaraujo/exome/5_sample/

BAM=/scratch/global/wcdaraujo/exome/bam/
VCF=/scratch/global/wcdaraujo/exome/vcf/
VARIANT=/scratch/global/wcdaraujo/exome/variantCalling/

# 1. Indexar genoma de referência no formato FASTA
bwa index -a bwtsw ${REF}

# 2. Alinhar com BWA mem
bwa mem -t 64 ${REF} ${SAMPLE1}1_R1.fastq.gz ${SAMPLE1}1_R2.fastq.gz | samtools view -h -b -S -F4 - > ${BAM}1_sample.mem.bam
bwa mem -t 64 ${REF} ${SAMPLE3}3_R1.fastq.gz ${SAMPLE3}3_R2.fastq.gz | samtools view -h -b -S -F4 - > ${BAM}3_sample.mem.bam
bwa mem -t 64 ${REF} ${SAMPLE4}4_R1.fastq.gz ${SAMPLE4}4_R2.fastq.gz | samtools view -h -b -S -F4 - > ${BAM}4_sample.mem.bam
bwa mem -t 64 ${REF} ${SAMPLE5}5_R1.fastq.gz ${SAMPLE5}5_R2.fastq.gz | samtools view -h -b -S -F4 - > ${BAM}5_sample.mem.bam
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

## Visualizar:
# samtools tview 1_sample.mem.sorted.bam ${REF}
# samtools tview 3_sample.mem.sorted.bam ${REF}
# samtools tview 4_sample.mem.sorted.bam ${REF}
# samtools tview 5_sample.mem.sorted.bam ${REF}


# 5. Samtools mpileup e VarScan: Chamada de Variantes com todas as amostras
samtools mpileup -f ${REF} \
                    ${BAM}1_sample.mem.sorted.bam \
                    ${BAM}3_sample.mem.sorted.bam \
                    ${BAM}4_sample.mem.sorted.bam \
                    ${BAM}5_sample.mem.sorted.bam \
                    | java -jar VarScan.v2.4.3.jar mpileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}allSamples.vcf
                   
                   
# 6. snpEff - Preditor
java -jar /home/wcdaraujo/snpEff/snpEff.jar ann -c /home/wcdaraujo/snpEff/snpEff.config GRCh38.86 ${VCF}allSamples.vcf > ${VARIANT}allSamples.eff.vcf


# 7. Samtools mpileup e VarScan: Chamada de Variantes em amostras individuais
samtools mpileup -f ${REF} ${BAM}1_sample.mem.sorted.bam | java -jar VarScan.v2.4.3.jar pileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}1_sample.vcf
samtools mpileup -f ${REF} ${BAM}3_sample.mem.sorted.bam | java -jar VarScan.v2.4.3.jar pileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}3_sample.vcf
samtools mpileup -f ${REF} ${BAM}4_sample.mem.sorted.bam | java -jar VarScan.v2.4.3.jar pileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}4_sample.vcf
samtools mpileup -f ${REF} ${BAM}5_sample.mem.sorted.bam | java -jar VarScan.v2.4.3.jar pileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}5_sample.vcf
 
 
# 8. snpEff - Preditor
java -jar /home/wcdaraujo/snpEff/snpEff.jar ann -c /home/wcdaraujo/snpEff/snpEff.config GRCh38.86 ${VCF}1_sample.vcf > ${VARIANT}1_sample.eff.vcf
java -jar /home/wcdaraujo/snpEff/snpEff.jar ann -c /home/wcdaraujo/snpEff/snpEff.config GRCh38.86 ${VCF}3_sample.vcf > ${VARIANT}3_sample.eff.vcf
java -jar /home/wcdaraujo/snpEff/snpEff.jar ann -c /home/wcdaraujo/snpEff/snpEff.config GRCh38.86 ${VCF}4_sample.vcf > ${VARIANT}4_sample.eff.vcf
java -jar /home/wcdaraujo/snpEff/snpEff.jar ann -c /home/wcdaraujo/snpEff/snpEff.config GRCh38.86 ${VCF}5_sample.vcf > ${VARIANT}5_sample.eff.vcf

