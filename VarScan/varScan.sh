#!/bin/bash
#SBATCH --job-name=varscan
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
##SBATCH --nodes=2
##SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=7-0:0

set -euo pipefail

MPILEUP=/scratch/global/wcdaraujo/exome/bwasw/mpileup/
VARIANT=/scratch/global/wcdaraujo/exome/bwasw/variantCalling/

## Uso de VarScan sob diferentes formas

# 1. VarScan: Chamada de Variantes em amostras individuais
java -jar VarScan.v2.4.3.jar pileup2snp ${MPILEUP}1_sample.mpileup --output-vcf 1 --strand-filter 0 > ${VCF}1_sample.vcf
java -jar VarScan.v2.4.3.jar pileup2snp ${MPILEUP}3_sample.mpileup --output-vcf 1 --strand-filter 0 > ${VCF}3_sample.vcf
java -jar VarScan.v2.4.3.jar pileup2snp ${MPILEUP}4_sample.mpileup --output-vcf 1 --strand-filter 0 > ${VCF}4_sample.vcf
java -jar VarScan.v2.4.3.jar pileup2snp ${MPILEUP}5_sample.mpileup --output-vcf 1 --strand-filter 0 > ${VCF}5_sample.vcf


# 2. Samtools mpileup e VarScan: Chamada de Variantes em amostras individuais (em pipes)
samtools mpileup -f ${REF} ${BAM}1_sample.mem.sorted.bam | java -jar VarScan.v2.4.3.jar pileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}1_sample.vcf
samtools mpileup -f ${REF} ${BAM}3_sample.mem.sorted.bam | java -jar VarScan.v2.4.3.jar pileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}3_sample.vcf
samtools mpileup -f ${REF} ${BAM}4_sample.mem.sorted.bam | java -jar VarScan.v2.4.3.jar pileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}4_sample.vcf
samtools mpileup -f ${REF} ${BAM}5_sample.mem.sorted.bam | java -jar VarScan.v2.4.3.jar pileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}5_sample.vcf


# 3. Samtools mpileup e VarScan: Chamada de Variantes com todas as amostras
java -jar VarScan.v2.4.3.jar mpileup2snp ${MPILEUP}allSamples.mpileup --output-vcf 1 --strand-filter 0 > ${VCF}allSamples.vcf


# 4. Samtools mpileup e VarScan: Chamada de Variantes com todas as amostras (em pipe)
samtools mpileup -f ${REF} \
                    ${BAM}1_sample.mem.sorted.bam \
                    ${BAM}3_sample.mem.sorted.bam \
                    ${BAM}4_sample.mem.sorted.bam \
                    ${BAM}5_sample.mem.sorted.bam \
                    | java -jar VarScan.v2.4.3.jar mpileup2snp --output-vcf 1 --strand-filter 0 > ${VCF}allSamples.vcf
                    

