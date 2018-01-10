#!/bin/bash
#SBATCH --job-name=SW
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=4
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=7-0:0

##-------------------------------------------------
# Indexação do genoma de referência e alinhamentos
# BWA bwasw --> samtools  #2
# Autor: Washington Candeia de Araujo
# Instituto de Medicina Tropical - IMT - UFRN
# Data: 09 Jan 2018
##--------------------------------------------------
set -euo pipefail

# Variáveis 
SAMPLE1=/scratch/global/wcdaraujo/exome/1_sample/
SAMPLE3=/scratch/global/wcdaraujo/exome/3_sample/
SAMPLE4=/scratch/global/wcdaraujo/exome/4_sample/
SAMPLE5=/scratch/global/wcdaraujo/exome/5_sample/

REF=/scratch/global/wcdaraujo/exome/ref/Homo_sapiens.GRCh37.75.dna.fa
ALN=/scratch/global/wcdaraujo/exome/bwasw/aln/
SAM=/scratch/global/wcdaraujo/exome/bwasw/sam/

# Mapear com BWASW
bwa bwasw -t 32 ${REF} ${SAMPLE1}1_R1.fastq.gz ${SAMPLE1}1_R2.fastq.gz -f ${SAM}1_sample.bwasw.sam
bwa bwasw -t 32 ${REF} ${SAMPLE3}3_R1.fastq.gz ${SAMPLE3}3_R2.fastq.gz -f ${SAM}3_sample.bwasw.sam
bwa bwasw -t 32 ${REF} ${SAMPLE4}4_R1.fastq.gz ${SAMPLE4}4_R2.fastq.gz -f ${SAM}4_sample.bwasw.sam
bwa bwasw -t 32 ${REF} ${SAMPLE5}5_R1.fastq.gz ${SAMPLE5}5_R2.fastq.gz -f ${SAM}5_sample.bwasw.sam
