#!/bin/bash
#SBATCH --job-name=ALN
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=4
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=7-0:0

# Utilizar este arquivo como base para os três scripts
##-------------------------------------------------
# Indexação do genoma de referência e alinhamentos
# BWA index --> aln --> sampe (#1)
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

# 1. Indexar genoma de referência no formato FASTA
bwa index -a bwtsw ${REF}

# 2a. Alinhamento em dois passos: aln e sampe
bwa aln ${REF} ${SAMPLE1}1_R1.fastq.gz > ${ALN}1_R1_aln.sai
bwa aln ${REF} ${SAMPLE1}1_R2.fastq.gz > ${ALN}1_R2_aln.sai

bwa aln ${REF} ${SAMPLE3}3_R1.fastq.gz > ${ALN}3_R1_aln.sai
bwa aln ${REF} ${SAMPLE3}3_R2.fastq.gz > ${ALN}3_R2_aln.sai

bwa aln ${REF} ${SAMPLE4}4_R1.fastq.gz > ${ALN}4_R1_aln.sai
bwa aln ${REF} ${SAMPLE4}4_R2.fastq.gz > ${ALN}4_R2_aln.sai

bwa aln ${REF} ${SAMPLE5}5_R1.fastq.gz > ${ALN}5_R1_aln.sai
bwa aln ${REF} ${SAMPLE5}5_R2.fastq.gz > ${ALN}5_R2_aln.sai

# 2b. Sampe para gerar alinhamentos paired end no formato SAM
bwa sampe ${REF} ${ALN}1_R1_aln.sai ${ALN}1_R2_aln.sai ${SAMPLE1}1_R1.fastq.gz ${SAMPLE1}1_R2.fastq.gz > ${SAM}1_sample.bwasw.sam
bwa sampe ${REF} ${ALN}3_R1_aln.sai ${ALN}3_R2_aln.sai ${SAMPLE3}3_R1.fastq.gz ${SAMPLE3}3_R2.fastq.gz > ${SAM}3_sample.bwasw.sam
bwa sampe ${REF} ${ALN}4_R1_aln.sai ${ALN}4_R2_aln.sai ${SAMPLE4}4_R1.fastq.gz ${SAMPLE4}4_R2.fastq.gz > ${SAM}4_sample.bwasw.sam
bwa sampe ${REF} ${ALN}5_R1_aln.sai ${ALN}5_R2_aln.sai ${SAMPLE5}5_R1.fastq.gz ${SAMPLE5}5_R2.fastq.gz > ${SAM}5_sample.bwasw.sam

