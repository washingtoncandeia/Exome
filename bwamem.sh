#!/bin/bash
#SBATCH --job-name=MEM
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=4
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=7-0:0

##-------------------------------------------------
# Indexação do genoma de referência e alinhamentos
# BWA mem --> samtools  #3
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

# Mapear com MEM
bwa mem -t 32 ${REF} ${SAMPLE1}1_R1.fastq.gz ${SAMPLE1}1_R2.fastq.gz -f ${SAM}1_sample.mem.sam
bwa mem -t 32 ${REF} ${SAMPLE3}3_R1.fastq.gz ${SAMPLE3}3_R2.fastq.gz -f ${SAM}3_sample.mem.sam
bwa mem -t 32 ${REF} ${SAMPLE4}4_R1.fastq.gz ${SAMPLE4}4_R2.fastq.gz -f ${SAM}4_sample.mem.sam
bwa mem -t 32 ${REF} ${SAMPLE5}5_R1.fastq.gz ${SAMPLE5}5_R2.fastq.gz -f ${SAM}5_sample.mem.sam

# BWA-MEM uses quite a different strategy. 
# It detects long exact matches between the read and the reference, 
# and then chains them into local or global alignments, based on what is more appropriate in that specific case. 
# Such an automatic local-global switching can be very powerful and BWA-MEM works well 
# with various types of data (short reads, long reads, low error rates, high error rates, etc.).
# (https://bioinformatics.stackexchange.com/questions/15/difference-between-bwa-backtrack-and-bwa-mem)
