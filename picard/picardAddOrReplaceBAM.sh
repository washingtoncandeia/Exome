#!/bin/bash
#SBATCH --job-name=AddOrRepBam
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=2
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=1-0:0

##-------------------------------------------------
# Indexação do genoma de referência e alinhamentos
# Picard Tools: AddOrReplaceReadGroups BAM
# Autor: Washington Candeia de Araujo
# Instituto de Medicina Tropical - IMT - UFRN
# Data: 01 FEV 2018
##--------------------------------------------------
set -euo pipefail

# Variaveis ambientais
BAM=/scratch/global/wcdaraujo/exome/bam/


# Picard tools: AddOrReplaceReadGroups
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

