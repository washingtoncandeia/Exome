#!/bin/bash
#SBATCH --job-name=ALN
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=1-0:0

java -jar picard.jar AddOrReplaceReadGroups \
      I=1_sample.mem.bam \
      O=1_sample.mem.picard.bam \
      RGID=1 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=sample1


java -jar picard.jar AddOrReplaceReadGroups \
      I=3_sample.mem.bam \
      O=3_sample.mem.picard.bam \
      RGID=3 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=sample3


java -jar picard.jar AddOrReplaceReadGroups \
      I=4_sample.mem.bam \
      O=4_sample.mem.picard.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=sample4

java -jar picard.jar picard AddOrReplaceReadGroups \
      I=5_sample.mem.bam \
      O=5_sample.mem.picard.bam \
      RGID=5 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=sample5


