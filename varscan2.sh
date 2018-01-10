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

# Procedendo ao uso de VarScan
java -jar VarScan.v2.4.3.jar pileup2snp ${MPILEUP}1_sample.mpileup --output-vcf --strand-filter 0 > ${VARIANT}1_sample.vcf
java -jar VarScan.v2.4.3.jar pileup2snp ${MPILEUP}3_sample.mpileup --output-vcf --strand-filter 0 > ${VARIANT}3_sample.vcf
java -jar VarScan.v2.4.3.jar pileup2snp ${MPILEUP}4_sample.mpileup --output-vcf --strand-filter 0 > ${VARIANT}4_sample.vcf
java -jar VarScan.v2.4.3.jar pileup2snp ${MPILEUP}5_sample.mpileup --output-vcf --strand-filter 0 > ${VARIANT}5_sample.vcf

