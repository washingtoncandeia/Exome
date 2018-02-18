#!/bin/bash
#SBATCH --job-name=Var
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=4
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --time=7-0:0

set -euo pipefail

VCF=/scratch/global/wcdaraujo/exome/vcf/
VARIANT=/scratch/global/wcdaraujo/exome/bwasw/variantCalling/

# snpEff - Preditor (Todas as amostras)
java -jar /home/wcdaraujo/snpEff/snpEff.jar ann -c /home/wcdaraujo/snpEff/snpEff.config GRCh38.86 ${VCF}allSamples.vcf > ${VARIANT}allSamples.eff.vcf

# snpEff - Preditor (Cada uma das amostras)
java -jar /home/wcdaraujo/snpEff/snpEff.jar ann -c /home/wcdaraujo/snpEff/snpEff.config GRCh38.86 ${VCF}1_sample.vcf > ${VARIANT}1_sample.eff.vcf
java -jar /home/wcdaraujo/snpEff/snpEff.jar ann -c /home/wcdaraujo/snpEff/snpEff.config GRCh38.86 ${VCF}3_sample.vcf > ${VARIANT}3_sample.eff.vcf
java -jar /home/wcdaraujo/snpEff/snpEff.jar ann -c /home/wcdaraujo/snpEff/snpEff.config GRCh38.86 ${VCF}4_sample.vcf > ${VARIANT}4_sample.eff.vcf
java -jar /home/wcdaraujo/snpEff/snpEff.jar ann -c /home/wcdaraujo/snpEff/snpEff.config GRCh38.86 ${VCF}5_sample.vcf > ${VARIANT}5_sample.eff.vcf

