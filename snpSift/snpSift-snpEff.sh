#!/bin/bash

##-----------------------------------------
# Análise de Exoma - snpSift
# Viviane Nogueira
# Washington Candeia
# Filtrar e manipular arquivos de anotacao
##-----------------------------------------



java -jar /home/vbnogueira/snpEff/SnpSift.jar \
     extractFields /scratch/global/vbnogueira/databases/Male_chr22_multi.vcf \
     AC CHROM POS ID REF ALT ANN[*].GENE  > /scratch/global/vbnogueira/databases/Male_chr22.final.vcf


-
for m in {1..21}
do
     /home/vbnogueira/vcftools/vcftools_0.1.13/bin/vcf-subset \
     -c /scratch/global/vbnogueira/databases/male_list.txt \
     /scratch/global/vbnogueira/databases/chr$m.eff.vcf > /scratch/global/vbnogueira/databases/Male_chr$m.vcf &
done

     /home/vbnogueira/vcftools/vcftools_0.1.13/bin/vcf-subset \
     -c /scratch/global/vbnogueira/databases/male_list.txt \
     /scratch/global/vbnogueira/databases/chr22.eff.vcf > /scratch/global/vbnogueira/databases/Male_chr22.vcf
wait


