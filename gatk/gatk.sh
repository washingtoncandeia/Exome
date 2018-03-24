#!/bin/bash

##-------------------------------
# Análise de Exoma - GATK
# José Eduardo Kroll
# Uso em análise de exoma
##-------------------------------


java -Djava.io.tmpdir=./TMP  -Xmx200G -jar /data2/jkroll/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar \
     -T BaseRecalibrator  \
     -nct 40 \
     -R /data2/jkroll/BWA_HG38/Homo_sapiens_assembly38.fasta \
     -I ARQUIVO.bam -o saida.bsqr --knownSites /data2/jkroll/KNOWNSITES/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
     --knownSites /data2/jkroll/KNOWNSITES/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

java -Djava.io.tmpdir=./TMP -Xmx200G -jar /data2/jkroll/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar \
     -T PrintReads \
     -nct 40 \
     -R /data2/jkroll/BWA_HG38/Homo_sapiens_assembly38.fasta \
     -I INPUT.bam \
     -o SAIDA.bam \
     -BQSR input.bsqr

java -Djava.io.tmpdir=./TMP -Xmx150G -jar /data2/jkroll/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -R /data2/jkroll/BWA_HG38/Homo_sapiens_assembly38.fasta \
     -maxReadsInRegionPerSample 10000 \
     --minPruning 100 \
     --dbsnp /data2/jkroll/KNOWNSITES/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
     -nct 40 \
     --genotyping_mode DISCOVERY \
     -I TOHAPLOTYPE.list  \
     -o RESULT.vcf


