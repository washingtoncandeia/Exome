#!/bin/bash
##-------------------------------------------------------------------------------
# Washington C. Araujo
# IMT - UFRN
# Usar antes de iniciar experimento 
# How to create a database for BWA and BWASW
# https://edwards.sdsu.edu/research/how-to-create-a-database-for-bwa-and-bwa-sw/
# Download GRCh37.75 Ensembl:
# ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/
# Download GRCh38.86 Ensembl:
# ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/
# ftp://ftp.ensembl.org/../pub/release-86/fasta/homo_sapiens/dna/
##-------------------------------------------------------------------------------

# 1. Download sequence database
for CHROM in {1..22} X Y MT;
do
   wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.${CHROM}.fa.gz;
done;

# 2. Unzip as sequências
for CHROM in {1..22} X Y MT;
do
   gzip -dvc Homo_sapiens.GRCh37.75.dna.chromosome.${CHROM}.fa.gz >> Homo_sapiens.GRCh37.75.dna.fa; rm Homo_sapiens.GRCh37.75.dna.chromosome.${CHROM}.fa.gz;
done;

# 3. Retirar do banco de dados (referência Homo sapiens NCBI, GRCh38.p7) os Ns,
# que levam à ambiguidades:

cat Homo_sapiens.GRCh37.75.dna.fa | perl -p -e 's/Nn/N/' | perl -p -e 's/^N+//;s/N+$//;s/N{200,}/n>splitn/' > Homo_sapiens.GRCh37.75.dna_split.fa  

# Após isso, usar biopython para corrigir arquivo fasta e levá-lo a ser usado por prinseq

# 4. Biopython corrigir fasta
# python fastaCorrect.py
# Lembrar de mudar os nomes dos arquivos. 

# 5. Prinseq para limpeza
perl prinseq-lite.pl -log -verbose -fasta Homo_sapiens.GRCh37.75.dna_split.fa -min_len 200 -ns_max_p 10 -derep 12345 -out_good Homo_sapiens.GRCh37.75.dna_split_prinseq -seq_id Homo_sapiens.GRCh37.75_ -rm_header -out_bad null

