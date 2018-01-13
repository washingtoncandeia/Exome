#!/bin/bash
##-------------------------------------------------------------------------------
# Washington C. Araujo
# IMT - UFRN
# Data: 12 jan 2018
# Usar antes de iniciar experimento 
# How to create a database for bwasw and mem
# https://edwards.sdsu.edu/research/how-to-create-a-database-for-bwa-and-bwa-sw/

# Download GRCh37.75 Ensembl:
# ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/

# Download GRCh38.86 Ensembl:
# ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/
# ftp://ftp.ensembl.org/../pub/release-86/fasta/homo_sapiens/dna/
##-------------------------------------------------------------------------------

# 1. Download sequence database Ensembl

for CHROM in {1..22} X Y MT;
do
   wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${CHROM}.fa.gz;
done;

# 2. Unzip as sequências
for CHROM in {1..22} X Y MT;
do
   gzip -dvc Homo_sapiens.GRCh38.dna.chromosome.${CHROM}.fa.gz; rm Homo_sapiens.GRCh38.dna.chromosome.${CHROM}.fa.gz;
done;

# 3. Concatenar com cat
for CHROM in {1.22} X Y MT; do cat Homo_sapiens.GRCh38.dna.chromosome.${CHROM}.fa >> Homo_sapiens.GRCh38.86.fa; done;

# 4. Retirar do banco de dados (referência Homo sapiens NCBI, GRCh38.p7) os Ns,
# que levam à ambiguidades:

cat Homo_sapiens.GRCh38.86.fa | perl -p -e 's/Nn/N/' | perl -p -e 's/^N+//;s/N+$//;s/N{200,}/n>splitn/' > Homo_sapiens.GRCh38.86.split.fa  

# Após isso, usar biopython para corrigir arquivo fasta e levá-lo a ser usado por prinseq

# 5. Biopython corrigir fasta
# python fastaCorrect.py
# Lembrar de mudar os nomes dos arquivos. 

# 5. Prinseq para limpeza
perl prinseq-lite.pl -log -verbose -fasta Homo_sapiens.GRCh38.86.split.fa \ 
    -min_len 200 \ 
    -ns_max_p 10 \ 
    -derep 12345 \ 
    -out_good Homo_sapiens.GRCh38.86.split_prinseq \ 
    -seq_id Homo_sapiens.GRCh38.86. \ 
    -rm_header \ 
    -out_bad null
    
## Opções de prinseq :
# -log --> Arquivo de log para acompanhar parâmetros, erros, etc. O nome do arquivo de log é opcional.
# -verbose --> Envia à tela mensagens de informação e status durante o processamento.
# -fasta --> Arquivo de entrada, contendo a sequência a se limpar.
# -min_len --> Filtra as sequências mais curtas do que o especificado em -min_len
# -ns_max_p --> Filtra sequências com valor de porcentagem de Ns maior que especificado em -ns_max_p.
# -derep 12345 -->  Tipos de replicatas a filtrar: 
#                   1 (exact duplicate), 2 (5' duplicate), 3 (3' duplicate), 
#                   4 (reverse complement exact duplicate), 5 (reverse complement 5'/3' duplicate).
# -out_good --> Para alterar o nome do arquivo de saída e a localização, especifique o nome do arquivo usando esta opção.
#               A extensão do arquivo será adicionada automaticamente (quer .fasta, .qual ou .fastq).            
# -seq_id --> Renomeia o identificador da sequência. Um contador é adicionado a cada identificador para garantir que este seja único.
# -rm_header --> Remove o cabeçalho (header) da sequência. Isso inclui tudo após o identificador de sequência (que é mantido inalterado).
# -out_bad null --> Use "-out_bad null" para evitar que o programa gere o (s) arquivo (s) de saída para dados que não passem por nenhum filtro. 


## Código original:
#perl prinseq-lite.pl -log -verbose -fasta hs_ref_GRCh37_p2_split.fa 
#                    -min_len 200 \ 
#                    -ns_max_p 10 \ 
#                    -derep 12345 \ 
#                    -out_good hs_ref_GRCh37_p2_split_prinseq \ 
#                    -seq_id hs_ref_GRCh37_p2_ \ 
#                    -rm_header \ 
#                    -out_bad null
                    
