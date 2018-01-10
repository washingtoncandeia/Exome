#!/bin/env python

# Programa para corrigir arquivo fasta.
# Importante ao usar prinseq

from Bio import SeqIO
SeqIO.convert("Homo_sapiens.GRCh37.75.dna_split.fa", "fasta", "Homo_sapiens.GRCh37.75.dna_split.fa", "fasta")

