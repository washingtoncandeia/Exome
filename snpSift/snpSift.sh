#!/bin/bash

##-----------------------------------------
# Análise de Exoma - snpSift
# José Eduardo Kroll
# Washington Candeia
# Filtrar e manipular arquivos de anotacao
##-----------------------------------------

db=/data3/jkroll/dbNSFP3.3a.txt.gz

# temp4.vcf --> arquivo VCF allSamples.vcf
# annot=(\
  #'ExAC_AF'\
  #'ExAC_AFR_AF'\
  #'ExAC_AMR_AF'\
  #'ExAC_EAS_AF'\
  #'ExAC_NFE_AF'\
  #'ExAC_SAS_AF'\
  #'SIFT_pred'\
  #'Polyphen2_HDIV_pred'\
  #'Polyphen2_HVAR_pred'\
  #'LRT_pred'\
  #'MutationTaster_pred'\
  #'c'\
  #'FATHMM_pred'\
  #'PROVEAN_pred'\
  #'MetaSVM_pred'\
  #'MetaLR_pred'\
  #'M-CAP_pred'\
  #'fathmm-MKL_coding_pred'\
  #'SiPhy_29way_logOdds'\
  #'GERP++_RS'\
)
#snpSift dbnsfp -f $( echo ${annot[@]} | sed 's/ /,/g' ) -v -a -db $db tmp4.vcf | gzip > $1.eff_new.gz

snpSift dbnsfp -f ExAC_AF, ExAC_AFR_AF, ExAC_AMR_AF, ExAC_EAS_AF, ExAC_NFE_AF, \
                  ExAC_SAS_AF, SIFT_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, LRT_pred, \
                  MutationTaster_pred, c, FATHMM_pred, PROVEAN_pred, MetaSVM_pred, \
                  MetaLR_pred, M-CAP_pred, fathmm-MKL_coding_pred, SiPhy_29way_logOdds, GERP++_RS \
                  -v -a -db ${db} tmp allSamples.vcf | gzip > allSamples.eff.gz
                  
                  
