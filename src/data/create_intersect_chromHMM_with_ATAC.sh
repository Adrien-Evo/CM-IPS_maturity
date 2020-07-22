bedtools intersect -a ../external/TSS_ensembl_havana_CCDS_transcriptsID.bed -b chromhmm_cmips.bed *cleaned.txt ATAC/DiffBind/IPS_diff_ATAC_fdr-0.01_fold-1.csv ATAC/DiffBind/CM_diff_ATAC_fdr-0.01_fold-1.csv -wb -names cmips *cleaned.txt iPS CMiPS | cut -f4,5,9 | sed 's/\.cleaned.txt//g' > tss_ccds_intersect_all_chromhmm_CM-IPS_atac.txt



bedtools intersect -a ../external/TSS_ensembl_havana_CCDS_transcriptsID.bed -b chromhmm_cmips.bed *cleaned.txt ATAC/DiffBind/CONTROL_diff_ATAC_block_fdr-0.1_fold-1.csv ATAC/DiffBind/BRUG_diff_ATAC_block_fdr-0.1_fold-1.csv  -wb -names cmips *cleaned.txt CONTROL BRUG | cut -f4,5,9 | sed 's/\.cleaned.txt//g' > tss_ccds_intersect_all_chromhmm_CONTROL-BRUG_atac.txt






