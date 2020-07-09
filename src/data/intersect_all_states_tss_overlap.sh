bedtools intersect -a ../external/TSS_ensembl_havana_CCDS_transcripts.bed -b chromhmm_cmips.bed *cleaned.txt -wb -names cmips *cleaned.txt | cut -f4,5,9 | sed 's/\.cleaned.txt//g' > tss_ccds_intersect_all_chromhmm.txt



bedtools intersect -a ../external/TSS_ensembl_havana_CCDS_transcriptsID.bed -b chromhmm_cmips.bed *cleaned.txt -wb -names cmips *cleaned.txt | cut -f4,5,9 | sed 's/\.cleaned.txt//g' > tss_ccds_intersect_all_chromhmmID.txt


