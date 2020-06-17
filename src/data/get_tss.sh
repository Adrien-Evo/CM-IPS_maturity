#This command query the ensemble gtf GRCh37, get all the transcript that are deteected both by the automaticcally annotated ensembl pipe and HAVAN manually curated. Transcripts are then filtered to only include the CCDS that are  member of the consensus CDS gene set, confirming coding regions between ENSEMBL, UCSC, NCBI and HAVANA (CCDS tag)

wget -q -O - "http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz" | gunzip -c | grep -v "#" | awk '($3=="transcript")' | awk '{if($2="ensembl_havana"){print $0}}' | grep -E "tag \"CCDS\"" | awk '{OFS="\t"};{if($7 == "+"){start = $4} else if($7 == "-"){start = $5}};{print $24,$1, start-1, start,"0"}' | sed -r 's/"|;//g' > data/external/TSS_ensembl_havana_CCDS_transcripts_peakfile_for_HOMER.txt


awk '{OFS="\t"}{print $2,$3,$4,$1}' data/external/TSS_ensembl_havana_CCDS_transcripts_peakfile_for_HOMER.txt > data/external/TSS_ensembl_havana_CCDS_transcripts.bed


wget -q -O - "http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz" | gunzip -c | grep -v "#" | awk '($3=="transcript")' | awk '{if($2="ensembl_havana"){print $0}}' | grep -E "tag \"CCDS\"" | awk '{OFS="\t"};{if($7 == "+"){start = $4} else if($7 == "-"){start = $5}};{print $14,$1, start-1, start,"0"}' | sed -r 's/"|;//g' > data/external/TSS_ensembl_havana_CCDS_transcripts_ENST_peakfile_for_HOMER.txt


awk '{OFS="\t"}{print $2,$3,$4,$1}' data/external/TSS_ensembl_havana_CCDS_transcripts_ENST_peakfile_for_HOMER.txt > data/external/TSS_ensembl_havana_CCDS_transcriptsID.bed
