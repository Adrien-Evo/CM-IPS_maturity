#Notes

## On the analyses restricted on the TSS of the ensembl_havana and CCDS

The scripts to get the tss is here:
src/data/get_tss.sh

Its gives 37309 tss for 18686 unique gene names (more than 1 transcript per gene). When looking for the chromHMM states overlapping those tss, we would expect to get 18686 regions (radio$chromhmm_homer_tss_ccds). But it's only 16595 :
```
cmips = read.table(radio$chromhmm_homer_tss_ccds,sep="\t", header=T, quote="")
overlap=(cmips$End- cmips$Start)/2 > abs(cmips$Distance.to.TSS)
cmips_overlap = cmips[overlap,]
```

But there is tss that are very close from each other but on different strand. Example : DHX38 . So some chromHMM annotation overlap more than one tss. That's why you have less lines that expected
