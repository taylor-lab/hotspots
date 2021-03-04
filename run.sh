# To capture everything and then filter later on
./hotspot_algo.R \
	--input-maf=ERBBs/Data/ERBBS.maf \
	--rdata=hotspot_algo.Rdata \
	--output-file=ERBBs/sig_hotspots.txt \
	--homopolymer=FALSE \
	--qval=1

mv putative_false_positive_mutations.txt ERBBs/.

# ./hotspot_algo.R \
# 	--input-maf=ERBBs/Data/ERBBS.maf \
# 	--rdata=hotspot_algo.Rdata \
# 	--gene-query=ERBBs/Data/genes.txt \
# 	--output-file=sig_hotspots.txt \
# 	--homopolymer=FALSE