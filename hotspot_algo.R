#! /usr/bin/env Rscript 

source('funcs.R')
args = commandArgs(TRUE)

usage="\n  Usage: ./hotspot_algo.R
    --input-maf=[Required: mutation file]
    --rdata=[Required: Rdata object with necessary files for algorithm]
    --true-positive=[Required: List of known true-positives]
    --output-file=[Required: output file to print statistically significant hotspots]
    --gene-query=[Optional (default=all genes in mutation file): List of Hugo Symbol in which to query for hotspots]
    --homopolymer=[Optional: BED file of homopolyer regions in hg19 for false positive filtering]
    --align100mer=[Optional: BED file of hg19 UCSC alignability track for 100-mer for false positive filtering]
    --align24mer=[Optional: BED file of hg19 UCSC alignability track for 24-mer for false positive filtering]
    --filter-centerbias=[Optional (default=FALSE): TRUE|FALSE Filter mutations based on mutation calling center bias]
    --qval=[Optional (default=0.01) q-value cut-off]\n
"

# Initialize default parameters
TRUEPOSITIVE_FILTER=FALSE
HOMOPOLYMER_FILTER=FALSE
ALIGN_FILTER=FALSE
CENTER_BIAS_FILTER=FALSE
GENES_INTEREST=FALSE

# Verify that at least required parameters and their values are passed 
if(length(args)<3) {
	cat(usage)
	stop("Incorrect or missing required input!")
}

# Verify mutation file
idx=grep("--input-maf=",args)
if(is.integer(idx)) {
	maf_fn=gsub("--input-maf=","",args[idx])
	if(!file.exists(maf_fn)) {
		stop("Unable to find mutation file!")
	} 
} else {
	stop("Missing required --input-maf parameter!")
}

# Verify mutation file
idx=grep("--rdata=",args)
if(is.integer(idx)) {
	rdata_file=gsub("--rdata=","",args[idx])
	if(!file.exists(rdata_file)) {
		stop("Unable to find required Rdata file!")
	} 
} else {
	stop("Missing required --rdata parameter!")
}

# Verify output file destination
idx=grep("--output-file=",args)
if(is.integer(idx)) {
	output_fn=gsub("--output-file=","",args[idx])
} else {
	stop("Missing required --output-file parameter!")
}


# Check for true positive file
idx=grep("--true-positive=",args)
if(is.integer(idx)) {
	true_pos_fn=gsub("--true-positive=","",args[idx])
	if(!file.exists(true_pos_fn)) {
		stop("Unable to find true-positive file!")
	} 
} else {
	stop("Missing required --true-positives parameter!")
}

# Check for gene of interest file
idx=grep("--gene-query",args)
if(length(idx)!=0) {
	gene_interest_fn=gsub("--gene-query=","",args[idx])
	if(!file.exists(gene_interest_fn)) {
		stop("Unable to find gene-query file!")
	} 
	GENES_INTEREST=TRUE
} 

# Check for homopolymer file
idx=grep("--homopolymer",args)
if(length(idx)!=0) {
	homopolymer_fn=gsub("--homopolymer=","",args[idx])
	if(!file.exists(homopolymer_fn)) {
		stop("Unable to find homopolymer file!")
	} 
	HOMOPOLYMER_FILTER=TRUE
} 

# Check for 24-mer alignability file
idx=grep("--align24mer",args)
if(length(idx)!=0) {
	align24_fn=gsub("--align24mer=","",args[idx])
	if(!file.exists(align24_fn)) {
		stop("Unable to find 24-mer alignability file!")
	} 
} 

# Check for homopolymer file
idx=grep("--align100mer",args)
if(length(idx)!=0) {
	align100_fn=gsub("--align100mer=","",args[idx])
	if(!file.exists(align100_fn)) {
		stop("Unable to find align100mer file!")
	} 
	if(file.exists(align24_fn) & file.exists(align100_fn)) ALIGN_FILTER=TRUE
}

idx=grep("--qval",args)
if(length(idx)!=0) {
	QVAL=as.numeric(as.character(gsub("--qval=","",args[idx])))
	if(QVAL >1 | QVAL < 0 | is.na(QVAL) ) {
		cat('Invalid q-value entered. Using default q-value cutf-off, 0.01\n')
		QVAL=0.01
	}  
}
if(length(idx)==0) {
	cat('No q-value entered. Using default q-value cut-off, 0.01\n')
	QVAL=0.01
}


# Load necessary packages
options(warn=-1)
required="\n  Required packages:
	data.table
	BSgenome.Hsapiens.UCSC.hg19
	IRanges\n
"
if(!suppressMessages(library(data.table,logical.return=TRUE)) | 
	!suppressMessages(library(IRanges,logical.return=TRUE)) | 
	 !suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19,logical.return=TRUE))) {
		cat(required)
		stop("Required packages are not installed!")
}

# read in necessary files
cat('\n\nReading in files...\n')
load(rdata_file)
d=read.csv(maf_fn,header=T,as.is=T,sep="\t",comment.char='#')

# add additional annotations
d$Amino_Acid_Change=gsub('p.','',d$HGVSp_Short)
d$Amino_Acid_Position=unlist(lapply(d$Protein_position,function(x) unlist(strsplit(x,"/"))[1] ))
d$Protein_Length=unlist(lapply(d$Protein_position,function(x) unlist(strsplit(x,'\\/'))[2]))
d$Reference_Amino_Acid=unlist(lapply(1:nrow(d), function(x) unlist(strsplit(d$Amino_acids[x],'/'))[1]))
d$Variant_Amino_Acid=unlist(lapply(1:nrow(d), function(x) unlist(strsplit(d$Amino_acids[x],'/'))[2]))
d$allele_freq=d$t_alt_count/(d$t_alt_count+d$t_ref_count)

if(!'TUMORTYPE' %in% colnames(d)) d$TUMORTYPE='none'
if(!'Master_ID' %in% colnames(d)) d$Master_ID=d$Tumor_Sample_Barcode

# get tri nucleotide
write.table(d,'___temp_maf.tm',row.names=F,quote=F,sep="\t")
system(command='python make_trinuc_maf.py ___temp_maf.tm ___temp_maf-tri.tm')
# clean up tmp files
d=read.csv('___temp_maf-tri.tm',header=T,as.is=T,sep="\t",comment.char='#')
system(command='rm ___t*')

d=prepmaf(d,expressiontb)

TOTAL_SAMPLES=length(unique(d$Master_ID))

genes=unique(d$Hugo_Symbol)
if(GENES_INTEREST) { 
	genes=read.csv(gene_interest_fn,header=F,as.is=T,sep="\t")[,1]
	genes=genes[ which(genes %in% d$Hugo_Symbol) ]
	# ii=which(!genes %in% d$Hugo_Symbol)
	# if(length(ii) > 0) stop(paste(genes[ii],collapse=", "),' not mutated in mutation file\n')
}
# set minimum probabily for binomial model
min_prob=quantile(all_prob,0.2)

# run algorithm
cat('Running algorithm...\n')
output=lapply(genes,binom.test_snp)
output=do.call('rbind',output)

# filter output to only the significant hits
output=output[ which(output$qvalue < QVAL), ]
output=output[ order(output$qvalue), ]

# annotating significant hotspots with various metrics
# annotate alignability 24- and 100-mers scores
if(ALIGN_FILTER) {
	cat('Annotating alignability tracks...\n')
	# align 100-mer track
	align100=fread(align100_fn,header=F)
	output=cbind(output,'align100'=annotate.bedGraph.tracks(output,align100))
	# align 24-mer track
	align24=fread(align24_fn,header=F)
	output=cbind(output,'align24'=annotate.bedGraph.tracks(output,align24))
}
# annotate mutation calling center bias
if(CENTER_BIAS_FILTER) {
	cat('Annotating mutation calling center bias...\n')
	output=cbind(output,'seq_bias'=as.logical(annotate.center.bias(sig=output,maf=d)))
}
# annotating homopolymer 
if(HOMOPOLYMER_FILTER) {
	output=cbind(output,annotate.homopolymers(sig=output,fn=homopolymer_fn))
	output$Is_repeat=as.logical(output$Is_repeat)
	output$seq=as.character(output$seq)
	output$length=as.numeric(as.character(output$length))
}
# annotating sequence entropy
output=cbind(output,annotate.entropy(output))
# annotating known true positives based on input file
output=cbind(output,'TP'=as.logical(annotate.true.positives(output,true_pos_fn)))

# removing putative false positives
cat('Post-hoc filtering...\n')
output=cbind(output,'reason'=as.character(annotate.filtering.judgement(output)))

putative_false_positives=output[ which(output$reason!=''), ]
if(nrow(putative_false_positives)>0) write.table(putative_false_positives,'putative_false_positive_mutations.txt',quote=F,row.names=F,sep="\t")

output=output[ which(output$reason==''), ]

# write output
cat('Writing output...\n')
write.table(output,output_fn,quote=F,row.names=F,sep="\t")








