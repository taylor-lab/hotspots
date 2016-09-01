
# functions for hotspot_algo.R

# returns baseline expected probability given gene length, gene mutation burden, and total samples
get.probability=function(gene,aa,total.muts,total.samples,aa.length) {
	top=sum(aa$count[which(aa$toohot)])
	aa.length=aa.length-length(which(aa$toohot))
	total.muts=total.muts-top
	return((1/aa.length)*(total.muts/total.samples))
}

# returns the mutability of the codon given the gene
get.alpha=function(k,aa,mu_prot){
	if(is.na(aa$mu_position[k])) return(1)
	else return(aa$mu_position[k]/mu_prot)
}

# returns log10 p-value of the binomial distribution
get.pvalues=function(k,total.muts,aa.length,mu_prot,aa,gene,min_prob,total.samples=TOTAL_SAMPLES) {
	prob=max(get.probability(gene,aa,total.muts,total.samples,aa.length),min_prob)
	alpha=get.alpha(k,aa,mu_prot)
	prob=prob*alpha
	pbinom(as.numeric(as.character(aa$count[k]))-1,size=total.samples,prob=prob,
		lower.tail=FALSE,log.p=TRUE)/log(10,exp(1))
}

# returns the table and column names of table into a single string 
combine=function(tb,sep=":") {
	out=paste(names(tb)[1],tb[1],sep=sep)
	if(length(tb)>1) for(i in 2:length(tb)) out=paste(out,paste(names(tb)[i],tb[i],sep=sep),sep="|")
	return(out)
}

# returns the number and types of mutant residues of the codon
get.variant.aa=function(position,mf) {
	mf=mf[ which(mf$Amino_Acid_Position==position), ]
	variant=rev(sort(table(mf$Variant_Amino_Acid)))
	variant=combine(variant)
	codon=rev(sort(table(mf$Codons)))
	codon=combine(codon)
	genomic_pos=rev(sort(table(paste(mf$Chromosome,mf$Start_Position,sep=":"))))
	genomic_pos=combine(genomic_pos,sep="_")
	return(c(variant,codon,genomic_pos))
}

# returns the number of cancer types with the mutation
get.cancer.type=function(position,mf,d=d) {
	mf=mf[ mf$Amino_Acid_Position==position, ]
	cnt=c()
	for(i in unique(mf$TUMORTYPE)) {
		tm=mf[ mf$TUMORTYPE==i, ]
		cnt=rbind(cnt,c(i,length(unique(d$Master_ID[ d$TUMORTYPE==i ])),nrow(mf[ mf$TUMORTYPE==i,])))
	}
	cnt=as.data.frame(cnt)
	colnames(cnt)=c('tumortype','samples_muted','mutations')
	cnt$mutations=as.numeric(as.character(cnt$mutations))
	cnt=cnt[ order(-cnt$mutations), ]
	output=paste(cnt$tumortype[1],cnt$samples_muted[1],cnt$mutations[1],sep=':')
	if(nrow(cnt)>1) for(i in 2:nrow(cnt)) output=paste(output,paste(cnt$tumortype[i],cnt$samples_muted[i],cnt$mutations[i],sep=':'),sep="|")
	return(output)
}

# returns the count of tri-nucleotide contexts of the codon
get.reference.tri=function(pos,maf) {
	table(maf$Ref_Tri[ which(maf$Amino_Acid_Position==pos) ])
}

# return the weighted average tri-nucleotide mutability of the codon
get.mu.score=function(pos,maf,mu) {
	# given amino acid position, get count of all the reference tri-nucleotide context
	t=get.reference.tri(pos,maf)
	# get the mutability scores of those tri-nucleotides
	ind=match(names(t),mu$tri)
	# return the weighted average 
	return(sum(mu$mu[ind]*t)/sum(t))
}

# return the counts and types of tri-nucleotide contexts 
get.all.tri=function(pos,maf){
	t=rev(sort(table(maf$Ref_Tri[ which(maf$Amino_Acid_Position==pos) ])))
	return(paste(names(t),collapse="|"))
}

# returns the reference amino acid 
get.reference.aa=function(pos,maf) {
	tb=table(maf$Reference_Amino_Acid[ which(maf$Amino_Acid_Position==pos) ])
	tb=tb[ order(tb,decreasing=T) ]
	return(combine(tb))
}

# returns dbSNP/1000G/NHLBI ids, if any
get.rsid=function(pos,maf) {
	tb=table(maf$dbSNP_RS[ which(maf$Amino_Acid_Position==pos) ])
	tb=tb[ order(tb,decreasing=T) ]
	if(length(tb)==1) if(names(tb)=='') return('')
	return(combine(tb))
}

# return perentage into decimal
convert.to.decimal=function(num){
	if(is.na(num)) return(NA)
	if(grepl('%',num)) num=gsub('%','',num)
	num=as.numeric(as.character(num))
	if(num > 1) num=num/100
	return(num)
}

# return the variant allele frequency (VAF) ranks of the mutation in the sample
get.sample.rank=function(sample,gene,pos,maf) {
	sample=d[ which(d$Master_ID==sample), ]
	sample$allele_freq=as.numeric(as.character(unlist(lapply(sample$allele_freq,convert.to.decimal))))
	af.rank=sort(unique(sample$allele_freq),decreasing=T)
	return(which(af.rank==max(sample$allele_freq[ which(sample$Hugo_Symbol==gene & sample$Amino_Acid_Position==pos) ]))/length(af.rank))
}

# return all the VAF ranks of the mutation in all samples with the mutation
get.af.rank=function(pos,maf) {
	samples=maf$Master_ID[ which(maf$Amino_Acid_Position==pos & !is.na(maf$allele_freq)) ]
	if(length(samples)==0) return(c(NA,NA))
	af.rank=unlist(lapply(samples,get.sample.rank,gene=unique(maf$Hugo_Symbol),pos=pos,maf=maf))
	return(c(median(af.rank),paste(length(af.rank),paste(af.rank,collapse=':'),sep="|")))
}

# return the cancer cell fractions (CCF) of the mutation
get.ccf=function(pos,maf) {
	ccfs=maf$ccf[ which(maf$Amino_Acid_Position==pos & !is.na(maf$ccf)) ]
	if(length(ccfs)==0) return(NA)
	return(paste(ccfs,collapse=":"))
}

# perform binomial test for all mutations observed in gene
binom.test_snp=function(gene) {

	# Reduce maf down to just of gene G
	maf=d[ which(d$Hugo_Symbol==gene), ]
	# Get count of mutations as each amino acid position
	aa=table(maf$Amino_Acid_Position)
	# Get amino acid length of the gene
	aa.length=max(maf$Protein_Length)

	# Find the mutability of the gene (pre-calculated)
	# If gene does not exist, use the average mutability of all genes 
	ii=which(p$gene==gene)
	if(any(ii)) mu_protein=p$score[ p$gene==gene ]	else mu_protein=mean(p$score)

	# Find the mutability of the amino acid position
	mu_position=unlist(lapply(as.numeric(names(aa)),get.mu.score,maf=maf,mu=mu))
	aa=cbind(as.data.frame(aa),mu_position)
	colnames(aa)=c('pos','count','mu_position')
	aa$mu_position=as.numeric(as.character(aa$mu_position))
	aa$pos=as.numeric(as.character(aa$pos))
	aa$count=as.numeric(as.character(aa$count))
	aa$tri=unlist(lapply(aa$pos,get.all.tri,maf=maf))

	af.ranks=as.data.frame(do.call('rbind',lapply(aa$pos,get.af.rank,maf=maf)))
	colnames(af.ranks)=c('Median_Allele_Freq_Rank','Allele_Freq_Rank')
	if('ccf' %in% colnames(maf)) ccfs=unlist(lapply(aa$pos,get.ccf,maf=maf))
	else ccfs=rep(NA,nrow(aa))

	aa$toohot=FALSE
	cutoff=max(as.numeric(quantile(aa$count,0.99)),20)
	aa$toohot[ which(aa$count > cutoff) ]=TRUE

	# cat('      Compiling output...\n')
	info=as.data.frame(do.call('rbind',lapply(aa$pos,get.variant.aa,mf=maf)))
	colnames(info)=c('Variant_Amino_Acid','Codon_Change','Genomic_Position')
	cancer.types=as.data.frame(do.call('rbind',lapply(aa$pos,get.cancer.type,mf=maf,d=d)))
	colnames(cancer.types)='Samples'
	output=cbind('Hugo_Symbol'=rep(gene,nrow(aa)),'Amino_Acid_Position'=aa$pos,
		'log10_pvalue'=unlist(lapply(1:nrow(aa),get.pvalues,total.muts=sum(aa$count),aa.length=aa.length,mu_prot=mu_protein,aa=aa,gene=gene,min_prob=min_prob)),
		'Mutation_Count'=aa$count,'Reference_Amino_Acid'=unlist(lapply(aa$pos,get.reference.aa,maf=maf)),
		'Total_Mutations_in_Gene'=rep(sum(aa$count),nrow(aa)),'Median_Allele_Freq_Rank'=af.ranks$Median_Allele_Freq_Rank,
		'Allele_Freq_Rank'=af.ranks$Allele_Freq_Rank,'SNP_ID'=unlist(lapply(aa$pos,get.rsid,maf=maf)),info,'Tumortypes'=cancer.types)
	output=cbind(output,'Tri-nucleotides'=aa$tri, 'Mutability'=aa$mu_position,'mu_protein'=rep(mu_protein,nrow(aa)))
	output=cbind(output,'ccf'=ccfs)
	output$Genomic_Position=as.character(output$Genomic_Position)

	pval=10^output$log10_pvalue
	ppval=c(pval,rep(1,aa.length-length(pval)))
	q2=p.adjust(ppval,method='BY')[1:length(pval)]
	output=cbind(output,'qvalue'=q2)

	# sorting by q-value
	output=output[ order(output$qvalue), ]
	return(output)
}

# return annotation track of hotspots with alignability/uniqueness scores (UCSC)
annotate.bedGraph.tracks=function(df,bedgraph) {
	# formating bedGraph file
	setnames(bedgraph,old=colnames(bedgraph),new=c('chromosome','start','stop','score'))
	bedgraph$chromosome=gsub('chr','',bedgraph$chromosome)
	bedgraph$chromosome[ bedgraph$chromosome=='X' ]=23
	bedgraph$chromosome[ bedgraph$chromosome=='Y' ]=24
	bedgraph$chromosome=suppressWarnings(as.numeric(as.character(bedgraph$chromosome)))
	bedgraph=bedgraph[ which(!is.na(bedgraph$chromosome)), ]
	bedgraphgrange=RangedData(IRanges(bedgraph$start,bedgraph$stop),space=bedgraph$chromosome,values=bedgraph$score)

	# collect all the genomic positions in hotspot list
	pos=matrix(unlist(lapply(df$Genomic_Position,function(x) unlist(strsplit(x,'[|:_]')))),ncol=3,byrow=TRUE)[,1:2]
	pos=as.data.frame(pos)
	colnames(pos)=c('Chromosome','Start_Position')
	pos$Chromosome=as.character(pos$Chromosome)
	pos$Chromosome[ which(pos$Chromosome=='X') ]=23
	pos$Chromosome[ which(pos$Chromosome=='Y') ]=24
	pos$Chromosome=as.numeric(as.character(pos$Chromosome))
	pos$Start_Position=as.numeric(as.character(pos$Start_Position))
	dirange=split(IRanges(pos$Start_Position,pos$Start_Position),pos$Chromosome)

	overlap=findOverlaps(dirange,bedgraphgrange)
	overlap=as.data.frame(as.matrix(overlap))
	overlap$queryHits=as.numeric(as.character(overlap$queryHits))
	overlap$subjectHits=as.numeric(as.character(overlap$subjectHits))
	dirange=as.data.frame(dirange)
	colnames(dirange)=c('n','chromosome','start','stop','width')
	dirange$width=NULL
	dirange$n=NULL
	dirange$stop=NULL
	dirange$chromosome=as.numeric(as.character(dirange$chromosome))
	dirange$start=as.numeric(as.character(dirange$start))

	bedgraphgrange=as.data.frame(bedgraphgrange)
	colnames(bedgraphgrange)=c('chromosome','start','end','width','score')
	bedgraphgrange$width=NULL

	# get alignability score from IRanges object
	get.score=function(i) {
		z=which(overlap$queryHits==i)
		if(length(z)==0) return(NA)
		s=max(bedgraphgrange$score[overlap$subjectHits[ z ]])
		return(s)
	}
	scores=unlist(lapply(1:nrow(dirange),get.score))
	dirange=cbind(dirange,scores)

	# returns matrix of genomic positions (chromosome, start position)
	get.genomic.pos=function(i,df) {
		return(	matrix(unlist(strsplit(df$Genomic_Position[i],'[|:_]')),ncol=3,byrow=TRUE))
	}

	# collect alignment scores
	final_scores=c()
	for(ii in 1:nrow(df)) {
		mx=get.genomic.pos(ii,df)
		mx[,1]=gsub('X',23,mx[,1])
		mx[,1]=gsub('Y',24,mx[,1])
		codon_score=c()
		for(p in 1:nrow(mx)) {
			yy=which(dirange$chromosome==mx[p,1] & dirange$start==mx[p,2])
			codon_score=c(codon_score,rep(dirange$scores[yy],mx[p,3]))
		}
		final_scores=c(final_scores,mean(codon_score))
	}
	return(final_scores)
}

# return mutation calling centers (TCGA) of a given hotspot
mut.center.breakdown=function(i,sig,maf) {

	get.genomic.position=function(x,hs) {
		pos=do.call('rbind',strsplit(unlist(strsplit(hs$Genomic_Position[x],'\\|')),'_'))
		pos=matrix(do.call('rbind',strsplit(pos[,1],':')),ncol=2)
		out=c()
		for(i in 1:nrow(pos)) out=c(out,paste(pos[i,],collapse=" "))
		return(paste(out,collapse=":"))
	}

	get.sequencing.centers=function(pos,maf) {
		b=pos
		chr=unlist(strsplit(pos,' '))[1]
		pos=as.numeric(unlist(strsplit(pos,' '))[2])
		maf=maf[ which( as.character(maf$Chromosome)==chr & maf$Start_Position==pos& maf$End_Position==pos), ]
		seqcenter=maf$Center[ which(!duplicated(maf$Tumor_Sample_Barcode)) ]
		if(length(seqcenter)==0) return(c(NA,NA,NA,NA,NA,NA,NA))
		return(seqcenter)
	}

	pos=get.genomic.position(i,hs=sig)
	pos=unlist(strsplit(pos,':'))
	out=unlist(lapply(pos,get.sequencing.centers,maf=maf))
	out=out[ which(out!='') ]
	tt=out
	out=table(out)
	return(c(sig$Hugo_Symbol[i],sig$Amino_Acid_Position[i],names(out)[ which(out==max(out))][1], max(out),sum(out),max(out)/sum(out),paste(tt,collapse=":")))
}

# return TRUE|FALSE of hotspots with putative mutation calling center bias
annotate.center.bias=function(sig,maf) {

	# initialize variables
	sig$Samples=as.character(sig$Samples)
	sig$Hugo_Symbol=as.character(sig$Hugo_Symbol)
	sig$Amino_Acid_Position=as.numeric(as.character(sig$Amino_Acid_Position))
	sig$tm=paste(sig$Hugo_Symbol,sig$Amino_Acid_Position)

	maf$tm=paste(maf$Hugo_Symbol,maf$Amino_Acid_Position)

	major.tumortype=function(x,hs) {
		samples=do.call('rbind',strsplit(unlist(strsplit(hs$Samples[x],'\\|')),':'))
		count=as.numeric(as.vector(samples[,3]))
		return(c(samples[1,1],count[1],sum(count)))
	}
	mjr=lapply(1:nrow(sig),major.tumortype,hs=sig)
	mjr=as.data.frame(do.call('rbind',mjr))
	colnames(mjr)=c('tt','count','tot')
	mjr$tt=as.character(mjr$tt)
	mjr$count=as.numeric(as.character(mjr$count))
	mjr$tot=as.numeric(as.character(mjr$tot))
	mjr$ratio=mjr$count/mjr$tot
	ind=which( mjr$tt %in% 'coadread' & mjr$ratio >= 0.5)

	tm=suppressWarnings(do.call('rbind',lapply(ind,mut.center.breakdown,sig=sig,maf=maf)))
	tm=as.data.frame(tm)
	colnames(tm)=c('gene','aa_position','major_calling_center','from_major','total','ratio','seq_centers')
	for(i in 1:ncol(tm)) tm[,i]=as.character(tm[,i])
	tm$ratio=as.numeric(as.character(tm$ratio))
	tm$from_major=as.numeric(as.character(tm$from_major))
	tm$total=as.numeric(as.character(tm$total))

	tm$ratio=tm$from_major/tm$total
	tm=tm[ order(tm$ratio,decreasing=T), ]
	tm=tm[ which(!grepl(',',tm$major_calling_center)), ]
	tm$idx=paste(tm$gene,tm$aa_position)
	tm$num_tt=unlist(lapply(tm$idx,function(x) length(unique(maf$TUMORTYPE[ which(maf$tm==x) ]))))

	centerbias=rep(FALSE,nrow(sig))
	# return TRUE if there are at least 7 mutations with known mutation calling centers and
	# more than 85% of mutations were called from a single mutation calling center or
	# if mutaitons are localized to only two tumor types whose majority comes from greater than 75%
	ii=which(((tm$ratio>0.85)|(tm$ratio==0.75&tm$num_tt<3)) &tm$from_major>5)
	centerbias[ which(sig$tm%in%tm$idx[ii]) ]=TRUE
	return(centerbias)

}

# return calculcated entropy around the site of each hotspot
annotate.entropy=function(sig) {
	
	# find genomic position of hotspot
	hspos=lapply(strsplit(sig$Genomic_Position,"[|]"),function(x) gsub("_.*$","",x))
	hspos=do.call(rbind,lapply(hspos,function(x) c(gsub(":.*$","",x[1]),range(as.numeric(gsub("^.*:","",x))))))
	hspos[,1]=paste("chr",hspos[,1],sep="")
	pcN=rep(1/4,4)
	names(pcN)=c("A","C","T","G")
	n=nrow(hspos)

	# make the intervals to get sequence
	regions=c(rep("up12",n),rep("up24",n),rep("up36",n),rep("dn12",n),rep("dn24",n), rep("dn36",n))
	starts=c(as.numeric(hspos[,2])-11,as.numeric(hspos[,2])-23,as.numeric(hspos[,2])-35,rep.int(as.numeric(hspos[,3]),3))
	ends=c(rep.int(as.numeric(hspos[,2]),3),as.numeric(hspos[,3])+11,as.numeric(hspos[,3])+23,as.numeric(hspos[,3])+35)
	contigs=rep.int(hspos[,1],6)
	for(pad in c(11,23,35)) {
		rcontigs=sample(seqnames(Hsapiens)[1:23],1000,replace=T)
		rstarts=unlist(lapply(rcontigs,function(x) sample(seqlengths(Hsapiens)[x]-(pad+1),1)))
		rends=rstarts+pad
		starts=c(starts,rstarts)
		ends=c(ends,rends)
		contigs=c(contigs,rcontigs)
		regions=c(regions,rep(paste("rand",pad+1,sep=""),1000))
	}
	rpts=GRanges(seqnames=contigs,ranges=IRanges(starts,ends),names=regions)
	rpts_seq=getSeq(Hsapiens,rpts)
	# calculate entropy of returned sequence
	pN=lapply(strsplit(as.character(rpts_seq),""),function(x) table(x)+(1/4))
	rpts_entropy=unlist(lapply(pN,function(x) -sum((x/sum(x))*log(x/sum(x)))))  
	for(pad in c(12,24,36)) {
		ui=which(rpts$names==paste("up",pad,sep=""))
		di=which(rpts$names==paste("dn",pad,sep=""))
		sig[[paste("pad",pad,"entropy",sep="")]]=pmin(rpts_entropy[ui],rpts_entropy[di])
	}
	return(sig[,(ncol(sig)-2):ncol(sig)])
}

# return TRUE|FALSE of whether hotspot is a true positive
annotate.true.positives=function(sig,dmp) {

	# return TRUE if hotspot mutation is in true-positive file
	tp=rep(FALSE,nrow(sig))
	sig$id=paste(sig$Hugo_Symbol,sig$Amino_Acid_Position,sep="-")
	tp[sig$id%in%unique(dmp$id)]=TRUE
	return(tp)
}

# return homopolyer region, if any, of a given hotspot
get.repeats=function(i,sig,rr) {

	# initialize homopolymer file
	gene_pos=as.data.frame(matrix(unlist(strsplit(sig$Genomic_Position[i],'[_:|]')),ncol=3,byrow=T))
	colnames(gene_pos)=c('chromosome','pos','count')
	gene_pos$chromosome=paste('chr',as.character(gene_pos$chromosome),sep="")
	gene_pos$chromosome=as.character(gene_pos$chromosome)
	gene_pos$pos=as.numeric(as.character(gene_pos$pos))
	gene_pos$count=as.numeric(as.character(gene_pos$count))

	# find overlap padded by 1 position
	pos_range=RangedData(IRanges(gene_pos$pos-1,gene_pos$pos+1),space=gene_pos$chromosome,values=gene_pos$count)
	hits=as.data.frame(as.matrix(findOverlaps(pos_range,rr)))
	if(nrow(hits)==0) return(c(FALSE,NA,NA))
	else {
		results=as.data.frame(rr)[as.numeric(as.character(hits$subjectHits[1:nrow(hits)])),]
		ret=unlist(strsplit(results$values[1],' '))
		return(c(TRUE,ret[1],ret[2]))
	}

}

# return homopolymer region metrics: TRUE|FALSE if it is a repeat region, the repeat nucleotide(s), and length of repeat region
annotate.homopolymers=function(sig,homopolymerbed) {

	# find overlap
	rep_range=RangedData(IRanges(homopolymerbed$start,homopolymerbed$end),
		space=homopolymerbed$chromosome,values=paste(homopolymerbed$rep,nchar(homopolymerbed$seq)))
	out=lapply(1:nrow(sig),get.repeats,sig=sig,rr=rep_range)
	bb=do.call('rbind',out)
	colnames(bb)=c('Is_repeat','seq','length')
	return(bb)

}

# return filtering judgement based on annotated metrics
annotate.filtering.judgement=function(sig) {

	reason=rep('',nrow(sig))
	sig$Hugo_Symbol=as.character(sig$Hugo_Symbol)

	# sub-clonality filter
	ccfs=lapply(strsplit(as.character(sig$ccf),":"),as.numeric)
	sig$ccf_subclonal_frac=unlist(lapply(ccfs,function(x) sum(x<0.8)/length(x)))
	sig$num_ccf=unlist(lapply(ccfs,function(x) length(x)))
	rem=which(sig$num_ccf>1 & sig$ccf_subclonal_frac==1)
	threshold=max(sig$ccf_subclonal_frac[sig$TP & sig$num_ccf>1])
	xi=which(sig$num_ccf>1 & sig$ccf_subclonal_frac>threshold)
	reason[xi]=ifelse(reason[xi]=="","SUBCLONAL_FRAC",paste(reason[xi],"SUBCLONAL_FRAC",sep="|"))

	# low af filter
	sig$Allele_Freq_Rank=as.character(sig$Allele_Freq_Rank)
	sig$num_af=as.numeric(as.character(unlist(lapply(sig$Allele_Freq_Rank,function(x) unlist(strsplit(x,'\\|'))[1] ))))
	sig$Allele_Freq_Rank=unlist(lapply(sig$Allele_Freq_Rank,function(x) unlist(strsplit(x,'\\|'))[2] ))
	afs=lapply(strsplit(sig$Allele_Freq_Rank,":"),as.numeric)
	sig$low_af_frac=unlist(lapply(afs,function(x) sum(x>0.8)/length(x)))
	threshold=max(sig$low_af_frac[which(sig$TP & sig$num_af>5)])
	xi=which(sig$num_af>5 & sig$low_af_frac>threshold)
	reason[xi]=ifelse(reason[xi]=="","LOW_AF_FRAC",paste(reason[xi],"LOW_AF_FRAC",sep="|"))


	# center bias filter
	if('seq_bias' %in% colnames(sig)) {
		xi=which(sig$seq_bias)
		reason[xi]=ifelse(reason[xi]=="","CENTER_BIAS",paste(reason[xi],"CENTER_BIAS",sep="|"))
	}

	# entropy filter
	xi=which(sig$pad12entropy<min(sig$pad12entropy[sig$TP]) | sig$pad24entropy<min(sig$pad24entropy[sig$TP]))
	reason[xi]=ifelse(reason[xi]=="","LOCAL_ENTROPY",paste(reason[xi],"LOCAL_ENTROPY",sep="|"))

	# homopolymer filter
	if('Is_repeat' %in% colnames(sig) & 'seq' %in% colnames(sig) & 'length' %in% colnames(sig)) {
		xi=which(sig$Is_repeat & sig$length>=10)
		reason[xi]=ifelse(reason[xi]=="","HOMOPOLYMER",paste(reason[xi],"HOMOPOLYMER",sep="|"))
	}

	# mapability filter
	if('align100'%in%colnames(sig) & 'align24'%in%colnames(sig)) {		
		ranked_align=rowSums(cbind(rank(sig$align100),rank(sig$align24)))
		xi=which(ranked_align<min(ranked_align[sig$TP]) | sig$align24<quantile(sig$align24,probs=0.125))
		reason[xi]=ifelse(reason[xi]=="","ALIGNABILITY",paste(reason[xi],"ALIGNABILITY",sep="|"))
	}

	# biologically-driven filters of putative false positives
	retained=table(sig$Hugo_Symbol[reason==""])
	excluded=table(sig$Hugo_Symbol[reason!=""])
	xi=names(excluded)[which(names(excluded)%in%names(retained))]
	de=list()
	de$excluded=excluded[xi]
	de$retained=retained[xi]
	de=as.data.frame(de)
	xi=which(sig$Hugo_Symbol%in%rownames(de)[de$retained<=(de$excluded*3.5)] & reason=="")
	reason[xi]=ifelse(reason[xi]=="","FP_RICH_GENE",paste(reason[xi],"FP_GENE",sep="|"))
	
	supp_table_7_mutsigcv_paper=c(
		"PCLO","FLG","BAGE2","TPTE","TTN","CSMD1","CSMD3","RYR2","RYR3","MUC16","MUC4",
		"MUC17","MUC5B","OR10G9","OR2G6","OR4C6","OR4M2","OR2T4","OR5L2","OR2T33"
	)
	xi=which(sig$Hugo_Symbol%in%supp_table_7_mutsigcv_paper)
	reason[xi]=ifelse(reason[xi]=="","FP_MUTSIG",paste(reason[xi],"FP_MUTSIG",sep="|"))

	return(reason)
}


split.double=function(ind,df) {
	ref=unlist(strsplit(df$Reference_Amino_Acid[ind],''))
	vart=unlist(strsplit(df$Variant_Amino_Acid[ind],''))
	out=c()
	for(i in 1:length(ref)) {
		tm=df[ind,]
		tm$Amino_Acid_Position=tm$Amino_Acid_Position+(i-1)
		tm$Reference_Amino_Acid=ref[i]
		tm$Variant_Amino_Acid=vart[i]
		if(ref[i]=='*') tm$Consequence='stop_lost'
		else if(vart[i]=='*') tm$Consequence='stop_gained'
		else if(ref[i]==vart[i]) tm$Consequence='synonymous_variant'
		else tm$Consequence='missense_variant'
		tm$Amino_Acid_Change=paste(tm$Reference_Amino_Acid,tm$Amino_Acid_Position,tm$Variant_Amino_Acid,sep="")
		out=rbind(out,tm)
	}
	return(out)
}

#returns other positions that are 1 bp away from pos
find.splice.neighbors=function(pos,df) {
	dif=abs(df$Amino_Acid_Position-pos)
	ind=which(dif==1)
	if(length(ind)==0) return(NA)
	return(unique(c(pos,df$Amino_Acid_Position[ind])))
}
#returns TRUE if there are other positions that are 1 bp away
any.neighbors=function(df) {
	tm=sort(unique(df$Amino_Acid_Position))
	tt=lapply(tm,find.splice.neighbors,df=df)
	if(length(which(!is.na(tt))) > 0) return(TRUE)
	return(FALSE)
}

#returns the first occurance of a neighbor
get.neighbors=function(df) {
	tm=sort(unique(df$Amino_Acid_Position))
	tt=lapply(tm,find.splice.neighbors,df=df)
	return(tt[[which(!is.na(tt))[1]]])
}
remove.unexpressed.genes=function(i,maf,express) {
	
	#coad samples
	if(maf$TUMORTYPE[i]=='coadread') {
		tt=which(colnames(express)=='coad' | colnames(express)=='read')
		g=which(express$gene_id==maf$Hugo_Symbol[i])
		if(!any(tt) | !any(g) ) return(1==2)
		else {
			explvl=as.character(express[ g, tt ])
			explvl=do.call('rbind',strsplit(explvl,';'))
			samples=as.numeric(explvl[,1])
			tot=as.numeric(explvl[,3])
			return(samples < 0.1*tot[1] && samples < 0.1*tot[2])
		}
	}
	
	tt=grep(maf$TUMORTYPE[i],colnames(express))
	g=which(express$gene_id==maf$Hugo_Symbol[i])
	if(!any(g)) return(1==2)
	#if tumortype is not included, then i will infer gene expressions based on other tumortypes
	#if gene is not expressed in 95% of all tumor samples, then i infer this gene is not expressed
	#in the tumor samples we do not have RNASeqV2 data for 
	else if(!any(tt)) {
		explvl=express[ g, ]
		explvl=do.call('rbind',strsplit(as.character(explvl[2:length(explvl)]),';'))
		samples=sum(as.numeric(as.character(explvl[,1])))
		tot=sum(as.numeric(as.character(explvl[,3])))
		return(samples < 0.05*tot)
	}
	else {
		explvl=express[ g, tt ]
		explvl=unlist(strsplit(explvl,';'))
		samples=as.numeric(explvl[1])
		tot=as.numeric(explvl[3])
		return(samples < 0.1*tot)
	}

}

# deprecated germline SNP filtering based on 1000/NHLBI
# putative germline SNPs are filtered based on ExAC with minor AF > 0.06%
# ExAC r0.2 has been preloaded 
remove.snps=function(maf) {
	return(maf[ which(!paste(maf$Chromosome,maf$Start_Position,maf$Tumor_Seq_Allele2) 
		%in% paste(exacr0_2snps$Chromosome,exacr0_2snps$Position,exacr0_2snps$Alt)), ])
}

remove.unexpressed.mutations=function(maf,expressiontb) {
	ind=unlist(lapply(1:nrow(maf),remove.unexpressed.genes,maf=maf,express=expressiontb))
	ind=which(ind)
	if(any(ind)) maf=maf[ -ind, ]
	return(maf)
}



prepmaf=function(maf,expressiontb) {

	cat('Prepping MAF for analysis ...\n')
	#only non-indel, coding mutations
	coding=c('initiator_codon_variant','missense_variant','splice_acceptor_variant',
		'splice_donor_variant','stop_gained','stop_lost','stop_retained_variant','synonymous_variant')
	cat(' ... Reducing MAF to protein-coding substitutions\n')
	maf=maf[ which(maf$Consequence %in% coding & maf$CANONICAL=='YES' & maf$BIOTYPE=='protein_coding'), ]
	#remove indels
	maf=maf[ which(!maf$Variant_Type%in%c('INS','DEL')), ]

	# add additional annotations
	if(!'TUMORTYPE' %in% colnames(maf)) maf$TUMORTYPE='none'
	if(!'Master_ID' %in% colnames(maf)) maf$Master_ID=maf$Tumor_Sample_Barcode
	maf$Amino_Acid_Change=gsub('p.','',maf$HGVSp_Short)
	maf$Amino_Acid_Position=unlist(lapply(maf$Protein_position,function(x) unlist(strsplit(x,"/"))[1] ))
	maf$Protein_Length=as.numeric(unlist(lapply(maf$Protein_position,function(x) unlist(strsplit(x,'\\/'))[2])))
	maf$Reference_Amino_Acid=unlist(lapply(1:nrow(maf), function(x) unlist(strsplit(maf$Amino_acids[x],'/'))[1]))
	maf$Variant_Amino_Acid=unlist(lapply(1:nrow(maf), function(x) unlist(strsplit(maf$Amino_acids[x],'/'))[2]))
	maf$allele_freq=maf$t_alt_count/(maf$t_alt_count+maf$t_ref_count)

	#subset splice_site mutations
	splice=c('splice_acceptor_variant','splice_donor_variant')
	ss=maf[ which(maf$Consequence %in% splice), ]
	maf=maf[ which(!maf$Consequence %in% splice), ]

	cat(' ... Re-annotating substitutions\n')
	# annotating oligo mutations
	i=which(nchar(maf$Reference_Amino_Acid)!=1 & nchar(maf$Variant_Amino_Acid)!=1 & !grepl('ext',maf$Reference_Amino_Acid))
	if(any(i)) {
		dd=maf[i,]
		dd$Reference_Amino_Acid[ which(is.na(dd$Reference_Amino_Acid))]='NA'
		dd$Variant_Amino_Acid[ which(is.na(dd$Variant_Amino_Acid))]='NA'
		dd$Amino_Acid_Position=unlist(lapply(dd$Amino_Acid_Position,function(x) as.numeric(as.character(unlist(strsplit(x,'-'))[1]))))
		tm=do.call('rbind',lapply(1:nrow(dd),split.double,df=dd))
		maf=maf[-i, ]
		maf=rbind(maf,tm)
	}
	# annotating splice site mutations
	ss=ss[ nchar(ss$Reference_Allele)<3, ]
	ss$Amino_Acid_Position=as.numeric(ss$Start_Position)
	ss$Reference_Amino_Acid='SS'
	ss$Variant_Amino_Acid='SS'
	out=c()
	for(i in unique(ss$Chromosome)) {
		tm=ss[ ss$Chromosome==i, ]
		while(any.neighbors(tm)) {
			neighbors=get.neighbors(tm)
			tm$Amino_Acid_Position[ tm$Amino_Acid_Position %in% neighbors ]=min(neighbors)
		}
		out=rbind(out,tm)
	}
	ss=out
	maf=rbind(maf,ss)

	# remove putative germline mutations
	cat(' ... Removing putative germline SNPs based on ExAC r0.2\n')
	maf=remove.snps(maf)

	# remove mutations in unexpressed genes
	cat(' ... Removing unexpressed genes\n')
	maf=remove.unexpressed.mutations(maf,expressiontb)
	maf$Amino_Acid_Position=as.numeric(maf$Amino_Acid_Position)
	maf$Variant_Type='SNP'

	return(maf)
}
