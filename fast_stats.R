getStats = function(out, popDF, extent) {

	freqs = alleleFreqs(out, TRUE)
	popid = strataNames(out)
	split_out = strataSplit(out)

	freqmat = matrix(data = NA, nrow = length(freqs), ncol = length(popid))
	for(p in 1:length(freqmat[1,])) {
		freqmat[,p] = sapply(freqs, FUN=function(x){a=x[,,popid[p]][3]})
	}
	colnames(freqmat) = popid

	variable_loci = matrix(data = 0, nrow = length(freqs), ncol = length(popid))
	for(p in 1:length(variable_loci[1,])) {
		variable_loci[freqmat[,p] > 0 & freqmat[,p] < 1 & !is.na(freqmat[,p]),p] = 1
	}
	colnames(variable_loci) = popid

	nind_mat = matrix(data = NA, nrow = length(freqs), ncol = length(popid))
	for(p in 1:length(nind_mat[1,])) {
		nind_mat[,p] = sapply(freqs, FUN=function(x){a=x[,,popid[p]][1]+x[,,popid[p]][2]})
	}
	colnames(nind_mat) = popid

	#Local SNP number
	localSNP = colSums(variable_loci)
	names(localSNP) = paste0("S.", popid)

	Smean = mean(localSNP)
	names(Smean) = "mean_S"

	Ssd = sd(localSNP)
	names(Ssd) = "sd_S"

	#Private SNP number
	privateSNP = c()
	for(p in 1:length(variable_loci[1,])) {
		privateSNP[p] = length(which(rowSums(variable_loci) == 1 & variable_loci[,p] == 1))
	}
	names(privateSNP) = paste0("pS.", popid)

	prSmean = mean(privateSNP)
	names(prSmean) = "mean_prS"

	prSsd = sd(privateSNP)
	names(prSsd) = "sd_prS"

	total_priv = sum(privateSNP)
	names(total_priv) = "tot_priv"

	#Pairwise private SNPs
	combos = combn(popid, 2)
	pwPriv = c()
	for(pair in 1:length(combos[1,])) {
		subset = which(nind_mat[,combos[1,pair]] > 0 & nind_mat[,combos[2,pair]] > 0)
		tmpvar = variable_loci[subset,c(combos[1,pair], combos[2,pair])]
		tmp = c(length(which(tmpvar[,1] == 1 & tmpvar[,2] == 0)), length(which(tmpvar[,2] == 1 & tmpvar[,1] == 0)))
		names(tmp) = c(paste0("pwPriv.",combos[1,pair],".",combos[2,pair]), paste0("pwPriv.",combos[2,pair],".",combos[1,pair]))
		pwPriv = c(pwPriv, tmp)
		rm(tmpvar, tmp, subset)
	}

	#He
	Hetmp = sapply(split_out, FUN=function(x){exptdHet(x)})
	He = colMeans(Hetmp, na.rm = TRUE)
	names(He) = paste0("He.",  colnames(Hetmp))

	HeALL = mean(exptdHet(out))
	names(HeALL) = "HeALL"

	#Simple Fst calculation...
	#pwFst = c()
	#for(pair in 1:length(combos[1,])) {
	#	subset = which(nind_mat[,combos[1,pair]] > 0 & nind_mat[,combos[2,pair]] > 0)
	#	Hs_tmp = (nind_mat[subset,combos[1,pair]]*Hetmp[subset,combos[1,pair]] + nind_mat[subset,combos[2,pair]]*Hetmp[subset,combos[2,pair]])/(nind_mat[subset,combos[1,pair]]+nind_mat[subset,combos[2,pair]])
	#	combo_freq = (nind_mat[subset,combos[1,pair]]*freqmat[subset,combos[1,pair]] + nind_mat[subset,combos[2,pair]]*freqmat[subset,combos[2,pair]])/(nind_mat[subset,combos[1,pair]]+nind_mat[subset,combos[2,pair]])
	#	Ht_tmp = 1-combo_freq^2-(1-combo_freq)^2
	#	Fst_tmp = (Ht_tmp-Hs_tmp)/Ht_tmp
	#	pwFst[pair] = mean(Fst_tmp, na.rm = TRUE)
	#	names(pwFst)[pair] = paste0(combos[1,pair],".",combos[2,pair])
	#	rm(subset, Hs_tmp, combo_freq, Ht_tmp, Fst_tmp)
	#}

	#Using hierfstat to do this instead
    genind_out = gtypes2genind(out)

	FstMat = pairwise.fst(genind_out, pop = genind_out$pop)
	pairFst = as.vector(FstMat)
	Fstnames = c()
	#This will be different...
	for(pid in 1:(length(popDF$id)-1)) {
		Fstnames = c(Fstnames, paste0(popDF$id[pid],".", popDF$id[(pid+1):length(popDF$id)]))
	}
	names(pairFst) = Fstnames

########### some spatially focused stats
###     get pop coords:
        pops=t(matrix(1:prod(extent),ncol=extent[1]))
        popcrd=data.frame(t(sapply(popDF$grid.cell,function(i) {which(pops==i,arr.ind=T)})))
        names(popcrd)=c("y","x")
        popDF=cbind(popDF,popcrd)
        
###     get pairwise distances
        pdist=as.matrix(dist(popDF[,c("x","y")]))
        colnames(pdist) <- popDF$id
        rownames(pdist) <- colnames(pdist)
        diag(pdist) <- NA
        pdist[upper.tri(pdist)] <- NA
        dsts = as.data.frame(as.table(pdist))
        dsts = dsts[complete.cases(dsts),]
        names(dsts) <- c("to","from","d")
        dsts$from <- as.character(dsts$from)
        dsts$to <- as.character(dsts$to)
        dsts <- dsts[order(dsts$to,dsts$from),]

print("this far")

        sum_stats_gi<-summary(genind_out)
        numAll=sum_stats_gi$pop.n.all
        numAll=data.frame(id=names(numAll),na=numAll,stringsAsFactors=F)

        he_by_pop<-as.vector(as.numeric(lapply(lapply(seppop(genind_out),summary),
                                               function(x) mean(x$Hexp))))

        hedf <- data.frame(he=he_by_pop, id=names(seppop(genind_out)),stringsAsFactors=F)

        popDF <- merge(merge(popDF,numAll),hedf)

print("this far2")
        
        na.lat.fit <- lm(na~y,popDF)
        na.lat.slope=c(coef(na.lat.fit)[2])
        na.lat.int=c(coef(na.lat.fit)[1])

        he.lat.fit <- lm(he~y,popDF)
        he.lat.slope=c(coef(he.lat.fit)[2])
        he.lat.int=c(coef(he.lat.fit)[1])

        fsts = data.frame(fst=pairFst)
        fsts = cbind(fsts,data.frame(t(sapply(strsplit(names(pairFst),"\\."),function(nms){c(from=nms[1],to=nms[2])})),stringsAsFactors=F))
        fsts = fsts[order(fsts$to,fsts$from),]

        dsts <- merge(dsts,fsts)

print("this far3")
        
        IBD <- lm(log(fst)~log(d),dsts)
        ibd.slope <- c(coef(IBD)[2])
        ibd.int <- c(coef(IBD)[1])
        bs <- segmentGLM(c(dsts$d),log(c(dsts$fst)))

	stats = c(localSNP, Smean, Ssd, privateSNP, prSmean, prSsd, total_priv, pwPriv, He, HeALL, pairFst,
                  ibd.slope=ibd.slope,ibd.int=ibd.int,bs.break=bs[1],bs.ll=bs[2],
                  na.lat.slope=na.lat.slope,
                  na.lat.int=na.lat.int, he.lat.slope=he.lat.slope, he.lat.int=he.lat.int)
        
	stats1 = matrix(data=stats, nrow = 1)
	colnames(stats1) = names(stats)
	stats = as.data.frame(stats1)
	stats
}



