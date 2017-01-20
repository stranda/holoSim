#holoStats.R 
#Script to define a function to calculate stats for holoSim ABC simulations

holoStats = function(out, popDF) {
	require(hierfstat)
	
	popid = strataNames(out)
	SNPs = as.numeric(table(summarizeLoci(out)[,3]))
	names(SNPs) = "tot_SNPs"
	split_out = strataSplit(out)
	localHo = sapply(split_out, FUN=function(x){summary(x)$strata.smry})[5,]
	localSNP = sapply(split_out, FUN=function(x){length(which(summarizeLoci(x)[,3] == 2))})

	#sample.ids = unique(unlist(strsplit(colnames(localsum), split = " ")))[-1]
	#sample.order = order(as.numeric(sample.ids))

	#localSNP = localSNP[sample.order]
	names(localSNP) = paste0("S.", popid)
	#localHo = localHo[sample.order]
	names(localHo) = paste0("Ho.", popid)

	#localSNP = sapply(split_out, FUN=function(x){(summary(x)$strata.smry[3]*360)-360}) 
	#localHo = sapply(split_out, FUN=function(x){summary(x)$strata.smry[5]})

	total_Ho = mean(obsvdHet(out))
	names(total_Ho) = "tot_Ho"

	privateSNP = colSums(privateAlleles(out))
	#privateSNP = privateSNP[sample.order]
	names(privateSNP) = paste0("pS.", popid)
	
	total_priv = sum(privateSNP)
	names(total_priv) = "tot_priv"
	
	genind_out = gtypes2genind(out)
	FstMat = pairwise.fst(genind_out, pop = genind_out$pop)
	pairFst = as.vector(FstMat)
	Fstnames = c()
	
	#This will be different...
	for(pid in 1:(length(popDF$id)-1)) {
		Fstnames = c(Fstnames, paste0(popDF$id[pid],".", popDF$id[(pid+1):length(popDF$id)]))
	}
	names(pairFst) = Fstnames
	
	tot_Fst = as.numeric(statFst(out)$result[1])
	names(tot_Fst) = "tot_Fst"
	
	stats = c(SNPs, localSNP, localHo, total_Ho, privateSNP, total_priv, pairFst, tot_Fst)
	stats1 = matrix(data=stats, nrow = 1)
	colnames(stats1) = names(stats)
	stats = as.data.frame(stats1)
	stats
}



get.nSNP = function(out) {

	polysites = which(summarizeLoci(out)[,3] == 2)
	nSNPloci = length(polysites)
	nSNPloci
}


#This gets off when <347 variable loci in the file...
#Need this to work!
mask.data = function(out, popDF, nSNP, mask=dm) {

	require(hierfstat)

	#Rename populations if dealing with simulated data...
	strat_id = c()
	for(popu in 1:length(popDF$id)) {
		if(popu == 1) {
			strat_id = rep(as.character(popDF$id[popu]), popDF$sample.size[popu])
		} else {
			strat_id = append(strat_id, rep(as.character(popDF$id[popu]), popDF$sample.size[popu]))
		}
	}

	strata(out) <- strat_id	

	tmp = as.matrix(out)
	all1 = seq(3,dim(tmp)[2],2)
	tmp2 = tmp[,1:2]
	keeploc = c()
	varloc = c()
	sampleme = c()
	for(loc in 1:nSNP) {
		locmat = tmp
		NArows = which(mask[,all1[loc]] == TRUE)
		locmat[NArows,] = NA
		sites = apply(locmat[-NArows,], 2, FUN=function(x){length(unique(x))})[all1]
		
		if(loc == 1) {
			varloc = which(sites == 2)
			#sampleme = varloc
		} else {
			sites[keeploc] = 0
			varloc = which(sites == 2)
			#sampleme = varloc[!(varloc %in% keeploc)]
		}
		

		if(length(varloc) > 0) {
			if(length(varloc) == 1) {
				keeploc[loc] = varloc
			} else {
				keeploc[loc] = sample(varloc,1,replace=FALSE)
			}
			tmp2 = cbind(tmp2,locmat[,((2*keeploc[loc])+1):((2*keeploc[loc])+2)])
			#tmp2[NArows,((2*loc)+1):((2*loc)+2)] = NA
		} 
	}

	#tmp2[mask == TRUE] = NA
	tmp2 = as.data.frame(tmp2)
	out2 = df2gtypes(tmp2, ploidy = 2)

	out2
} 

