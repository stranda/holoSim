#holoStats.R 
#Script to define a function to calculate stats for holoSim ABC simulations

holoStats.r = function(out, popDF, extent, cores=1) {
    
    allMAF = mafreq(out)
    totalHe = 2*allMAF*(1-allMAF)
    popid = strataNames(out)
    
    SNPs = sum(allMAF<1)
    names(SNPs) = "tot_SNPs"
    
    split_out = strataSplit(out) #list of strataG objects for each pop
###    locMAF = sapply(split_out,function(o){mafreq(o)})
    locMAF = do.call(cbind,mclapply(split_out,mc.cores=cores,function(o){mafreq(o)}))
    colnames(locMAF) <- names(split_out)
    locHe = colMeans(2*locMAF*(1-locMAF))
    varlocHe = apply(2*locMAF*(1-locMAF),2,var)
    
    locN = sapply(split_out,function(o){length(o@data$ids)})
    names(locN) <- popid

    localSNP = apply(locMAF,2,function(x){sum(x<1)})
    
    names(localSNP) = paste0("S.", popid)
    
    ##not sure we need this, does not seem to be variable and costs some time
    privateSNP = colSums(privateAlleles(out))
                                        #privateSNP = privateSNP[sample.order]
    names(privateSNP) = paste0("pS.", popid)

    total_priv = sum(privateSNP)
    names(total_priv) = "tot_priv"

    pwhet=pwise.het(locMAF,locN,cores)
    
    FstMat.loc = as.dist(pwise.fst.loc(locMAF,allMAF,locN,pwhet))
    neimat     = pwise.nei(locMAF,cores)
    
    pairnei = as.vector(neimat)
    neinames = c()
                                        #This will be different...
    for(pid in 1:(length(popDF$id)-1)) {
        neinames = c(neinames, paste0(popDF$id[pid],".", popDF$id[(pid+1):length(popDF$id)]))
    }
    
    names(pairnei) = neinames
    
    pairFst.loc = as.vector(FstMat.loc)
    Fstnames.loc = c()
                                        #This will be different...
    for(pid in 1:(length(popDF$id)-1)) {
        Fstnames.loc = c(Fstnames.loc, paste0(popDF$id[pid],".", popDF$id[(pid+1):length(popDF$id)]))
    }
    
    names(pairFst.loc) = Fstnames.loc
    
    tot_Fst = as.numeric(statFst(out)$result[1])
    names(tot_Fst) = "tot_Fst"

    eucdist = dist(t(locMAF))
    paireuc = as.vector(eucdist)
    eucnames = c()
    dnames = attr(eucdist,"Labels")
     for(pid in 1:(length(popDF$id)-1)) {
        eucnames = c(eucnames, paste0(popDF$id[pid],".", popDF$id[(pid+1):length(popDF$id)]))
    }
    names(paireuc) <- eucnames
########### some spatially focused stats
###     get pop coords:
        pops=t(matrix(1:prod(extent),ncol=extent[1]))
        popcrd=data.frame(t(sapply(popDF$grid.cell,function(i) {which(pops==i,arr.ind=T)})))
        names(popcrd)=c("y","x")
        popDF=cbind(popDF,popcrd)
        
###     get pairwise geo distances
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
        
        #sum_stats_gi<-summary(genind_out)
        #numAll=sum_stats_gi$pop.n.all
        #numAll=data.frame(id=names(numAll),na=numAll,stringsAsFactors=F)

#        he_by_pop<-as.vector(as.numeric(lapply(lapply(seppop(genind_out),summary),
#                                               function(x) mean(x$Hexp))))
    he_by_pop <- colMeans(2*locMAF*(1-locMAF))

    hedf <- data.frame(he=he_by_pop, id=unique(out@data$strata),stringsAsFactors=F)
    pr=prcomp(t(locMAF))
    pcadf <- data.frame(id=rownames(predict(pr)),pc1=predict(pr)[,1],pc2=predict(pr)[,2],pc3=predict(pr)[,3])
    
    popDF <- merge(pcadf,merge(popDF,hedf))

print("this far2")
    polyfit <- function(p,resp="he",ind="y",ord=1)
    {
        fit <- lm(p[,resp]~poly(p[,ind],ord))
        c(coef(fit))
    }

    he.lat.stats <- polyfit(popDF,"he","y",ord=2)
    names(he.lat.stats) <- paste0("helat.",c("int","frst","scnd"))
    he.long.stats <- polyfit(popDF,"he","x",ord=2)
    names(he.long.stats) <- paste0("helong.",c("int","frst","scnd"))
    pc1.lat.stats <- polyfit(popDF,"pc1","y",ord=2)
    pc1.long.stats <- polyfit(popDF,"pc1","x",ord=2)
    names(pc1.lat.stats) <- paste0("pc1lat.",c("int","frst","scnd"))
    names(pc1.long.stats) <- paste0("pc1long.",c("int","frst","scnd"))
    pc2.lat.stats <- polyfit(popDF,"pc2","y",ord=2)
    pc2.long.stats <- polyfit(popDF,"pc2","x",ord=2)
    names(pc2.lat.stats) <- paste0("pc2lat.",c("int","frst","scnd"))
    names(pc2.long.stats) <- paste0("pc2long.",c("int","frst","scnd"))

    pc3.lat.stats <- polyfit(popDF,"pc3","y",ord=2)
    pc3.long.stats <- polyfit(popDF,"pc3","x",ord=2)
    names(pc3.lat.stats) <- paste0("pc3lat.",c("int","frst","scnd"))
    names(pc3.long.stats) <- paste0("pc3long.",c("int","frst","scnd"))

    
    fsts = data.frame(fst=pairFst.loc)
    fsts = cbind(fsts,data.frame(t(sapply(strsplit(names(pairFst.loc),"\\."),function(nms){c(from=nms[1],to=nms[2])})),stringsAsFactors=F))
    fsts = fsts[order(fsts$to,fsts$from),]

    neis = data.frame(nei=pairnei)
    neis = cbind(neis,data.frame(t(sapply(strsplit(names(pairnei),"\\."),function(nms){c(from=nms[1],to=nms[2])})),stringsAsFactors=F))
    neis = neis[order(neis$to,neis$from),]



    
    edist = data.frame(edist=paireuc)
    edist = cbind(edist,data.frame(t(sapply(strsplit(names(paireuc),"\\."),function(nms){c(from=nms[1],to=nms[2])})),stringsAsFactors=F))
    edist = edist[order(edist$to,edist$from),]
    
    dsts <- merge(merge(merge(dsts,fsts),edist),neis)
    
    print("this far3")
        
    IBDfst <- lm(log(fst+1)~log(d),dsts)
    ibdfst.slope <- c(coef(IBDfst)[2])
    ibdfst.int <- c(coef(IBDfst)[1])
    bsfst <- segmentGLM(c(dsts$d),log(c(dsts$fst+1)))

    IBDnei <- lm(log(nei+1)~log(d),dsts)
    ibdnei.slope <- c(coef(IBDnei)[2])
    ibdnei.int <- c(coef(IBDnei)[1])
    bsnei <- segmentGLM(c(dsts$d),log(c(dsts$nei+1)))

       
    IBDedist <- lm(log(edist+1)~log(d),dsts)
    ibdedist.slope <- c(coef(IBDedist)[2])
    ibdedist.int <- c(coef(IBDedist)[1])
    bsedist <- segmentGLM(c(dsts$d),log(c(dsts$edist+1)))

    
    stats = c(SNPs, localSNP, privateSNP, total_priv, pairFst.loc, pairnei,tot_Fst,
              ibdfst.slope=ibdfst.slope,ibdfst.int=ibdfst.int,bsfst.break=bsfst[1],bsfst.ll=bsfst[2],
              ibdedist.slope=ibdedist.slope,ibdedist.int=ibdedist.int,bsedist.break=bsedist[1],bsedist.ll=bsedist[2],
              ibdnei.slope=ibdnei.slope,ibdnei.int=ibdnei.int,bsnei.break=bsnei[1],bsnei.ll=bsnei[2],
              
              he.lat.stats, he.long.stats,
              pc1.lat.stats, pc1.long.stats,
              pc2.lat.stats, pc2.long.stats,
              pc3.lat.stats, pc3.long.stats)
    
    stats1 = matrix(data=stats, nrow = 1)
    colnames(stats1) = names(stats)
    stats = as.data.frame(stats1)
    stats
}

holoStats <- cmpfun(holoStats.r)  #requires compiler


#This gets off when <347 variable loci in the file...
#Need this to work!
mask.data_OLD = function(out, popDF, nSNP, mask=dm) {

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

