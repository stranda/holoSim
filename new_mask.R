mask.data = function(tmpout, popDF, nSNP, mask) {

#	This isn't particularly fast... lots of as.data.frame calls
#	strat_id will be the new strata IDs, looping over samples
#	requires popDF$id to be in the same order as in the FSC output
#	that order follows the grid.cell (which is the sorting vector for popDF currently)
	strat_id = c()
	for(pop in 1:length(popDF$id)) {
		if(pop == 1) {
			to.mask = as.data.frame(tmpout[,,paste("Sample",pop)])
			strat_id = rep(as.character(popDF$id[pop]), popDF$sample.size[pop])
		} else {
			to.mask = rbind(to.mask, as.data.frame(tmpout[,,paste("Sample",pop)]))
			strat_id = append(strat_id, rep(as.character(popDF$id[pop]), popDF$sample.size[pop]))
		}
	}

#	IDvec will be individual IDs.  Needed something that wouldn't lead to odd sorting
#	e.g., 1, 10, 11, 12, 13, 14, 15, 2, 3, 4, 5, 6, 7, 8, 9 
	IDvec = c(1001:(1000+length(to.mask$ids)))
	to.mask$strata = strat_id
	to.mask$ids = paste0("ind_",IDvec,"_",to.mask$strata)

#	Turn this on and run up to this point for a check on whether pops were being shuffled
#	named_out = df2gtypes(to.mask,ploidy=2,id.col=1,strata.col=2)  

#	all1 gives the column number where the first allele of all loci sits
	all1 = seq(3,dim(to.mask)[2],2)
	
#	First attempt at masking used the first nSNP loci (two columns per locus in to.mask)
#	However, with the infinite sites model we want to break up potential associations between SNPs on the same chrom
#	version 2 samples to form the initial dataset
	starters = sample(all1, nSNP, replace = FALSE)
	masked = to.mask[,c(1:2,sort(c(starters, starters+1)))]
#	masked = to.mask[,1:(2*nSNP+2)]
#	masked and mask should have the same dimensions, strata order 
	masked[mask==TRUE] = NA



#	onecoldata is a one-column genetic dataset, only used to check for polymorphism following mask
	onecoldata = matrix(data = NA, nrow = 2*dim(masked)[1], ncol = nSNP)
	for(x in 1:nSNP) {
		onecoldata[,x] = c(masked[,all1[x]], masked[,all1[x]+1])
	}

#	loci that have only one unique state (excluding NAs) after masking need to be replaced
	fixme = apply(onecoldata, 2, FUN=function(x){length(unique(x[is.na(x)==FALSE]))})
	to.fix = which(fixme == 1)

# 	inuse gives the columns for the loci (both alleles) that are currently in the dataset
	inuse = sort(c(starters, starters+1))

#	looping over the list of loci that need replacing
	if(length(to.fix) > 0) {
		for(fl in 1:length(to.fix)) {

#			to.mask has all data, start from there, but turn inuse columns to NA
			remaining = to.mask
			remaining[,inuse] = NA

#			which column in the dataset needs to be replaced this time?
			mat.col.leaving = all1[to.fix[fl]]  
#			grab the masking pattern for that particular locus (it's the same for both alleles)
			currmask = mask[,all1[to.fix[fl]]]  

#			mask ALL available loci using this pattern
			remaining[currmask==TRUE,] = NA

#			check for polymoprhism in the masked loci
			onecoldata = matrix(data = NA, nrow = 2*dim(remaining)[1], ncol = length(all1))
			for(x in 1:length(all1)) {
				onecoldata[,x] = c(remaining[,all1[x]], remaining[,all1[x]+1])
			}

#			choices will be a list of the number of alleles in each locus
			choices = apply(onecoldata, 2, FUN=function(x){length(unique(x[is.na(x)==FALSE]))})

#			If no available loci contain polymorphisms, have to try again (back to FSC)
			if(sum(choices == 2) > 0) {
				if(sum(choices == 2) == 1) {
#					Don't use sample if there's only 1 possible value, returns random number	
					chose = which(choices == 2)
				} else {
#					If there are multiple options, sample those options
					chose = sample(which(choices==2), 1, replace = FALSE)
				}

#				Change the appropiate column (that needed to be fixed) to the masked data for the chosen available locus
				masked[,mat.col.leaving] = remaining[,all1[chose]]
				masked[,mat.col.leaving+1] = remaining[,all1[chose]+1]

#				Add the chosen locus' column to the list of inuse loci
				inuse = c(inuse, all1[chose], all1[chose]+1)

			} else {
				chose = c()
				return("BUMMER! too few SNPs due to missing data")
				stop("BUMMER! too few SNPs due to missing data")
			}
		}
	}

#	Convert data back to gtypes object, specify ploidy, id.col, strata.col
	postmask = df2gtypes(masked, ploidy = 2, id.col = 1, strata.col = 2)

	postmask
}


