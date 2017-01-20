#runFSC.R
#Script to run Fastsimcoal using strataG

runFSC = function(pops = NULL, rland = NULL, parms=NULL, sample_pops = NULL, sample_n = NULL, label = "test", marker = "snp", nloci = 500, delete.files = TRUE, num.cores = 1, exec="fsc25", growth.model = "step") {
	

	#ID the empty pops
	empty_pops = pops$pop[is.na(pops$arrive) == TRUE]

#POP.INFO 
	#Definitions for pop.info
	h = parms$xdim*parms$ydim   
	pop.size = rep(parms$Ne,h)
	sample.size = rep(0,h)
	sample.size[sample_pops] = sample_n
	
	#Remove populations that are never founded or that are outside the species range
	pop.size=pop.size[-empty_pops]
	sample.size = sample.size[-empty_pops]

	#Calculate growth rates -- use step change or exponential, both work (!?!)
	if(growth.model == "step") {
		growth.rate = rep(0,h)
		growth.rate = growth.rate[-empty_pops]
	} else if(growth.model == "exp") {
		init_size = rep(parms$found_Ne,h)
		init_size[pops$arrive == 0] = parms$ref_Ne
		init_size = init_size[-empty_pops]
		growth.rate = log(init_size/parms$Ne)/(parms$texp-pops$arrive[-empty_pops])
	}

	pop.info = fscPopInfo(pop.size = pop.size, sample.size = sample.size, growth.rate = growth.rate)	

#LOCUS.PARAMS 
	locus.params = fscLocusParams(locus.type = marker, num.loci = nloci, mut.rate = parms$seq.length*parms$mu)

#MIGRATION MATRICES
	full_IDs = pops$pop[-empty_pops]
	full_pops = pops[-empty_pops,]
	full_pops$pop = c(1:length(full_pops$pop))-1

	for(p in 1:length(full_IDs)) {
		if(full_pops$arrive[p] > 0) {
			full_pops$source[p] = which(full_IDs == full_pops$source[p])
			full_pops$source[p] = full_pops$source[p]-1
		}
	}
	sorted.pops = full_pops[order(full_pops$arrive, decreasing = TRUE),]
	sorted.events = sorted.pops[sorted.pops$arrive > 0,]

	oddrow = seq(2,dim(rland$demography$epoch[[1]]$M)[1],by=2)
	migmat = vector("list", (length(unique(sorted.events$arrive))+1))
	migmat[[1]] = rland$demography$epoch[[1]]$M[oddrow,oddrow]
	migmat[[1]] = migmat[[1]][-empty_pops,]
	migmat[[1]] = migmat[[1]][,-empty_pops]

	#Instead of looping over populations, loop over unique time points in the sorted.events object?
	sorted.events$migmat = "nomig"
	sorted.events$botmigmat = "nomig"
	leaving = vector("list", length(unique(sorted.events$arrive)))
	for(m in 1:length(unique(sorted.events$arrive))) {
		leaving[[m]] = sorted.events$pop[sorted.events$arrive == unique(sorted.events$arrive)[m]]
		leaving[[m]] = leaving[[m]]+1
		migmat[[m+1]] = migmat[[m]]
		migmat[[m+1]][leaving[[m]],] = 0
		migmat[[m+1]][,leaving[[m]]] = 0
		sorted.events$migmat[sorted.events$arrive == unique(sorted.events$arrive)[m]] = m
		sorted.events$botmigmat[sorted.events$arrive == unique(sorted.events$arrive)[m]] = m-1
	}

	#Check the migration matrix -- leaving this here in case 
	#zerorows = vector("list", length(migmat)) 
	#newgone = vector("list", length(migmat))
	#zerorows[[1]] = which(rowSums(migmat[[1]]) == 0)
	#for(l in 2:length(migmat)) {
	#	zerorows[[l]] = which(rowSums(migmat[[l]]) == 0)
	#	newgone[[l]] = zerorows[[l]][which(zerorows[[l]] %in% zerorows[[l-1]] == FALSE)]
	#}
	#unlist(newgone)
	#unlist(leaving)
	#unlist(newgone)[-250] - unlist(leaving)

#HISTORICAL EVENTS 
	hist.ev = fscHistEv(num.gen = parms$texp-sorted.events$arrive, source.deme = sorted.events$pop, sink.deme = sorted.events$source, prop.migrants = 1, new.sink.size = 1, new.sink.growth = "keep", new.mig.mat = sorted.events$migmat)
	hist.ev1 = fscHistEv(num.gen = parms$texp, source.deme = full_pops$pop[full_pops$arrive == 0], sink.deme = full_pops$pop[full_pops$arrive == 0], prop.migrants = 0, new.sink.size = 1, new.sink.growth = 0, new.mig.mat = "nomig")
	hist.ev2 = fscHistEv(num.gen = parms$preLGM_t[1], source.deme = full_pops$pop[full_pops$arrive == 0], sink.deme = full_pops$pop[full_pops$arrive == 0], prop.migrants = 0, new.sink.size = round(parms$preLGM_Ne/parms$ref_Ne,5), new.sink.growth = 0, new.mig.mat = "nomig")
	hist.ev = rbind(hist.ev, hist.ev1, hist.ev2)
	hist.ev[length(sorted.events[,1]),7] = "nomig"

	#Step change and exponential both require extra lines
	#Step change extra line is a bottleneck 1 generation before (after) population founding, looking backwards (forward) in time
	#Exponential extra line turns population decline off before population fusion, looking backward in time (fixes TDemeCollection error)
	if(growth.model == "step") {
		bott.ev = fscHistEv(num.gen = parms$texp-sorted.events$arrive-1, source.deme = sorted.events$pop, sink.deme = sorted.events$pop, prop.migrants = 0, new.sink.size = round(parms$found_Ne/parms$Ne,5), new.sink.growth = 0, new.mig.mat = sorted.events$botmigmat)
		bott.ev[bott.ev[,1] < 0,1] = 0
		hist.ev = rbind(hist.ev, bott.ev)
		hist.ev = hist.ev[order(as.numeric(hist.ev[,1]), as.numeric(hist.ev[,4]), decreasing = FALSE),]
	} else if(growth.model == "exp") {
		stop.grow = fscHistEv(num.gen=parms$texp-sorted.events$arrive, source.deme = sorted.events$pop, sink.deme = sorted.events$pop, prop.migrants = 1, new.sink.size = 1, new.sink.growth = 0, new.mig.mat = sorted.events$botmigmat)
		hist.ev = rbind(hist.ev, stop.grow)
		hist.ev = hist.ev[order(as.numeric(hist.ev[,1]), as.numeric(hist.ev[,6]), decreasing = FALSE),]
	}


	#out = fastsimcoal(label = label, pop.info = pop.info, locus.params = locus.params, mig.rates = migmat, hist.ev = hist.ev)
	out = fastsimcoal(label = label, pop.info = pop.info, locus.params = locus.params, mig.rates = migmat, hist.ev = hist.ev, num.cores = num.cores, delete.files = delete.files, exec = exec)

	out
}


