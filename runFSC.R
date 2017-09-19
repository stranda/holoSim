#runFSC.R
#Script to run Fastsimcoal using strataG

runFSC = function(pops = NULL, rland = NULL, parms=NULL, sample_pops = NULL, sample_n = NULL, label = "test", marker = "dna", nloci = 500, delete.files = TRUE, num.cores = 1, exec="fsc25", growth.model = "step") {
	

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
	attr(fscmodel[[3]], "ploidy") = 2
	attr(fscmodel[[3]], "opts") = "-I -s"
	
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
	out = myfsc(label = label, pop.info = pop.info, locus.params = locus.params, mig.rates = migmat, hist.ev = hist.ev, num.cores = num.cores, delete.files = delete.files, exec = exec)

	out
}

#myfsc()
#Changed from strataG's fastsimcoal() fxn (turning off -S option, only output polymorphic sites)
myfsc = function (pop.info, locus.params, mig.rates = NULL, hist.ev = NULL, 
    label = NULL, quiet = TRUE, delete.files = TRUE, exec = "fsc252", 
    num.cores = NULL, label.haplotypes = TRUE) 
{
    if (is.null(label)) 
        label <- "fsc.run"
    label <- make.names(label)
    if (file.exists(label)) 
        for (f in dir(label, full.names = T)) file.remove(f)
    if (!quiet) 
        cat("fastsimcoal: writing input file\n")
    infile <- fscWrite(pop.info = pop.info, locus.params = locus.params, 
        mig.rates = mig.rates, hist.ev = hist.ev, label = label)
    if (!quiet) 
        cat("fastsimcoal: running\n")
    cores.spec <- if (!is.null(num.cores)) {
        num.cores <- max(1, num.cores)
        num.cores <- min(num.cores, min(detectCores(), 12))
        if (num.cores == 1) 
            ""
        else paste(c("-c", "-B"), num.cores, collapse = " ")
    }
    else ""
    cmd.line <- paste(exec, "-i", infile, "-n 1", ifelse(quiet, 
        "-q", ""), cores.spec, attr(locus.params, "opts"))
    err <- if (.Platform$OS.type == "unix") {
        system(cmd.line, intern = F)
    }
    else {
        shell(cmd.line, intern = F)
    }
    if (err == 0) {
        if (!quiet) 
            cat("fastsimcoal exited normally\n")
    }
    else {
        stop("fastsimcoal exited with error ", err, "\n")
    }
    arp.file <- file.path(label, paste(label, "_1_1.arp", sep = ""))
    if (!file.exists(arp.file)) 
        stop("fastsimcoal did not generate output")
    if (!quiet) 
        cat("fastsimcoal: parsing output to gtypes\n")
    g <- myfscRead(arp.file, locus.params, label.haplotypes)
    if (delete.files) {
        if (!quiet) 
            cat("fastsimcoal: removing output files\n")
        unlink(label, recursive = TRUE, force = TRUE)
        file.remove(infile)
        file.remove("seed.txt")
    }
    return(g)
}

#myfscRead()
#Changed from fscRead() function in strataG
#This version greps the number of polymorphic sites from the arlequin file
#Allows simulation of infinite sites model (i.e., strataG no longer expects a dataset of num.chrom * num.markers bp)
myfscRead = function (file, locus.params, label.haplotypes = FALSE) 
{
    .formatGenotypes <- function(x, ploidy) {
    	require(swfscMisc)
        nloci <- ncol(x) - 2
        loc.end <- seq(ploidy, nrow(x), by = ploidy)
        gen.data <- do.call(rbind, lapply(loc.end, function(i) {
            allele.i <- (i - ploidy + 1):i
            loci <- as.vector(x[allele.i, -(1:2)])
            id <- paste(x[allele.i, 2], collapse = ".")
            pop <- x[allele.i[1], 1]
            c(id, pop, loci)
        }))
        locus_names <- paste("Locus", zero.pad(1:nloci), sep = "_")
        locus_names <- paste(rep(locus_names, each = ploidy), 
            1:ploidy, sep = ".")
        colnames(gen.data) <- c("id", "pop", locus_names)
        gen.data
    }

    .formatDNA <- function(dna.seq, pop, locus.params, label) {
        #arp.file <- file.path(label, paste(label, "_1_1.arp", sep = ""))
        #num.chrom <- attr(locus.params, "num.chrom")
        num.chrom = system(paste("grep 'Total number of polymorphic sites:'", file, "| cut -f 2 -d : | cut -f 2 -d ' '"), intern = TRUE)
        chrom.pos <- if (is.null(num.chrom)) {
            tapply(locus.params$num.markers, locus.params$chromosome, 
                sum)
        }
        else {
            rep(sum(1), num.chrom)
        }
        chrom.pos <- cumsum(chrom.pos)
        chrom.pos <- cbind(start = c(1, chrom.pos[-length(chrom.pos)] + 
            1), end = chrom.pos)
        rownames(dna.seq) <- pop
        dna.seq <- tolower(dna.seq)
        new("multidna", lapply(1:nrow(chrom.pos), function(i) {
            as.matrix(dna.seq)[, chrom.pos[i, "start"]:chrom.pos[i, 
                "end"]]
        }))
    }
    #print("reading file")
    f <- readLines(file)
    start <- grep("SampleData=", f) + 1
    end <- which(f == "}") - 2
    pos <- cbind(start, end)
    .compileMatrix <- function(i, pos) {
        f.line <- f[pos[i, 1]:pos[i, 2]]
        f.line <- gsub("[[:space:]]+", "--", f.line)
        result <- do.call(rbind, strsplit(f.line, "--"))[, -2]
        cbind(rep(paste("Sample", i), nrow(result)), result)
    }
    #print("compiling matrix")
    data.mat <- do.call(rbind, lapply(1:nrow(pos), .compileMatrix, 
        pos = pos))
    ploidy <- attr(locus.params, "ploidy")
    data.type <- f[grep("DataType=", f)]
    data.type <- gsub("\tDataType=", "", data.type)
    switch(data.type, DNA = {
    	#print("splitting strings")
        dna.seq <- do.call(rbind, strsplit(data.mat[, 3], ""))
        if (attr(locus.params, "ploidy") == 2) {
            gen.data <- .formatGenotypes(cbind(data.mat[, 1:2], 
                dna.seq), ploidy)
            df2gtypes(gen.data, ploidy, description = file)
        } else {
        	#print("formatting DNA")
            dna.seq <- .formatDNA(dna.seq, data.mat[, 2], locus.params)
            #print("seq2gtype - SOMETHING IS WRONG HERE!!!  Has to do with creation of gen.data")
            g <- sequence2gtypes(dna.seq, strata = data.mat[,1], description = file)
            if (label.haplotypes) labelHaplotypes(g)$gtype else g
        }
    }, MICROSAT = {
        gen.data <- .formatGenotypes(data.mat, ploidy)
        df2gtypes(gen.data, ploidy, description = file)
    }, NULL)
}
