#holoSim-Ash.R
#Script to run simulations for the holoSim project, Ash dataset

################################
#Set number of replicates and node (CPU)
args = commandArgs(TRUE)
node = as.numeric(args[1])   #This is the computer being used
nrep = as.numeric(args[2])   #This is the number of replicates run
#use -
# Rscript holoSim_Ash.R 1 100 
# to run 100 replicates on the first CPU
# output will be Ash_<node>.out

#node = 0
#nrep = 10
#Rsource = "~/Google Drive/MSU/HOLOSIM/Software/SimTemplate/github/"
#SIMS="~/Desktop/AshSims"
Rsource = "~/GoogleDrive/src/holoSim/"
SIMS="~/tmp/AshSims"
#Rsource = "/mnt/research/ABC/holoSim/Rsource/"
#SIMS = "/mnt/research/ABC/holoSim/Ash/"

setwd(Rsource)
source("colrate.R") #Fxns to calculate colonization rate measures
source("runFSC.R")  #Uses strataG to run fastsimcoal, given pops, rland, parms, etc.
source("empty_or_sampled.R") #Vectors of empty grid cells, sampled grid cells, sample sizes
source("holoStats.R") #Functions to calculate stats from simulations, mask data, count SNPs
#Load stuff from Allan's fastsim folder
source("make-landscape2.R") #changed so pops grow, fecundity increased (1.05 -> 1.2), also changed popcoords fxn to fix numbering
source("getpophist3.R") #changed to track population abundance at beginning and end of window specified for getC4(), also new logical test!
source("plothist.R")  #Unchanged
source("segment-regression.R")
source("fast_mask.R") #New data masking step, shorter loop than previous version


dm = read.csv("data_mask.csv", header = TRUE)

setwd(SIMS)

#Load libraries
library(rmetasim)
library(strataG)
library(hierfstat)
library(parallel)

#Set fixed aspects of simulations
time_convert = 30 #Assumes 30 years per time click (=1 generation)
col_window = c(0,round((2000/time_convert)))
dens.scale = 0.05
xdim = 27 
ydim = 20
marker = "snp"
nloci = 750 
nSNP = 347 
seq.length = 80  #Length of individual RAD contigs -- scales mutation rate for fscLocusParams()

#I think these may be important... may require some thought / adjustment
#Won't total fecundity per generation (=per time click) be important for the pace of range expansion??
local_N = 20  #These aren't drawn as parameters -- coal_Ne will be though
ref_N = 20    #Not a parameter?

simout = vector("list",nrep)

repl = 1
while(repl <= nrep) {
	
	#Draw forward sim parameter values from (eventually hyper)priors 
	texp = round(runif(1, (10000/time_convert),(20000/time_convert)))
	refs = sample(c(148,148,148),1,replace = TRUE)   #Right now only allowing a refuge in Southeast Texas
	refconfig = paste(refs,sep="", collapse=".")
	
	mix = runif(1, 0.001, 0.1)  #0.1% to 10% LDD?
	longmean = runif(1, 1, 2)  #1-2 cell widths 
	shortscale = 0.1 #This is a fixed parameter now, 0.1 

	#Draw coalescent sim parameter values from (eventually hyper)priors
	#These Ne's are Ne (number of diploids), same with sample sizes!! 
	Ne = round(runif(1, 500, 5000))
	found_Ne = round(runif(1, 10, 50))
	ref_Ne = round(runif(1, 500, Ne))
	preLGM_Ne = round(runif(1,10000,100000))
	preLGM_t = round(runif(1,(100000/time_convert),(150000/time_convert)))
	mu = runif(1,1e-8, 1e-7)

	#Bind all parameters into a dataframe
	parms = data.frame(node, repl, time_convert, col_window[1], col_window[2], dens.scale, xdim, ydim, seq.length, nloci, nSNP, seq.length, local_N, ref_N, texp, refconfig, mix, longmean, shortscale, Ne, found_Ne, ref_Ne, preLGM_Ne, preLGM_t, mu)
	#Set up the rmetasim landscape based on these parameters
	h = parms$xdim*parms$ydim
	local_k = rep(local_N, h)
	local_k[empty_cells] = 0

	l = recolonizeLandscape(dens.scale = dens.scale,
                            h=h,
                            xdim=parms$xdim,
                            ydim=parms$ydim,
                            refs=refs,
                            k=local_k,
                            sizeref=ref_N,
                            mix=mix,
                            longmean=longmean,
                            shortscale=shortscale
                            )

	missing_sampled = length(samp_pops)
	pops = c()

	#May need some kind of kill switch for this loop...
	forward=0
	while(missing_sampled > 0 & forward < 10) {
		forward = forward+1
		while(is.data.frame(pops) == FALSE) {
			pops = try(getpophist(l, maxtime=texp, window = col_window))
		}
		missing_sampled = length(which(is.na(pops$arrive[samp_pops]) == TRUE))
	}

#This is the first big check point.  Possible that parameters don't allow colonization of the whole landscape
#If that happens, move on to the next replicate
#Allows 10 tries to fill all populations that we've sampled
	if(missing_sampled == 0) {
	#Calculate movement rate over the specified window
		parms$C1 = getC(pops, window = col_window)
		parms$C1n = getCn(pops, window = col_window)
		parms$C1s = getCs(pops, window = col_window)
	
		parms$C2 = getC2(pops, window = col_window)
	
		parms$C3 = getC3(pops, window = col_window)
		parms$C3n = getC3n(pops, window = col_window)
		parms$C3s = getC3s(pops, window = col_window)
	
		parms$C4 = getC4(pops, window = col_window)
		parms$C4n = getC4n(pops, window = col_window)
		parms$C4s = getC4s(pops, window = col_window)

	#Run the simulation
		SNPloci = 0
		FSCtries = 0
		while(SNPloci < nSNP & FSCtries < 5) {
			FSCtries = FSCtries+1
			out = runFSC(pops=pops, rland = l, parms =parms, sample_pops = samp_pops, sample_n = sampn, label = paste0("Ash_", node, "-", repl), marker = marker, nloci = nloci, delete.files = TRUE, num.cores = 1, exec="fsc25", growth.model = "step")
			out2 = mask.data2(out, popDF, nSNP, mask = dm)
			SNPloci = get.nSNP(out2)
		}

#This is the second checkpoint.  If you can't make a dataset of 347 SNPs with the missing data structure, move on
#Allows 5 tries with fastsimcoal
	#If the right number of SNPs come out of the simulation, write a line, otherwise replicate is skipped
		if(SNPloci == nSNP) {
			
                    print("about to run holostats")
			stats = holoStats(out2, popDF, extent=c(ydim,xdim))
                    print("done running holostats")
                    
###risky to wait till end of loop to save, but makes for really convenient data structure
###I propose not running too many reps per node, but using lots of nodes?
###still saving the original way for stats output.
###                     save(parms, out2, file = paste0("Ash_", node,"-",repl,".Rda"))
                    simout[[repl]] <- list(out=out2, parms=parms, popDF=popDF, stats=stats)

			done = format(Sys.time(), "%H:%M:%S")   #Better time?  t2-t1 would need a tweak to always measure on the same scale
			#Print all parms and stats to a file
			if(!file.exists(paste0("Ash_", node,".out"))) {
				write.table(cbind(parms,stats,done), paste0("Ash_", node,".out"), quote = FALSE, row.names = FALSE)
			} else {
				write.table(cbind(parms,stats,done), paste0("Ash_", node,".out"), col.names = FALSE, quote = FALSE, row.names = FALSE, append = TRUE)
			}

			#Advance the replicate counter
			repl = repl+1

		}
	}
}

save(file=paste0("Ash_",node,"_",gsub(" ","",date()),".rda"),simout)

