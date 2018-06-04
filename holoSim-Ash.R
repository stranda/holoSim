#holoSim-Ash.R
#Script to run simulations for the holoSim project, Ash dataset

################################
#Set number of replicates and node (CPU)
args = commandArgs(TRUE)
#if (length(args)>0)
if (FALSE)
{
node = as.numeric(args[1])   #This is the computer being used
nrep = as.numeric(args[2])   #This is the number of replicates run
#use -
# Rscript holoSim_Ash.R 1 100 
# to run 100 replicates on the first CPU
# output will be Ash_<node>.out
} else {
   node = 2
   nrep = 3
}

Rsource = "."
SIMS=paste0("/var/tmp/Ash")

#Rsource = "~/GoogleDrive/src/holoSim/"
#SIMS="~/tmp/AshSims"
#Rsource = "/mnt/research/ABC/holoSim/Rsource/"
#SIMS = "/mnt/research/ABC/holoSim/Ash/"

#setwd(Rsource)


#Load libraries
library(rmetasim)
library(strataG)
library(compiler)
library(sys)
library(parallel)


source("colrate.R") #Fxns to calculate colonization rate measures
source("runFSC.R")  #Uses strataG to run fastsimcoal, given pops, rland, parms, etc.
source("empty_or_sampled.R") #Vectors of empty grid cells, sampled grid cells, sample sizes
source("holoStats.R") #Functions to calculate stats from simulations, mask data, count SNPs
#Load stuff from Allan's fastsim folder
source("make-landscape2.R") #changed so pops grow, fecundity increased (1.05 -> 1.2), also changed popcoords fxn to fix numbering
source("landscape-functions.R")
source("getpophist4.R") #changed to track population abundance at beginning and end of window specified for getC4(), also new logical test!
source("plothist.R")  #Unchanged
source("segment-regression.R")
#source("new_mask.R")
source("helper-functions.R")

#if(!dir.exists(SIMS)) {dir.create(SIMS)}
#setwd(SIMS)

cores = 1

#enableJIT(3)

###FastSimcoal can 'run away' under certain conditions. We set the next parameter to the time we wait
###in seconds before giving up on these params
fsctimeout = 690

####Set fixed aspects of simulations
time_convert = 30 #Assumes 30 years per time click (=1 generation)
col_window = c(0,round((2000/time_convert)))
dens.scale = 0.05
xdim = 27 
ydim = 20
marker = "snp"
nloci = 3350
nSNP = 2909
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
	refs = sample(c(9,148,378,98),1,prob=c(0.33,0,0.33,0.34))   # 9=texas, 148=TN, 378=NB. 98=SC
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
	mu = runif(1,1e-7, 1e-6)

	#Bind all parameters into a dataframe
	parms = data.frame(node, repl, time_convert, col_window[1], col_window[2], dens.scale, xdim, ydim, seq.length, nloci, nSNP, local_N, ref_N, texp, refconfig, mix, longmean, shortscale, Ne, found_Ne, ref_Ne, preLGM_Ne, preLGM_t, mu)
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
    print (date())
	missing_sampled = length(samp_pops)
	pops = c()


    print(parms)
    
	#May need some kind of kill switch for this loop...
	forward=0
	while(missing_sampled > 0 & forward < 10) {
		forward = forward+1
		while(is.data.frame(pops) == FALSE) {
                    print("trying getpophist")
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
            print(date())
		SNPloci = 0
		FSCtries = 0
		while((SNPloci < nSNP) & (FSCtries < 2)) {
                    FSCtries = FSCtries+1
                    
                    out=NULL
                    out = tryCatch(
                    {
                        eval_safe( #allows timeout
                        {
                            runFSC(pops=pops, rland = l, parms =parms, sample_pops = samp_pops, sample_n = sampn, label = paste0("Ash_", node, "-", repl), marker = marker, nloci = nloci, delete.files = T, num.cores = cores, exec="fsc251", growth.model = "step")
                        },timeout=fsctimeout)
                    },
                                   warning=function(w){print(w)},
                                   error=function(e){print(e); NULL})

                    if (!is.null(out))
                    {
                        out2 = out
                        SNPsum = get.gSum(out2)
                        SNPloci = SNPsum$nvar
                    } else {
                        print("fsc timed out")
                        SNPsum=list(NULL,NULL)
                        SNPloci = 0
                    }
                    }
                    
                    print(date())
            
#This is the second checkpoint.  If you can't make a dataset of 2911 SNPs with the missing data structure, move on
#Allows 5 tries with fastsimcoal
	#If the right number of SNPs come out of the simulation, write a line, otherwise replicate is skipped
            if (SNPloci >= nSNP)
            {
                ## assign proper ids to FSC output
               
                out2=fsc.rename(out=out2,popDF=popDF,dip=TRUE)

                if (SNPloci>nSNP) out2=out2[,which(SNPsum$nall>1)[1:nSNP],]  #subsets to just the loci needed

#                save(file="tmpout.rda",out2)
                print("about to run holostats")
                stats = holoStats(out2, popDF, extent=c(ydim,xdim),cores=cores)
                print("done running holostats")

                done = format(Sys.time(), "%H:%M:%S")   #Better time?  t2-t1 would need a tweak to always measure on the same scale
                simout[[repl]] <- cbind(parms,stats,done)
                repl = repl+1

            }
	}
}

simoutdf <- do.call(rbind, simout)

write.table(file=paste0("Ash_",floor(runif(1,0,100000)),"_",gsub(" ","",date()),".csv"),sep=",",row.names=F,simoutdf)



print(date())


