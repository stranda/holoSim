source("colrate.R") #Fxns to calculate colonization rate measures
source("runFSC.R")  #Uses strataG to run fastsimcoal, given pops, rland, parms, etc.
source("empty_or_sampled.R") #Vectors of empty grid cells, sampled grid cells, sample sizes
source("holoStats.R") #Functions to calculate stats from simulations, mask data, count SNPs
#Load stuff from Allan's fastsim folder
source("make-landscape2.R") #changed so pops grow, fecundity increased (1.05 -> 1.2), also changed popcoords fxn to fix numbering
source("getpophist4.R") #changed to track population abundance at beginning and end of window specified for getC4(), also new logical test!
source("plothist.R")  #Unchanged
source("segment-regression.R")
source("new_mask.R")
source("landscape-functions.R")

library(rmetasim)
library(strataG)
library(hierfstat)
library(parallel)
library(dplyr)

node=1
repl=1
time_convert=30
refs = 2
refconfig = paste(refs,sep="", collapse=".")

mix = runif(1, 0.001, 0.1)  #0.1% to 10% LDD?
longmean = runif(1, 1, 3)  #1-2 cell widths 
shortscale = 0.1 #This is a fixed parameter now, 0.1 

                                        #Draw coalescent sim parameter values from (eventually hyper)priors
                                        #These Ne's are Ne (number of diploids), same with sample sizes!! 
Ne = round(runif(1, 500, 5000))
found_Ne = round(runif(1, 10, 50))
ref_Ne = round(runif(1, 500, Ne))
preLGM_Ne = round(runif(1,10000,100000))
preLGM_t = round(runif(1,(100000/time_convert),(150000/time_convert)))
mu = runif(1,1e-8, 1e-7)
seq.length=80
nloci=2911
dens.scale = 0.05
xdim = 6
ydim = 6
local_N=200
ref_N=1000
                                        #Bind all parameters into a dataframe
parms = data.frame(node, repl, time_convert, dens.scale, xdim, ydim, seq.length, nloci)
                                        #Set up the rmetasim landscape based on these parameters
h = parms$xdim*parms$ydim
local_k = rep(local_N, h)
#

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

ph = getpophist(l,maxtime=100,window=c(0,25))

