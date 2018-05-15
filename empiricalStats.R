###
### takes calls puts them in a strataG gtypes and then runs holostats on it
###

fn = "calls_from_angsd_imputed_0.2_0.3.rda"  #R binary file with ash calls
xdim = 27 
ydim = 20

numloci=300  #if zero then take all loci.  Positive numbers are sampled without replacement, so they
           #need to be smaller or equal to the actual number of loci


load(fn)

library(strataG)
library(parallel)
library(compiler)
source("empty_or_sampled.R")
source("holoStats.R")
source("helper-functions.R")
source("segment-regression.R")
calls  = t(imputed)

id = as.numeric(gsub("fp","",c(t(matrix(rep(rownames(calls),2),ncol=2)))))
strata = sapply(id,function(i){popmap[popmap[,1]==i,2]})

if (numloci>0)
{
    calls=calls[,sort(sample(1:dim(calls)[2],numloci,FALSE))]
}

l = do.call(cbind,lapply(1:dim(calls)[2],function(g){
    do.call(c,strsplit(calls[,g],'')) }))
l <- as.data.frame(l)

names(l) <- colnames(calls)

sclose = as.data.frame(cbind(id=rownames(l),strata,l))
sclose = sclose[sclose$strata%in%popDF$id,]
sclosel = lapply(1:dim(sclose)[2],function (i) {as.character(sclose[,i])})
names(sclosel) = names(sclose)

out=df2gtypes(sclose,ploidy=1,strata.col=2,loc.col=3)
out@data$ids <- gsub("fp","",gsub(".$","",out@data$ids))
out@ploidy=2L  #snps are in columns correct number of chromosomes collected.
ashstats=holoStats(out,popDF=popDF,extent=c(ydim,xdim))
save(file=paste0("empirical_stats_",numloci,".rda"),ashstats)
