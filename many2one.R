###
###
### convert a directory of files 
###
###
library(parallel)
many = "statsout"
one = "processed/ash_sum_stats.csv"
archive="processed"

cores=4

if (file.exists(one)) onedf <- read.csv(one) else onedf <- NULL

files = list.files(path=many,pattern="*.csv")
if (length(files)>0)
{
    pfiles = paste0(many,"/",files)
#    print(pfiles)
    lst <- mclapply(pfiles,mc.cores=cores,function(fn){df <- read.csv(fn,header=T);df })

    onedf <- rbind(onedf,do.call(rbind,lst))
    write.csv(file=one,row.names=F,onedf)
    
    system(paste0("tar zcvf ",archive,"/Ash_stats_arch_",round(runif(1,1,1000),4),".tar.gz statsout")) 

    sapply(pfiles,unlink)
}


