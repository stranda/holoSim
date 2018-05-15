getpophist.Rversion <- function(l,maxtime=500, window = c(0,100),force=T)
{
  if (FALSE)
  {
    l <- recolonizeLandscape(dens.scale=0.05,
                             h=225,
                             refs=c(15,7,1),
                             sizeref=c(25,25,25),
                             mix=0.01,
                             longmean=3,
                             shortscale=0.35
                             )
  }
    pops <- landscape.popcoord(l)
    pops$arrive <- NA
    pops$source <- NA
    pops$init.abun <- NA
    pops$abun <- NA
    initpops <- unique(landscape.populations(l))
    pops$arrive[pops$pop%in%initpops] <- 0
    pops$source[pops$pop%in%initpops] <- 0
    
    if(window[1] == 0) {
      occpops = unique(landscape.populations(l))
      pops$init.abun[pops$pop%in%occpops] <- as.numeric(table(landscape.populations(l)))
    }

    cnt=1
    extinct=F
    while ((!extinct) & (cnt<=maxtime) & (length(which(is.na(pops$arrive) == TRUE)) > length(which(l$demography$epochs[[1]]$Carry == 0))))
    {
#        print(cnt)
        potmoth <- data.frame(unique(cbind(pop=landscape.populations(l),matid=l$individuals[,4])))

        l <- landscape.reproduce.cmp(l)
        l <-  landscape.survive.cmp(l)
	if (dim(l$individuals)[1]>0) l <- landscape.carry.cmp(l)
	if (dim(l$individuals)[1]>0) l <- landscape.advance.Rversion(l)
        
	if (dim(l$individuals)[1]>0) lpops = landscape.populations(l)
	if (dim(l$individuals)[1]>0)
	{	
        	if(window[1] == cnt) {
          	occpops = unique(lpops)
          	pops$init.abun[pops$pop%in%occpops] <- as.numeric(table(lpops))
	       	}
        
		if(window[2] == cnt) {
          	occpops = unique(lpops)
          	pops$abun[pops$pop%in%occpops] <- as.numeric(table(lpops))
        	}

        	ocpops <-  unique(lpops) #occupied pops
		newpops <- ocpops[which(ocpops %in% pops$pop[is.na(pops$arrive)])] #ones that were not occupied last gen
        if (length(newpops)>0) #figure out who colonized
        {
#            print(newpops)
            pops$arrive[pops$pop%in%newpops] <- l$intparam$currentgen
            for (p in newpops)
            {
                mothers <- potmoth[potmoth[,2] %in% unique(l$individuals[lpops==p,5]),]
                if (dim(mothers)[1]==1) #just one mother
                {
                    ret <- mothers$pop
                } else { #choose mother at random
                    ret <- sample(mothers$pop,1,F)
                }
                pops$source[pops$pop==p] <- ret
            }
        }
        cnt <- cnt+1
	} else { extinct=T;  }
    }
    pops
}


getpophist <- cmpfun(getpophist.Rversion)
