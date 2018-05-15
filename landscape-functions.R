###
### quick way to pick multinomial data
###
quick_multinom.apply.r <- function(n,prob)
{
    v=rmultinom(n,1,prob)
    if (n>1)
    {
        apply(v==1,2,which)
    } else {
        which(v==1)
    }
}



quick_multinom.r <- function(n,prob)
{
    v=rmultinom(1,n,prob)
    inverse.rle(list(values=1:length(prob),lengths=v))
}

quick_multinom.apply <- cmpfun(quick_multinom.apply.r)
quick_multinom <- cmpfun(quick_multinom.r)



####
# fast reproduction using vectorized functions
# (trying to avoid rmetasim)
#
landscape.reproduce.Rversion <- function(rland)
{
  
###put the repro matrix in a nice form
    
    Rloc <- rland$demography$localdem[[1]]$LocalR
    R <-  rland$demography$epochs[[1]]$R #landscape reproduction
    md <- dim(Rloc)[1]
    for (i in seq(1,dim(R)[1],by=md))
        R[i:(i+(md-1)),i:(i+(md-1))] <- Rloc
    
    Rlong <- data.frame(cbind((which(R>0,arr.ind=T)-1),R[which(R>0,arr.ind=T)]))
    colnames(Rlong)[3] <- "rate"

    Rlong.tot=with(Rlong,aggregate(cbind(rate=rate),list(col=col),sum))

#    print(dim(Rlong.tot))
#print(str(Rlong.tot))
#    print(Rlong.tot$col)
#    print(fac2num(Rlong.tot$col))
    
#    Rlong.tot$col=fac2num(Rlong.tot$col) #convert factor to number as fast as possible
#    Rlong.tot$col=as.numeric(as.character(Rlong.tot$col))
  
    
###calculate the number of offspring per ind
  repind <- rland$individuals[rland$individuals[,1]%in%Rlong[,2],]
  if (dim(repind)[1]>0)
      {
          repindmn <- merge(Rlong.tot,data.frame(col=0:(dim(R)[1]-1)),all.y=T)[repind[,1]+1,2]
          repnoff <- repindmn
          for (i in unique(repindmn))
          {
              repnoff[repindmn==i] <- rpois(length(repindmn[repindmn==i]),i)
          }
          mother.indx <- inverse.rle(list(lengths=repnoff,values=1:dim(repind)[1]))
          newind <- matrix(0,nrow=length(mother.indx),ncol=6)
###       newind[,1] <- merge(Rlong[,2:1],data.frame(col=0:(dim(R)[1]-1)),all.y=T)[repind[mother.indx,1]+1,2]

          fromcols=repind[mother.indx,1]
          rlst = lapply(unique(fromcols),function(x) #divide the individuals up based on among pop probs
          {
###              v=rmultinom(sum(fromcols==x),1,R[,x+1]/sum(R[,x+1]))
###              if (is.matrix(v))
###              {
###                   apply(v==1,2,which)-1
###              } else {
###                  which(v==1)-1
###              }
              quick_multinom(sum(fromcols==x),prob=R[,x+1]/sum(R[,x+1]))-1
          })
          names(rlst) <- unique(fromcols)
          newrows = fromcols #just makingin a vector to overwrite
          for (nm in names(rlst))
          {
            newrows[fromcols==as.numeric(nm)]=rlst[[nm]]
          }
          
          newind[,1] <- newrows
          
          newind[,5] <- repind[mother.indx,4]
          newind[,6] <- repind[1,6] #not tracking fathers id
          newind[,4] <- rland$intparam$nextid:(rland$intparam$nextid+dim(newind)[1]-1) #set id #s
          rland$intparam$nextid <- (rland$intparam$nextid+dim(newind)[1])

#now decide where the offspring go:

          
          
          rland$individuals <- rbind(rland$individuals,
                                     cbind(newind,
                                           matrix(0,nrow=dim(newind)[1],
                                                  ncol=dim(rland$individuals)[2]-dim(newind)[2])
                                           )
                                     )
      }
  rland
}

#####################################################################################################


landscape.survive.Rversion <- function(rland)
  {
    Sloc <- rland$demography$localdem[[1]]$LocalS
    S <- rland$demography$epoch[[1]]$S
    md <- dim(Sloc)[1]
    for (i in seq(1,dim(S)[1],by=md))
        S[i:(i+(md-1)),i:(i+(md-1))] <- Sloc
    inds <- rland$individuals
    for (col in (unique(inds[,1])))
    {
        if (sum(S[,col+1])>=1)
        {
            p = S[,col+1]/sum(S[,col+1])
        } else {
            p=c(S[,col+1],1-sum(S[,col+1]))
            }
###        v=rmultinom(sum(inds[,1]==col),1,p)
###        if (is.matrix(v))
###        {
###            fates=apply(v==1,2,which)-1
###        } else {fates=which(v==1)-1}
###        inds[inds[,1]==col,1] <- fates
        inds[inds[,1]==col,1] <- quick_multinom(sum(inds[,1]==col),p)-1
    }
    inds <- inds[inds[,1]<(dim(S)[1]),]
    
    rland$individuals <- inds
    rland
  }

landscape.advance.Rversion <- function(rland)
{
    rland$intparam$currentgen=rland$intparam$currentgen+1
    rland
}

landscape.carry.Rversion <- function(rland)
{
    
    k = rland$demography$epoch[[1]]$Carry
    names(k) <- 1:length(k)
    lpops <- landscape.populations(rland)
    ps = table(lpops)
    lrg=c(ps)>k[names(ps)]
#    print(lrg)
#    print(ps)
#    inds=rland$individuals
        
    keeplst.big = unlist(sapply(names(lrg)[lrg],function(nm) 
    {
#        print(k[nm])
        unlist(sample(which(lpops==as.numeric(nm)),k[nm],replace=F))
    }))
    keeplst.small = unlist(sapply(names(lrg)[!lrg],function(nm) 
    {
        which(lpops==as.numeric(nm))
    }))

    keeplst <- c(keeplst.big,keeplst.small)
    rland$individuals = rland$individuals[keeplst,]
    rland
}


#"typelookup" <-
#function(type)
#  {
#    type <- type + 251
#    type
#  }

fac2num <- function(f)
{
    if (!is.numeric(f))
        {
            if((is.character(f)))
            {
                as.numeric(as.character(v))
            } else {
                l=levels(f)
                i=unclass(f)
                as.numeric(l[i])
            }
        } else {f}
                         
}

landscape.reproduce.cmp <- cmpfun(landscape.reproduce.Rversion)
landscape.survive.cmp <- cmpfun(landscape.survive.Rversion)
landscape.carry.cmp <- cmpfun(landscape.carry.Rversion)
