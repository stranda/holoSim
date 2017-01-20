#assuming a square landscape, what are the row,col coords for each pop
landscape.popcoord <- function(rland)
    {
        h <- rland$intparam$h
        nr=rland$intparam$ydim
        nc=rland$intparam$xdim
        if (sum(c(nr,nc)>0)!=2)  #assumes that xdim and ydim are correct
        {
          nr=sqrt(h)
          nc=sqrt(h)
        }
        df <- expand.grid(col=c(1:nc), row=c(1:nr))
        df$pop=1:h
        df[,c(3,1,2)]
    }

distancePDF <- function(x, ssh=1,ssc=1,lmn=100,lsd=100,mix=0)
    {
        (1-mix)*dweibull(x,shape=ssh,scale=ssc)+mix*dnorm(x,mean=lmn,sd=lsd)
    }

#creates a square landscape.  As a result, the square-root of 'h' must be a whole number
recolonizeLandscape <- function(
    h=100,
    s=2,
    k = rep(25,h),
    e =  rep(0,h),
    ncp = 1,
    cpmut=0.0,
    nssr = 1,
    ssrmut=0.0,
    niam = 1,
    iamut=0.0,
    refs=c(1),
    sizeref=c(25),
    shortshape=1,
    shortscale=0.3,
    longmean=2,
    mix = 0,       #ranges 0-1 0=all short, 1=all.long
    dens.scale = 0.05,
    xdim=0,
    ydim=0
    )

    {
  if (sum(c(xdim,ydim)>0)>0)  #theres a non-zero
    if (!(xdim*ydim)==h)
      {stop("the extent of the landscape should be the product of the x and y dimensions")}
  rland <- landscape.new.empty()
  rland <- landscape.new.intparam(rland, h = h, s = s, xdim=xdim, ydim=ydim)

  rland <- landscape.new.switchparam(rland, mp = 0)
  rland <- landscape.new.floatparam(rland)

  S <- matrix(c(0, 0, 1, 0), byrow = TRUE, nrow = 2)
  R <- matrix(c(0, 1.2, 0, 0), byrow = TRUE, nrow = 2)
  M <- matrix(c(0, 0, 0, 1), byrow = TRUE, nrow = 2)
  rland <- landscape.new.local.demo(rland, S, R, M)

  S <- matrix(0, nrow = s*h, ncol = s*h)
# R <- matrix(0, nrow = s*h, ncol = s*h)
  migmat <- landscape.mig.matrix(h=h,s=s,h.dim=c(ydim,xdim),mig.model="distance",distance.fun=distancePDF,
                                   ssc=shortscale,ssh=shortshape,
                                   lmn=longmean,lsd=longmean,mix=mix)
  R <- migmat$R*dens.scale
  M <- matrix(0, nrow = s*h, ncol = s*h)
  for (row in 1:dim(migmat$R.int)[1])
      for (col in 1:dim(migmat$R.int)[2])
          {
              M[2+(row-1)*2,2+(col-1)*2] <- migmat$R.int[row,col]
          }
  M <- M*dens.scale

  rland <- landscape.new.epoch(rland, S = S, R = R, M = M,
                               carry = k,
                               extinct = e)
  for (l in 1:ncp)
      rland <- landscape.new.locus(rland, type = 0, ploidy = 1,
                                   mutationrate = cpmut, transmission = 1, numalleles = 5)
  for (l in 1:nssr)
      rland <- landscape.new.locus(rland, type = 1, ploidy = 2,
                                   mutationrate = ssrmut, transmission = 0, numalleles = 5)
  for (l in 1:niam)
      rland <- landscape.new.locus(rland, type = 0, ploidy = 2,
                                   mutationrate = iamut, transmission = 0, numalleles = 5)
  popsizes <- rep(c(0,0),h)

  for (i in 1:length(refs))
      {
          p=refs[i]
          ps=sizeref[i]
          popsizes[(1+((p-1)*2)):(2+((p-1)*2))] <- c(ps,0)
      }
  rland <- landscape.new.individuals(rland, popsizes)
  rland
}
