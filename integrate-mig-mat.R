##
## calculate a "more correct" dispersal rate into rectangular landscape 
##

integratedMigMat <- function(xnum=4,ynum=4,xsz=100,ysz=100,sshp=1,ssc=10,mix=0.1,nmean=100,nvar=nmean)
{
  require(cubature)
 
  ynum1 = ynum+1
  xnum1 = xnum+1
  
  qmat = matrix(0,nrow = ynum1, ncol = xnum1)

 
  for (x in 1:xnum1)
  {
    if (x==1) lft=0.5*xsz else lft = (x-1)*xsz
    rgt = lft+xsz
    for (y in 1:ynum1)
    {
      if (y==1) bt = 0.5*ysz else bt = (y-1)*ysz
      tp = bt+ysz
 #     print(paste(lft,rgt))
 #     print(paste(bt,tp))
      qmat[y,x] <- adaptIntegrate(
        function(l)
          {
            (1-mix)*dweibull(l[1],shape=sshp,scale=ssc)+mix*dnorm(l[1],mean=nmean,sd=sqrt(nvar))+
            (1-mix)*dweibull(l[2],shape=sshp,scale=ssc)+mix*dnorm(l[2],mean=nmean,sd=sqrt(nvar))
          },
        c(lft,bt),c(rgt,tp))$integral
    }
  }
  
  qmat[1,1]=4*qmat[1,1]  #becaus we only integrated across 1/4 of central cell
  
  #make the single quadrant reflect across a 2d surface
  fmat <- matrix(0,nrow=1+ynum*2,ncol=1+xnum*2)
  fmat[(ynum)+1,] <-qmat[1,c(ynum1:1,2:ynum1)]
  fmat[,(xnum)+1] <-qmat[c(xnum1:1,2:xnum1),1]
  fmat[1:ynum,1:xnum] = qmat[ynum1:2,xnum1:2]
  fmat[1:ynum,(xnum+2):(2*xnum+1)] = qmat[ynum1:2,2:xnum1]
  fmat[(ynum+2):(ynum*2+1),] = fmat[ynum:1,]

  
  
  fmat
}
