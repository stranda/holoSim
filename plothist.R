plothist <- function(pops)
{
    layout(matrix(c(1,1,1,1,2,3),nrow=3,ncol=2,byrow=T))
    rows <- min(pops$row):max(pops$row)
    cols <- min(pops$col):max(pops$col)
    pops$coldist=NA
    plot(1,type="n",ylab="",xlab="",ylim=range(rows),xlim=range(cols),axes=F)
    points(row~col,pops,pch=16,cex=1.6)
    for (i in 1:dim(pops)[1])
    {
        if (pops$arrive[i]!=0)
        {
            
            x1=pops[i,"col"]
            y1=pops[i,"row"]
            x0= pops[pops$pop==pops[i,"source"],"col"]
            y0= pops[pops$pop==pops[i,"source"],"row"]
#            print(c(x0,y0,x1,y1))
            arrows(x0,y0,x1,y1,
                   col=heat.colors(max(pops$arrive))[pops$arrive[i]],
                   lwd=2,
                   length=0.1
               )
            pops$coldist[i] <- sqrt((y0-y1)^2 + (x0-x1)^2)
        }
    }

    hist(pops$coldist,xlab="Colonization distance",main="")
    hist(pops$arrive,xlab="Colonization time",main="")
    layout(matrix(1))
}
