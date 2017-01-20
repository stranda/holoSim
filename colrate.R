#Maximum Euclidean Distance from the source population 
#	More challenging for more than 1 refuge?
#How to extend this to 2 or more sources?
#	take max of the average Euclidean distance from sources?
#	take max of min Euclidean distance from sources?
#Right now this only works from the start of the simulation -- would need to modify to look at other windows (r1,c1 would change)
getC = function(p=NULL, window = c(0,100)) {
	r1 = p$row[which(p$arrive == 0)]
	c1 = p$col[which(p$arrive == 0)]
	tmp = p[p$arrive %in% seq(window[1],window[2]),]
	dis = c()
	for(i in 1:length(tmp[,1])) {
		dis[i] = sqrt((tmp$row[i]-r1)^2+(tmp$col[i]-c1)^2)
	}
	C = max(dis)
	C
}

#Northward movement only
getCn = function(p=NULL, window = c(0,100)) {
	r1 = p$row[which(p$arrive == 0)]
	tmp = p[p$arrive %in% seq(window[1],window[2]),]
	if(length(which(tmp$row > r1)) > 0) {
		tmp = tmp[tmp$row > r1,]
		dis = c()
		for(i in 1:length(tmp[,1])) {
			dis[i] = sqrt((tmp$row[i]-r1)^2)
		}
		C = max(dis)
	} else {
		C = 0
	}
	C
}

#Southward movement only (this should always be 0 if TX refuge)
getCs = function(p=NULL, window = c(0,100)) {
	r1 = p$row[which(p$arrive == 0)]
	tmp = p[p$arrive %in% seq(window[1],window[2]),]
	if(length(which(tmp$row < r1)) > 0) {
		tmp = tmp[tmp$row < r1,]
		dis = c()
		for(i in 1:length(tmp[,1])) {
			dis[i] = sqrt((tmp$row[i]-r1)^2)
		}
		C = max(dis)
	} else {
		C = 0
	}
	C
}

#Focus on growth in area occupied instead of distance from source? - overly simplistic?
#	no distinguishing between situations where all occupied cells are very close to source vs. much more spread out
#This will be more about the number of cells occupied relative to the starting point
getC2 = function(p=NULL, window = c(0,100)) {
	end = p[p$arrive <= window[2] & p$arrive >= window[1] & is.na(p$arrive) == FALSE,]
	start = p[p$arrive <= window[1] & is.na(p$arrive) == FALSE,]
	C = length(end[,1]) - length(start[,1])
	C
}

#No way to make this N / S specific


#Centroid of occupied grid cells (assuming equal density) 
#	probably does a better job of capturing long distance directional movement than C2?
#	also may be inherently better for more than one time slice?
#How much has this point changed?
getC3 = function(p=NULL, window = c(0,100)) {
	end = p[p$arrive <= window[2] & is.na(p$arrive) == FALSE,]
	start = p[p$arrive <= window[1] & is.na(p$arrive) == FALSE,]
	r1 = median(start$row)
	r2 = median(end$row)
	c1 = median(start$col)
	c2 = median(end$col)

	C = sqrt((r1-r2)^2+(c1-c2)^2)
	C
}

#Northward
getC3n = function(p=NULL, window = c(0,100)) {
	start = p[p$arrive <= window[1] & is.na(p$arrive) == FALSE,]
	r1 = median(start$row)

	end = p[p$arrive <= window[2] & is.na(p$arrive) == FALSE,]
	if(length(which(end$row > r1)) > 0) {
		end = end[end$row > r1,]
		r2 = median(end$row)
		C = sqrt((r1-r2)^2)
	} else {
		C = 0
	}
	C
}

#Southward
getC3s = function(p=NULL, window = c(0,100)) {
	start = p[p$arrive <= window[1] & is.na(p$arrive) == FALSE,]
	r1 = median(start$row)

	end = p[p$arrive <= window[2] & is.na(p$arrive) == FALSE,]
	if(length(which(end$row < r1)) > 0) {
		end = end[end$row < r1,]
		r2 = median(end$row)
		C = sqrt((r1-r2)^2)
	} else {
		C = 0
	}
	C
}

#Weighting range by population density
#Centroid of individual locations?
#This requires a change in getpophist.R -- pops object needs to record abundance at time t
#Change will introduce 2 new columns into pops object -- abundance - 'abun' & 'init.abun'
getC4 = function(p=NULL, window = c(0,100)) {
	end = p[p$arrive <= window[2] & is.na(p$arrive) == FALSE,]
	start = p[p$arrive <= window[1] & is.na(p$arrive) == FALSE,]

	start.weights = start$init.abun/sum(start$init.abun)
	end$abun[is.na(end$abun) == TRUE]=0		#I don't know that I like this, but at least it gives a value.  What about extinctions?!? 
	end.weights = end$abun/sum(end$abun)

	start.row = sum(start$row*start.weights)	#These average x,y positions on the grid weighted by abundance
	start.col = sum(start$col*start.weights)
	end.row = sum(end$row*end.weights)
	end.col = sum(end$col*end.weights)

	C = sqrt((start.row-end.row)^2 + (start.col-end.col)^2)
	C

}

#Northward
getC4n = function(p=NULL, window = c(0,100)) {
	start = p[p$arrive <= window[1] & is.na(p$arrive) == FALSE,]
	start.weights = start$init.abun/sum(start$init.abun)
	start.row = sum(start$row*start.weights)	#These average x,y positions on the grid weighted by abundance

	end = p[p$arrive <= window[2] & is.na(p$arrive) == FALSE,]
	if(length(which(end$row > start.row)) > 0) {
		end = end[end$row > start.row,]
		end$abun[is.na(end$abun) == TRUE]=0		#I don't know that I like this, but at least it gives a value.  What about extinctions?!? 
		end.weights = end$abun/sum(end$abun)
		end.row = sum(end$row*end.weights)
		C = sqrt((start.row-end.row)^2)
	} else {
		C = 0
	}
	C
}

#Southward
getC4s = function(p=NULL, window = c(0,100)) {
	start = p[p$arrive <= window[1] & is.na(p$arrive) == FALSE,]
	start.weights = start$init.abun/sum(start$init.abun)
	start.row = sum(start$row*start.weights)	#These average x,y positions on the grid weighted by abundance

	end = p[p$arrive <= window[2] & is.na(p$arrive) == FALSE,]
	if(length(which(end$row < start.row)) > 0) {
		end = end[end$row < start.row,]
		end$abun[is.na(end$abun) == TRUE]=0		#I don't know that I like this, but at least it gives a value.  What about extinctions?!? 
		end.weights = end$abun/sum(end$abun)
		end.row = sum(end$row*end.weights)	
		C = sqrt((start.row-end.row)^2)
	} else {
		C = 0
	}
	C

}