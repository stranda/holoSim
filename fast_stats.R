getStats = function(out) {

	freqs = alleleFreqs(out, TRUE)
	popid = strataNames(out)
	split_out = strataSplit(out)

	freqmat = matrix(data = NA, nrow = length(freqs), ncol = length(popid))
	for(p in 1:length(freqmat[1,])) {
		freqmat[,p] = sapply(freqs, FUN=function(x){a=x[,,popid[p]][3]})
	}
	colnames(freqmat) = popid

	variable_loci = matrix(data = 0, nrow = length(freqs), ncol = length(popid))
	for(p in 1:length(variable_loci[1,])) {
		variable_loci[freqmat[,p] > 0 & freqmat[,p] < 1 & !is.na(freqmat[,p]),p] = 1
	}
	colnames(variable_loci) = popid

	nind_mat = matrix(data = NA, nrow = length(freqs), ncol = length(popid))
	for(p in 1:length(nind_mat[1,])) {
		nind_mat[,p] = sapply(freqs, FUN=function(x){a=x[,,popid[p]][1]+x[,,popid[p]][2]})
	}
	colnames(nind_mat) = popid

	#Local SNP number
	localSNP = colSums(variable_loci)
	names(localSNP) = paste0("S.", popid)

	Smean = mean(localSNP)
	names(Smean) = "mean_S"

	Ssd = sd(localSNP)
	names(Ssd) = "sd_S"

	#Private SNP number
	privateSNP = c()
	for(p in 1:length(variable_loci[1,])) {
		privateSNP[p] = length(which(rowSums(variable_loci) == 1 & variable_loci[,p] == 1))
	}
	names(privateSNP) = paste0("pS.", popid)

	prSmean = mean(privateSNP)
	names(prSmean) = "mean_prS"

	prSsd = sd(privateSNP)
	names(prSsd) = "sd_prS"

	total_priv = sum(privateSNP)
	names(total_priv) = "tot_priv"

	#Pairwise private SNPs
	combos = combn(popid, 2)
	pwPriv = c()
	for(pair in 1:length(combos[1,])) {
		subset = which(nind_mat[,combos[1,pair]] > 0 & nind_mat[,combos[2,pair]] > 0)
		tmpvar = variable_loci[subset,c(combos[1,pair], combos[2,pair])]
		tmp = c(length(which(tmpvar[,1] == 1 & tmpvar[,2] == 0)), length(which(tmpvar[,2] == 1 & tmpvar[,1] == 0)))
		names(tmp) = c(paste0("pwPriv.",combos[1,pair],".",combos[2,pair]), paste0("pwPriv.",combos[2,pair],".",combos[1,pair]))
		pwPriv = c(pwPriv, tmp)
		rm(tmpvar, tmp, subset)
	}

	#He
	Hetmp = sapply(split_out, FUN=function(x){exptdHet(x)})
	He = colMeans(Hetmp, na.rm = TRUE)
	names(He) = paste0("He.",  colnames(Hetmp))

	HeALL = mean(exptdHet(out))
	names(HeALL) = "HeALL"

	#Simple Fst calculation...
	pwFst = c()
	for(pair in 1:length(combos[1,])) {
		subset = which(nind_mat[,combos[1,pair]] > 0 & nind_mat[,combos[2,pair]] > 0)
		Hs_tmp = (nind_mat[subset,combos[1,pair]]*Hetmp[subset,combos[1,pair]] + nind_mat[subset,combos[2,pair]]*Hetmp[subset,combos[2,pair]])/(nind_mat[subset,combos[1,pair]]+nind_mat[subset,combos[2,pair]])
		combo_freq = (nind_mat[subset,combos[1,pair]]*freqmat[subset,combos[1,pair]] + nind_mat[subset,combos[2,pair]]*freqmat[subset,combos[2,pair]])/(nind_mat[subset,combos[1,pair]]+nind_mat[subset,combos[2,pair]])
		Ht_tmp = 1-combo_freq^2-(1-combo_freq)^2
		Fst_tmp = (Ht_tmp-Hs_tmp)/Ht_tmp
		pwFst[pair] = mean(Fst_tmp, na.rm = TRUE)
		names(pwFst)[pair] = paste0(combos[1,pair],".",combos[2,pair])
		rm(subset, Hs_tmp, combo_freq, Ht_tmp, Fst_tmp)
	}

	stats = c(localSNP, Smean, Ssd, privateSNP, prSmean, prSsd, total_priv, pwPriv, He, HeALL, pwFst)
	stats
}

