
.subsetDivStats <- function(phy) {

	timeTree <- .getFile(target = 'medianTimeTree')
	divStats <- .getFile(target = 'divTimeStats')
	
	nn <- (ape::Ntip(phy) + 1) : ape::Nnode(phy, internal.only = FALSE)
	
	patdist <- ape::cophenetic.phylo(phy)
	patdist[upper.tri(patdist)] <- 0
	
	newDivStats <- as.data.frame(matrix(nrow = length(nn), ncol = 6))
	colnames(newDivStats) <- colnames(divStats)
	for (i in 1:length(nn)) {
		nodeTips <- ape::extract.clade(phy, node = nn[i])$tip.label
		n1 <- ape::getMRCA(timeTree, nodeTips)
		if (abs(ape::branching.times(timeTree)[as.character(n1)] - ape::branching.times(phy)[as.character(nn[i])]) > 1e10) stop()
		newDivStats[i, 1] <- nn[i]
		newDivStats[i, 4:6] <- divStats[which(divStats[, 1] == n1), 4:6]
		maxdist <- which(patdist[nodeTips, nodeTips] == max(patdist[nodeTips, nodeTips]), arr.ind = TRUE)
		newDivStats[i, 2:3] <- c(rownames(patdist)[maxdist[1, 1]], colnames(patdist)[maxdist[1, 2]])
	}
	
	return(newDivStats)
}
