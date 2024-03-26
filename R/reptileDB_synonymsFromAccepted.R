##' Get Reptile Database synonyms from accepted
##' 
##' Returns the synonyms that are associated with a given accepted taxon name.
##' 
##' @param sp genus and species (separated by underscore or space)
##' @param version which version of Reptile Database to use. Can be \code{'latest'} or a combination of month and year, such as \code{'3_2024'}. See \code{\link{reptileDB_versions}} for available versions.

##' 
##' @return Returns a list with a dataframe of synonyms for each queried accepted name.

##' @examples
##' 
##' reptileDB_synonymsFromAccepted(c('Phrynosoma_coronatum', 'Crotalus atrox'))
##' 
##' @export

reptileDB_synonymsFromAccepted <- function(sp, version = 'latest') {
	
	# is the requested version valid for ReptileDB?
	if (!isVersionValid(version)) {
		stop('Invalid version for Reptile DB.')
	}

	# check if Reptile DB data have already been made available via custom environment
	## if not, or if new version is requested, load it
	if (version == 'latest') {
		versionList <- reptileDB_versions()
		version <- versionList$version[1]
		# version <- paste(versionList$month[1], versionList$year[1], sep = '_')
	}
	if (length(names(.repDBvar)) == 0 | !identical(version, .repDBvar$version)) {
		.getRepDBfiles(version)
	}

	synonymTable <- .repDBvar$synonyms	
	
	components <- strsplit(sp, split = '_|\\s+')
	
	ret <- vector('list', length(sp))
	names(ret) <- sp
	
	for (i in 1:length(sp)) {
		
		ind <- which(synonymTable$acceptedGenus == components[[i]][1] & synonymTable$acceptedSpecies == components[[i]][2])
		
		if (length(ind) > 0) {
			xx <- synonymTable[ind, c('synonymGenus', 'synonymSpecies', 'synonymSubspecies', 'authority', 'year')]
			xx$synonym <- paste(xx$synonymGenus, xx$synonymSpecies, xx$synonymSubspecies, sep = '_')
			xx$synonym <- gsub('_NA$', '', xx$synonym)
			xx <- xx[, c('synonym', 'authority', 'year')]
			rownames(xx) <- NULL
		} else {
			xx <- as.data.frame(matrix(nrow = 0, ncol = 3))
			colnames(xx) <- c('synonym', 'authority', 'year')
		}
		
		ret[[i]] <- xx
	}

	return(ret)
}	
