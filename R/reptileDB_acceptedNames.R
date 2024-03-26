##' Get Reptile Database accepted names
##' 
##' Returns the full list of accepted names
##'
##' @param version which version of Reptile Database to use. Can be \code{'latest'} or a combination of month and year, such as \code{'3_2024'}. See \code{\link{reptileDB_versions}} for available versions.

##' 
##' 
##' @return
##' Returns the list of accepted species names in the
##' database along with higher taxonomy.
##' 
##' @examples
##' 
##' head(reptileDB_acceptedNames())
##' 
##' @export

reptileDB_acceptedNames <- function(version = 'latest') {
	
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

	.repDBvar$acceptedNames

}