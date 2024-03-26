#' @title List available versions of Reptile Database

#' @description As major updates to Reptile Database are released, they are made avaible via this R package.
#' Use this function to list the available versions. 

#' @return Returns a table of available versions, according to the month/year of download. If Reptile Database data have already been loaded, this table will also indicate which one is currently loaded.

#' @examples
#' reptileDB_versions()
#' @export

# function to provide list of available versions
reptileDB_versions <- function() {
		
	xx <- .getRepDBfiles(version = '')
	
	xx$date <- as.Date(paste(1, xx$month, xx$year, sep = '-'), format = '%e-%m-%Y')
	xx <- xx[order(xx$date, decreasing = TRUE),]
	xx <- xx[, c('month', 'year', 'version')]
	
	if (length(names(.repDBvar)) == 4) {
		xx$loaded <- xx$version == .repDBvar$version
	} else {
		xx$loaded <- FALSE
	}
	
	xx <- xx[!duplicated(xx),]
	rownames(xx) <- NULL
	return(xx)
}


isVersionValid <- function(version) {
	
	xx <- .getRepDBfiles(version = '')

	if (version == 'latest') {
		versionList <- reptileDB_versions()
		version <- versionList$version[1]
	}
	
	version %in% paste(xx$month, xx$year, sep = '_')
}
