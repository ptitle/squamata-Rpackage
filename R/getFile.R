## not exported

# @title Retrieve data files
# @description This function will download data files to a temporary directory, for use in this R package. If the file has already been downloaded and is found in the temporary directory, then it will not be re-downloaded.

# @param target which object is being retrieved
# @param seqClass if retrieving an alignment, which class to read in, DNAbin or DNAStringSet
# @return Returns the requested file. Could be phylo objects, data.frames or sequence data.



# this function will check if the file is already present in a certain directory, and if it is not, it will download it to that directory. 

.getFile <- function(target, seqClass = 'DNAbin') {
	
	targetOptions <- c('squamdat', 'ultrametric', 'molecular', 'unconstrained', 'pseudo-posterior', 'genomic', 'medianTimeTree',  'divTimeStats', 'timeTreesPost', 'alignment', 'geneIndex')
	
	target <- match.arg(target, choices = targetOptions)
	
	if (length(target) != 1) stop(paste0("target must be one of the following: ", paste0(targetOptions, collapse = ', ')))
	
	if (!requireNamespace("Biostrings", quietly = TRUE)) {
		seqClass <- 'DNAbin'
	}
	
	# define URL's and filenames
	
	if (target == 'squamdat') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/alldat.csv.xz"
		targetFile <- "alldat.csv.xz"
		
	} else if (target == 'ultrametric') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/squamates_Title_Science2024_ultrametric_constrained.tre.xz"
		targetFile <- "squamates_Title_Science2024_ultrametric_constrained.tre.xz"
		
	} else if (target == 'molecular') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/squamates_Title_Science2024_molecular_constrained.tre.xz"
		targetFile <- "squamates_Title_Science2024_molecular_constrained.tre.xz"
		
	} else if (target == 'unconstrained') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/squamates_Title_Science2024_unconstrained.tre.xz"
		targetFile <- "squamates_Title_Science2024_unconstrained.tre.xz"
		
	} else if (target == 'pseudo-posterior') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/pseudoposterior.100.trees.xz"
		targetFile <- "pseudoposterior.100.trees.xz"
		
	} else if (target == 'genomic') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/genomicTree.tre.xz"
		targetFile <- "genomicTree.tre.xz"
	
	} else if (target == 'medianTimeTree') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/medianTimeTree.tre.xz"
		targetFile <- "medianTimeTree.tre.xz"
	
	} else if (target == 'divTimeStats') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/divergenceTimeStats.csv.xz"
		targetFile <- "divergenceTimeStats.csv.xz"
	
	} else if (target == 'timeTreesPost') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/timeTreesPP.tre.xz"
		targetFile <- "timeTreesPP.tre.xz"
	
	} else if (target == 'alignment' & seqClass == 'DNAbin') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/concatenated.aln.xz"
		targetFile <- "concatenated.aln.xz"
		
	} else if (target == 'alignment' & seqClass == 'DNAStringSet') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/concatenated.aln.gz"
		targetFile <- "concatenated.aln.gz"
		
	} else if (target == 'geneIndex') {
		targetURL <- "https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/geneIndex.csv.xz"
		targetFile <- "geneIndex.csv.xz"

	}
		
	# Directory is a temporary session-specific location. It will be deleted at the end of the R session.
	targetFile <- file.path(tempdir(), targetFile)
	
	# Download the file if it has not yet been downloaded.
	if (!file.exists(targetFile)) {
		xx <- try(suppressWarnings(utils::download.file(url = targetURL, destfile = targetFile, quiet = TRUE)), silent = TRUE)
		if (inherits(xx, 'try-error')) stop('Download failed.')
	} 
	
	if (target %in% c('squamdat', 'divTimeStats', 'geneIndex')) {
		output <- utils::read.csv(xzfile(targetFile))
	} else if (target == 'alignment' & seqClass == 'DNAbin') {
		output <- ape::read.dna(xzfile(targetFile), format = 'fasta', as.matrix = TRUE)
	} else if (target == 'alignment' & seqClass == 'DNAStringSet') {
		output <- Biostrings::readDNAStringSet(targetFile)

	} else {
		output <- ape::read.tree(xzfile(targetFile))
	}

	return(output)
}


.getRepDBfiles <- function(version = 'latest') {
	
	# get the links table
	tableURL <- 'https://github.com/ptitle/squamata-Rpackage/raw/main/compressedDataFiles/reptileDB/repDBlinks.csv.xz'
	tableFile <- 'repDBlinks.csv.xz'
	
	# Directory is a temporary session-specific location. It will be deleted at the end of the R session.
	tableFile <- file.path(tempdir(), tableFile)
	
	# Download the table if it has not yet been downloaded.
	if (!file.exists(tableFile)) {
		xx <- try(suppressWarnings(utils::download.file(url = tableURL, destfile = tableFile, quiet = TRUE)), silent = TRUE)
		if (inherits(xx, 'try-error')) stop('Download failed.')
	}
	
	linksTable <- utils::read.csv(xzfile(tableFile))
	
	if (version != '') {
	
		# version must be either 'latest' or month_year
		if (version != 'latest' & !grepl('[0-9]{1,2}_[0-9]{4}', version)) {
			stop('version not recognized.')
		}
		
		linksTable$date <- as.Date(paste(1, linksTable$month, linksTable$year, sep = '-'), format = '%e-%m-%Y')
		linksTable <- linksTable[order(linksTable$date, decreasing = TRUE),]
		
		if (version == 'latest') {
			ind <- which(linksTable$date == linksTable$date[1])
			version <- paste(linksTable$month[ind[1]], linksTable$year[ind][1], sep = '_')
		} else {
			ind <- which(linksTable$month == strsplit(version, '_')[[1]][1] & linksTable$year == strsplit(version, '_')[[1]][2])
		}
		
		# download files if they are not already present
		urls <- linksTable[ind, 'link']
		files <- linksTable[ind, 'filename']
		files <- file.path(tempdir(), files)
		
		# Download the file if it has not yet been downloaded.
		for (i in 1:length(urls)) {
			if (!file.exists(files[i])) {
				xx <- try(suppressWarnings(utils::download.file(url = urls[i], destfile = files[i], quiet = TRUE)), silent = TRUE)
				if (inherits(xx, 'try-error')) stop('Download failed.')
			} 
		}
	
		# list.files(tempdir())
		
		# now read in those files
		acceptedNamesTable <- utils::read.csv(xzfile(files[grep('accepted', basename(files), ignore.case = TRUE)]))
		subspeciesTable <- utils::read.csv(xzfile(files[grep('subspecies', basename(files), ignore.case = TRUE)]))
		synonymTable <- utils::read.csv(xzfile(files[grep('synonym', basename(files), ignore.case = TRUE)]))
		
		# place in special environment to make them available
		## include version in environment too
		assign('acceptedNames', acceptedNamesTable, envir = .repDBvar)
		assign('subspecies', subspeciesTable, envir = .repDBvar)
		assign('synonyms', synonymTable, envir = .repDBvar)
		
		## how to add to the package environment
		assign('version', version, envir = .repDBvar)
	} else {
		return(linksTable)
	}
}






