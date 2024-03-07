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
		targetURL <- "https://www.dropbox.com/scl/fi/24kf94a8w28zoz7u1azi5/alldat.csv.xz?rlkey=rfaiwuvlgtk4axwhhfoxg7o39&dl=1"
		targetFile <- "alldat.csv.xz"
		
	} else if (target == 'ultrametric') {
		targetURL <- "https://www.dropbox.com/scl/fi/xl1c1iliosuyb36exjf2l/best_ultrametric_fulltree_ddBD_revision.tre.xz?rlkey=bpy8da3v1uc5k370psrv5zs76&dl=1"	
		targetFile <- "best_ultrametric_fulltree_ddBD_revision.tre.xz"
		
	} else if (target == 'molecular') {
		targetURL <- "https://www.dropbox.com/scl/fi/0mldz0e2595uo5s7cxy5h/fulltree_default_con_1_raxmlOpt.raxml.final.tre.xz?rlkey=yugr83f1v69m68mqm5deqxb68&dl=1"
		targetFile <- "fulltree_default_con_1_raxmlOpt.raxml.final.tre.xz"
		
	} else if (target == 'unconstrained') {
		targetURL <- "https://www.dropbox.com/scl/fi/umqfx81oivfozijbhvmr6/fulltree_default_uncon_2_raxmlOpt.raxml.bestTree.xz?rlkey=srbxkg0rhl57yp9ug2rbp665q&dl=1"
		targetFile <- "fulltree_default_uncon_2_raxmlOpt.raxml.bestTree.xz"
		
	} else if (target == 'pseudo-posterior') {
		targetURL <- "https://www.dropbox.com/scl/fi/7eg59x22egq7gufdjjw8j/pseudoposterior.100.trees.xz?rlkey=8mt9yvz1mhm4llrcquis4hq8u&dl=1"
		targetFile <- "pseudoposterior.100.trees.xz"
		
	} else if (target == 'genomic') {
		targetURL <- "https://www.dropbox.com/scl/fi/vnxzljf47rjkbs21wltbg/genomicTree.tre.xz?rlkey=8t2xwitjvsbf5y2qhqekd26eu&dl=1"
		targetFile <- "genomicTree.tre.xz"
	
	} else if (target == 'medianTimeTree') {
		targetURL <- "https://www.dropbox.com/scl/fi/iah6vh0j85bl37gzqlp0w/medianTimeTree.tre.xz?rlkey=gmxnmnf31yopvz3kvcsoe7j5v&dl=1"
		targetFile <- "medianTimeTree.tre.xz"
	
	} else if (target == 'divTimeStats') {
		targetURL <- "https://www.dropbox.com/scl/fi/sjqznnjh2iyy8pnb7fumx/divergenceTimeStats.csv.xz?rlkey=1fym2z7wunq3a3kr5phyhgk8i&dl=1"
		targetFile <- "divergenceTimeStats.csv.xz"
	
	} else if (target == 'timeTreesPost') {
		targetURL <- "https://www.dropbox.com/scl/fi/1ojn1e6x7oqxlg0d6793y/timeTreesPP.tre.xz?rlkey=q717viarcg39ekywdamjwod7p&dl=1"
		targetFile <- "timeTreesPP.tre.xz"
	
	} else if (target == 'alignment' & seqClass == 'DNAbin') {
		targetURL <- "https://www.dropbox.com/scl/fi/nifxnuiz9k1jyb3er9o6n/concatenated.aln.xz?rlkey=px8ee3en8d7yw5iwaqgkk6by1&dl=1"
		targetFile <- "concatenated.aln.xz"
		
	} else if (target == 'alignment' & seqClass == 'DNAStringSet') {
		targetURL <- "https://www.dropbox.com/scl/fi/gulhua7m1bwzz1ob2juhc/concatenated.aln.gz?rlkey=53dylz5j6ia0m7qn9mpcb1mpf&dl=1"
		targetFile <- "concatenated.aln.gz"
		
	} else if (target == 'geneIndex') {
		targetURL <- "https://www.dropbox.com/scl/fi/ij2nahoa1rk2d6gy9p7w8/geneIndex.csv.xz?rlkey=qwsdr0cck44uzgxzumdjprafe&dl=1"
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
	tableURL <- 'https://www.dropbox.com/scl/fi/vtmqc3ml2wqp63s23jlv0/repDBlinks.csv.xz?rlkey=wedvjv6ynq3t4fjfta6rrlpdb&dl=1'
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






