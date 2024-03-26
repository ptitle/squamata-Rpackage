#' @title Resolve taxonomy to match the Reptile Database

#' @description This function will attempt to assign taxon names to a set of Reptile Database accepted names, using both strict and fuzzy matching, and via synonymy.

#' @param taxon a single taxon name
#' @param fuzzyThresh the acceptable number of differences *per* component (genus, sp, subsp) for fuzzy matching. 
#' 	If zero, no fuzzy matching. See \code{\link[utils]{adist}}. 
#' @param alternateEndings allow for testing of alternate latin suffixes when matching to accepted names
#' @param retainSubspecies If subspecies provided, should subspecies be retained in result. 
#' 	This only pertains to cases where the taxon matches to a known subspecies. As synonyms are mapped 
#' 	to Genus_species accepted names (with no subspecies), if a taxon matches to a synonym with subspecies, 
#' 	it will still map back to a Genus_species.
#' @param yearCutoff When matching to synonymys, what is the earliest publication year to consider. 
#' 	Can be set to NULL to be ignored
#' @param version which version of Reptile Database to use. Can be \code{'latest'} or a combination of month and year, such as \code{'3_2024'}. See \code{\link{reptileDB_versions}} for available versions.

#' @return A dataframe with matches to either accepted taxon names or synonyms. 
#' 	Authority and year pertain to the accepted name if directly matched, otherwise pertain to the synonym. 
#' 	If a synonym was listed multiple times with different years under the same accepted name, then the oldest authority/year is listed.

#' @details
#' This function makes the following attempts (assuming all options enabled):
#' 
#'	- strict matching to accepted names
#'	- strict matching to accepted names, but with alternate latin suffixes
#'	- fuzzy matching to accepted names
#'	- strict matching to synonyms
#'	- fuzzy matching to synonyms
#' 
#' It is possible to match to both an accepted name and to the synonym of another accepted name.

#' @examples
#'
#' # exact matching to accepted name
#' reptileDB_match('Lachesis stenophrys')
#'
#' # fuzzy matching to accepted name
#' reptileDB_match('Anolis saxitilis')
#'
#' # exact matching to synonyms
#' reptileDB_match('Japalura planidorsata')
#'
#' # fuzzy matching to synonyms
#' reptileDB_match('Amphisbaena mensa')
#'
#' # matching via alternate latin suffixes
#' reptileDB_match('Magliophis exiguum')
#'
#' # inclusion of subspecies
#' reptileDB_match('Phrynosoma coronatum blainvillii')
#'
#' # multiple hits
#' reptileDB_match('Chironius flavolineatus')
#'
#' # querying multiple taxa
#' lapply(c('Coluber_constrictor', 'Anolis saxitilis', 'Xenochrophis flavipunctatus'), reptileDB_match) 
#'
#' # query a taxon that did not yet exist (described in 2022)
#' reptileDB_match('Abronia zongolica', version = 'latest')
#' reptileDB_match('Abronia zongolica', version = '1_2020')

# # c('Coluber_constrictor', 'Anolis saxitilis', 'Xenochrophis flavipunctatus')

# # exact mapping via accepted
# c('Liolaemus capillitas', 'Lachesis stenophrys', 'Phymaturus laurenti', 'Langaha madagascariensis', 'Sceloporus virgatus')

# # fuzzy mapping via accepted 
# c('Xenochrophis bellula', 'Aspidoscelis neomexicana', 'Magliophis exiguum') 

# # mapping via synonymy
# c('Japalura planidorsata', 'Mabuya ficta', 'Tricheilostoma dimidiatum', 'Amphisbaena mensae', 'Lycodon osmanhilli')

# # multiple matches
# c('Chironius flavolineatus', 'Vipera ursinii', 'Storeria dekayi', 'Anolis singularis', 'Sceloporus graciosus')

# # no match
# c('Acanthodactylus cantoris', 'Proablepharus tenuis', 'Calyptotis thortonensis', 'Macroprotodon cucullatus', 'Calyptotis thortonensis')

# # with subspecies
# c('Eumeces schneideri barani', 'Lacerta trilineata hansschweizeri', 'Nerodia sipedon pleuralis', 'Chilabothrus strigilatus mccraniei', 'Tarentola boettgeri boettgeri')
#'
#' @export

reptileDB_match <- function(taxon, fuzzyThresh = 1, alternateEndings = TRUE, retainSubspecies = FALSE, yearCutoff = 1950, version = 'latest') {
	
	# fuzzyThresh = 1; alternateEndings = TRUE; retainSubspecies = FALSE; yearCutoff = 1950; version = 'latest'
	
	# notes
	## if matches to synonym, but that synonym is identical to its accepted name, then ignored
	## if multiple synonyms matched, but they map to the same accepted species, then only keep one - the oldest one?
	
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
	# make more convenient by separating out the tables
	acceptedNamesTable <- .repDBvar$acceptedNames[, c('genus', 'species', 'authority', 'year')]
	subspeciesTable <- .repDBvar$subspecies
	synonymTable <- .repDBvar$synonyms
	
	
	# data check
	if (anyNA(synonymTable[, c('acceptedGenus', 'acceptedSpecies', 'synonymGenus', 'synonymSpecies')])) {
		stop("Synonym table should not contain NA's.")
	}
	if (anyNA(acceptedNamesTable[, c('genus', 'species')])) {
		stop("Accepted names table should not contain NA's.")
	}
	if (anyNA(subspeciesTable[, c('genus', 'species', 'subspecies')])) {
		stop("Accepted names table should not contain NA's.")
	}
	
	
	if (is.null(yearCutoff)) {
		yearCutoff <- min(synonymTable$year, na.rm = TRUE)
	}
		
	# split into components
	components <- strsplit(taxon, split = '_|\\s+')[[1]]
	genus <- components[1]
	species <- components[2]
	subspecies <- ifelse(length(components) == 3, components[3], NA)
		
	
	## --------------------------------------
	## EXACT MATCHING TO ACCEPTED TAXON NAMES
	## --------------------------------------

	# Is there an exact match with a reptile DB accepted name?
	if (!is.na(subspecies)) {
		acceptedMatch <- which(subspeciesTable$genus == genus & subspeciesTable$species == species & subspeciesTable$subspecies == subspecies)
		acceptedMatch <- subspeciesTable[acceptedMatch,]
				
	} else {	
		acceptedMatch <- which(acceptedNamesTable$genus == genus & acceptedNamesTable$species == species)
		acceptedMatch <- acceptedNamesTable[acceptedMatch,]
				
	}
	
	if (nrow(acceptedMatch) > 0) {
		acceptedMatch$approach <- 'exact match to accepted name'
	}


	## ---------------------------------------------------------------
	## ALTERNATE LATIN ENDINGS, EXACT MATCHING TO ACCEPTED TAXON NAMES
	## ---------------------------------------------------------------
	if (nrow(acceptedMatch) == 0 & alternateEndings) {
		
		## Consider errors in latin endings
		speciesOptions <- c(species, 
							gsub('a$', 'um', species),
							gsub('a$', 'is', species),
							gsub('a$', 'us', species),
							gsub('um$', 'us', species),
							gsub('us$', 'um', species),
							gsub('um$|is$|us$', 'a', species),
							gsub('[i]{2}$', 'i', species),
							gsub('[i]{1}$', 'ii', species))
		
		speciesOptions <- setdiff(unique(speciesOptions), species)
									
		acceptedMatchVec <- vector('list', length(speciesOptions))
		if (length(speciesOptions) > 0) {
			for (i in 1:length(speciesOptions)) {
			
				# Is there an exact match with a reptile DB accepted name?
				if (!is.na(subspecies)) {
					acceptedMatch <- which(subspeciesTable$genus == genus & subspeciesTable$species == speciesOptions[i] & subspeciesTable$subspecies == subspecies)
					acceptedMatchVec[[i]] <- subspeciesTable[acceptedMatch,]
										
				} else {	
					acceptedMatch <- which(acceptedNamesTable$genus == genus & acceptedNamesTable$species == speciesOptions[i])
					acceptedMatchVec[[i]] <- acceptedNamesTable[acceptedMatch,]
								
				}
			}
		}
		
		acceptedMatchVec <- do.call(rbind, acceptedMatchVec)
		
		# if any alternate endings helped, use that match.
		if (length(acceptedMatchVec) > 0) {
			acceptedMatch <- acceptedMatchVec
		} else {
			acceptedMatch <- acceptedNamesTable[integer(0),]
		}
		
		if (nrow(acceptedMatch) > 0) {
			acceptedMatch$approach <- 'exact match to accepted name using alternate latin endings'
		}
	}
	
	
	## ---------------------------------------------------------------
	## FUZZY MATCHING TO ACCEPTED TAXON NAMES
	## ---------------------------------------------------------------
	if (nrow(acceptedMatch) == 0 & fuzzyThresh > 0) {
			
		# What is the minimum distance to an accepted name and is it below threshold?
		if (!is.na(subspecies)) {
			acceptedDist_genus <- utils::adist(genus, subspeciesTable$genus)[1,]
			acceptedDist_species <- utils::adist(species, subspeciesTable$species)[1,]
			acceptedDist_subspecies <- utils::adist(subspecies, subspeciesTable$subspecies)[1,]
			
			if (any(acceptedDist_genus <= fuzzyThresh & acceptedDist_species <= fuzzyThresh & acceptedDist_subspecies <= fuzzyThresh)) {
				acceptedMatch <- which(acceptedDist_genus <= fuzzyThresh & acceptedDist_species <= fuzzyThresh & acceptedDist_subspecies <= fuzzyThresh)
			} else {
				acceptedMatch <- integer(0)
			}
			acceptedMatch <- subspeciesTable[acceptedMatch,]
						
		} else {	
			acceptedDist_genus <- utils::adist(genus, acceptedNamesTable$genus)[1,]
			acceptedDist_species <- utils::adist(species, acceptedNamesTable$species)[1,]
			
			if (any(acceptedDist_genus <= fuzzyThresh & acceptedDist_species <= fuzzyThresh)) {
				acceptedMatch <- which(acceptedDist_genus == min(acceptedDist_genus) & acceptedDist_species == min(acceptedDist_species))
			} else {
				acceptedMatch <- integer(0)
			}
			acceptedMatch <- acceptedNamesTable[acceptedMatch,]
					
		}
		
		if (nrow(acceptedMatch) > 0) {
			acceptedMatch$approach <- 'fuzzy match to accepted name'
		}	
	}
	
	## ---------------------------------------------------------------
	## EXACT MATCHING TO SYNONYMS
	## ---------------------------------------------------------------
	
	# Regardless of whether the taxon matched an accepted name, is it also present in synonyms for a different taxon?
	if (!is.na(subspecies)) {
		synonymMatch <- which(synonymTable$synonymGenus == genus & synonymTable$synonymSpecies == species & synonymTable$synonymSubspecies == subspecies & (synonymTable$year >= yearCutoff | is.na(synonymTable$year)))
		synonymMatch <- synonymTable[synonymMatch, ]
		
	} else {
		synonymMatch <- which(synonymTable$synonymGenus == genus & synonymTable$synonymSpecies == species & (synonymTable$year >= yearCutoff | is.na(synonymTable$year)))
		synonymMatch <- synonymTable[synonymMatch, ]
	}
	
	# exclude synonyms that exactly match their associated accepted names
	ident <- synonymMatch$acceptedGenus == synonymMatch$synonymGenus & synonymMatch$acceptedSpecies == synonymMatch$synonymSpecies
	synonymMatch <- synonymMatch[!ident,]
	
	synonymMatch <- synonymMatch[, c('acceptedGenus', 'acceptedSpecies', 'authority', 'year')]
	colnames(synonymMatch)[1:2] <- c('genus', 'species')
	
	if (nrow(synonymMatch) > 0) {
		synonymMatch$approach <- 'exact match to synonym'
	}

	## ---------------------------------------------------------------
	## FUZZY MATCHING TO SYNONYMS
	## ---------------------------------------------------------------
	
	if (nrow(synonymMatch) == 0 & fuzzyThresh > 0) {

		# What is the minimum distance to a synonym and is it below threshold?
		if (!is.na(subspecies)) {
			synonymDist_genus <- utils::adist(genus, synonymTable$synonymGenus)[1,]
			synonymDist_species <- utils::adist(species, synonymTable$synonymSpecies)[1,]
			synonymDist_subspecies <- utils::adist(subspecies, synonymTable$synonymSubspecies)[1,]
			
			if (any(synonymDist_genus <= fuzzyThresh & synonymDist_species <= fuzzyThresh & synonymDist_subspecies <= fuzzyThresh)) {
				synonymMatch <- which(synonymDist_genus <= fuzzyThresh & synonymDist_species <= fuzzyThresh & synonymDist_subspecies <= fuzzyThresh & (synonymTable$year >= yearCutoff | is.na(synonymTable$year)))
			} else {
				synonymMatch <- integer(0)
			}
						
			synonymMatch <- synonymTable[synonymMatch,]

			# it's possible that the matches have a range of fuzzy similarity.
			## If true, keep only the one(s) with the lowest score. 
			if (length(synonymMatch) > 1) {
				matchDist <- utils::adist(genus, synonymMatch$synonymGenus)[1,] + utils::adist(species, synonymMatch$synonymSpecies)[1,] + utils::adist(subspecies, synonymMatch$synonymSubspecies)[1,]
				synonymMatch <- synonymMatch[matchDist == min(matchDist), ]
			}
				
		} else {	
			synonymDist_genus <- utils::adist(genus, synonymTable$synonymGenus)[1,]
			synonymDist_species <- utils::adist(species, synonymTable$synonymSpecies)[1,]
			
			if (any(synonymDist_genus <= fuzzyThresh & synonymDist_species <= fuzzyThresh)) {
				synonymMatch <- which(synonymDist_genus <= fuzzyThresh & synonymDist_species <= fuzzyThresh & (synonymTable$year >= yearCutoff | is.na(synonymTable$year)))
			} else {
				synonymMatch <- integer(0)
			}
			synonymMatch <- synonymTable[synonymMatch,]
			
			# it's possible that the matches have a range of fuzzy similarity.
			## If true, keep only the one(s) with the lowest score. 
			if (nrow(synonymMatch) > 1) {
				matchDist <- utils::adist(genus, synonymMatch$synonymGenus)[1,] + utils::adist(species, synonymMatch$synonymSpecies)[1,]
				synonymMatch <- synonymMatch[matchDist == min(matchDist), ]
			}
		}
		
		# exclude synonyms that exactly match their associated accepted names
		ident <- synonymMatch$acceptedGenus == synonymMatch$synonymGenus & synonymMatch$acceptedSpecies == synonymMatch$synonymSpecies
		synonymMatch <- synonymMatch[!ident,]					

		synonymMatch <- synonymMatch[, c('acceptedGenus', 'acceptedSpecies', 'authority', 'year')]
		colnames(synonymMatch)[1:2] <- c('genus', 'species')
				
		if (nrow(synonymMatch) > 0) {
			synonymMatch$approach <- 'fuzzy match to synonym'
		}
	}

	
	# --------------------------------------
	# Organize results
	
	# add in handling of subspecies
	if (!'subspecies' %in% colnames(acceptedMatch) & nrow(acceptedMatch) > 0) {
		acceptedMatch$subspecies <- NA
	}
		
	# if multiple identical synonym mappings, just keep the oldest one
	if (nrow(synonymMatch) > 0) {
		synonymMatch$subspecies <- NA
		synonymMatch <- do.call(rbind, lapply(split(synonymMatch, paste(synonymMatch$genus, synonymMatch$species, synonymMatch$subspeciess, sep = '_')), \(x) x[which(x$year == min(x$year, na.rm = TRUE)),]))
	}
	
	# did taxon get matched to accepted name and then again to that same name via synonymy?
	## If so, we can remove that synonym.
	ind <- which(synonymMatch$genus %in% acceptedMatch$genus & synonymMatch$species %in% acceptedMatch$species)
	if (length(ind) > 0) {
		# stop()
		synonymMatch <- synonymMatch[setdiff(1:nrow(synonymMatch), ind), ]
	}
		
	res <- rbind.data.frame(acceptedMatch, synonymMatch[order(synonymMatch$year, decreasing = TRUE),])
	res$acceptedName <- paste(res$genus, res$species, res$subspecies, sep = '_')
	res$acceptedName <- gsub('_NA$', '', res$acceptedName)
	if (!retainSubspecies) {
		res$acceptedName <- gsub('(^[A-z]+)_([a-z]+)_([a-z]+$)', '\\1_\\2', res$acceptedName)
	}
	
	res$matchType <- ifelse(grepl('accepted', res$approach), 'accepted', 'synonym')

	if (sum(res$matchType == 'accepted') > 1) stop('more than one accepted match.')
			
	if (nrow(res) > 0) {
		res$query <- taxon	
		res <- res[, c('query', 'acceptedName', 'authority', 'year', 'matchType', 'approach')]
		res <- res[!duplicated(res),]
		rownames(res) <- NULL
	
	} else {
		
		res <- as.data.frame(matrix(nrow = 0, ncol = 6))
		colnames(res) <- c('query', 'acceptedName', 'authority', 'year', 'matchType', 'approach')		
	}
	
	return(res)
}

	


	
