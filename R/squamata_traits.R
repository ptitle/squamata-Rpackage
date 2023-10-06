#' @title Retrieve eco-morphological data

#' @description This function will return species-level ecological, morphological, 
#' 	environmental, geographical data compiled for the study of squamate macroevolution.
#' 	Trait tip rates and net innovation index are also available. 

#' @param species a vector of species names that the tree will be subset to, Default: NULL
#' @param genus a vector of genus names that the tree will be subset to, Default: NULL
#' @param family a vector of family names that the tree will be subset to, Default: NULL
#' @param subclade a vector of major squamate subclades that the tree will be subset to, Default: NULL
#' @param type which traits to return; defaults to \code{'all'}

#' @return a data.frame with taxonomic information and trait data
#' @details 
#'	The list of major subclades that can be used for subsetting are listed under \code{\link{squamata}}.
#'	
#'	Argument \code{type} can be 'all', 'morph', 'env', 'eco', 'geog', 'TR', 'innovationIndex'.
#' 
#'	'TR' refers to tip rates of trait evolution, 'innovationIndex' refers to a measure of net
#'	innovation of species compared to the inferred squamate ancestral state. 
#'	 
#'	See this page of our website for a detailed explanation of the traits dataset. 
#' 
#'	Note that this dataset contains trait data for species that are not in the phylogeny.
#'	The column \code{inTree} defines which taxa are included in the tree.
#' @examples 
#' \dontrun{
#' # get the full dataset with no subsetting
#' head(squamata_traits())
#' 
#' # get the full dataset for Crotalus
#' head(squamata_traits(genus = 'Crotalus'))
#'  
#' # return only morphological traits
#' head(squamata_traits(genus = 'Crotalus', type = 'morph'))
#' 
#' # return tip rates for all of Iguania
#' head(squamata_traits(subclade = 'Iguania', type = 'TR'))
#' 
#' }
#' @rdname squamata_traits
#' @export 




squamata_traits <- function(species = NULL, genus = NULL, family = NULL, subclade = NULL, type = 'all') {
	
	if (sum(c(is.null(species), is.null(genus), is.null(family), is.null(subclade))) < 3) {
		stop("You cannot subset by more than one taxonomic rank.")
	}

	subCladeOptions <- c('Acrodonta', 'Alethinophidia', 'Amphisbaenia', 'Anguiformes', 'Caenophidia', 'Colubriformes', 'Colubrinae', 'Colubroidea', 'Colubroides', 'Dipsadinae', 'Episquamata', 'Gekkota', 'Iguania', 'Lacertoidea', 'Pleurodonta', 'Scincoidea', 'Scolecophidia', 'Serpentes', 'Teioidea', 'Toxicofera', 'tropicalDipsadines', 'Unidentata')
	
	if (!is.null(subclade)) {
		subclade <- match.arg(subclade, subCladeOptions, several.ok = TRUE)
	}

	typeOptions <- c('all', 'morph', 'env', 'eco', 'geog', 'TR', 'innovationIndex')
	
	type <- match.arg(type, typeOptions)
	if (length(type) != 1) stop(paste0("type must be only one of the following: ", paste0(typeOptions, collapse = ', ')))

	dat <- .getFile(target = 'squamdat')

	taxCols <- c('treename', 'genus', 'subfamily', 'family', 'inTree')
	divRateCols <- c('dr' ,'meanImputedDR', 'clads', 'meanImputedCLADS', 'bamm')
	morphCols <- c('mass', 'completeSVL', 'elongationIndex', 'numberDigits', 'numberLimbs', 'numberPresacralVert', 'numberCaudalVert', 'skullPC1', 'skullPC2', 'mesokinesis', 'metakinesis', 'hypokinesis', 'hyperkinesis', 'combinedKinesis')
	envCols <- c('centroidLong', 'centroidLat', 'minLong', 'maxLong', 'minLat', 'maxLat', 'bio1' ,'bio12', 'bio7', 'cmi', 'npp', 'elev', 'tri', 'rangeSize', 'climPC1', 'climPC2', 'climPC3', 'climPC4', 'climPC5', 'climPC6')
	ecoCols <- c('parity', 'prehensionMechanism', 'dietBreadth', 'dietPC1', 'dietPC2', 'foragingMode', 'chemosensory_index')
	geogCols <- grep('geog_', colnames(dat), value = TRUE)
	tipRateCols <- grep('rate', colnames(dat), value = TRUE, ignore.case = TRUE)
	ancDistCols <- grep('ancDist', colnames(dat), value = TRUE, ignore.case = TRUE)
		
	# if genus or family are provided for subsetting, get the associated species
	if (!is.null(genus) | !is.null(family) | !is.null(subclade)) {
		
		if (!is.null(genus)) {
			genus <- gsub('^\\s+|\\s+$', '', genus)
			species <- dat[which(dat$genus %in% genus), 'treename']
			species <- intersect(species, dat$treename)
		}
		
		if (!is.null(family)) {
			family <- gsub('^\\s+|\\s+$', '', family)
			species <- dat[which(dat$family %in% family), 'treename']
			species <- intersect(species, dat$treename)
		}

		if (!is.null(subclade)) {
			subclade <- gsub('^\\s+|\\s+$', '', subclade)
			species <- lapply(subclade, function(x) which(dat[, paste0('clade_', x)] == 1))
			species <- unique(unlist(species))
			species <- dat[species, 'treename']
			species <- intersect(species, dat$treename)
		}
	}
	
	# Now subset by species, whether it be derived from genus/family or requested directly
	if (!is.null(species)) {
		species <- gsub('^\\s+|\\s+$', '', species)
		species <- gsub('\\s+', '_', species)
	
		if (!all(species %in% dat$treename)) {
			notInDat <- setdiff(species, dat$treename)
			notInDat <- paste0('\t', notInDat)
			stop(paste0('The following taxa did not match any tip labels: \n', paste0(notInDat, collapse = '\n')))
		}
		
		dat <- dat[dat$treename %in% species, ]
	}
	
	# subset to the dataset of interest
	if (type == 'all') {
		dat <- dat[, c(taxCols, morphCols, envCols, ecoCols, geogCols, tipRateCols, ancDistCols)]
	} else if (type == 'morph') {
		dat <- dat[, c(taxCols, morphCols)]
	} else if (type == 'env') {
		dat <- dat[, c(taxCols, envCols)]
	} else if (type == 'eco') {
		dat <- dat[, c(taxCols, ecoCols)]
	} else if (type == 'geog') {
		dat <- dat[, c(taxCols, geogCols)]
	} else if (type == 'TR') {
		dat <- dat[, c(taxCols, tipRateCols)]
	} else if (type == 'innovationIndex') {
		dat <- dat[, c(taxCols, ancDistCols)]
	}
	
	
	rownames(dat) <- NULL
	return(dat)

}




