#' @title Retrieve speciation rates

#' @description This function will retrieve species-specific speciation tip rates, including 
#' 	DR, BAMM and CLaDS.

#' @param species a vector of species names that the tree will be subset to, Default: NULL
#' @param genus a vector of genus names that the tree will be subset to, Default: NULL
#' @param family a vector of family names that the tree will be subset to, Default: NULL
#' @param subclade a vector of major squamate subclades that the tree will be subset to, Default: NULL

#' @return a data.frame with taxonomic information and speciation rates

#' @details 
#'	The list of major subclades that can be used for subsetting are listed under \code{\link{squamata}}.
#'	
#'	For DR and CLaDS, we generated 100 fully-sampled trees via phylogenetic imputation.
#'	The tip rates returned here are averages across those 100 trees for species that are present
#'	in our molecular phylogeny. In other words, we do not provide tip rates for those species that were 
#'	placed via imputation. For BAMM, we modeled diversification dynamics with family-level sampling fractions.
#' @examples 
#' \dontrun{
#' # get all rates for all speciation in the tree
#' head(squamata_speciationRates())
#' 
#' # get rates for Gekkkonid geckos
#' head(squamata_speciationRates(family = 'Gekkonidae'))
#' 
#' # get rates for all snakes
#' head(squamata_speciationRates(subclade = 'Serpentes'))
#' 
#' }
#' @rdname squamata_speciationRates
#' @export 




squamata_speciationRates <- function(species = NULL, genus = NULL, family = NULL, subclade = NULL) {
	
	if (sum(c(is.null(species), is.null(genus), is.null(family), is.null(subclade))) < 3) {
		stop("You cannot subset by more than one taxonomic rank.")
	}

	subCladeOptions <- c('Acrodonta', 'Alethinophidia', 'Amphisbaenia', 'Anguiformes', 'Caenophidia', 'Colubriformes', 'Colubrinae', 'Colubroidea', 'Colubroides', 'Dipsadinae', 'Episquamata', 'Gekkota', 'Iguania', 'Lacertoidea', 'Pleurodonta', 'Scincoidea', 'Scolecophidia', 'Serpentes', 'Teioidea', 'Toxicofera', 'tropicalDipsadines', 'Unidentata')
	
	if (!is.null(subclade)) {
		subclade <- match.arg(subclade, subCladeOptions, several.ok = TRUE)
	}

	dat <- .getFile(target = 'squamdat')
		
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
	
	dat <- dat[, c('treename', 'genus', 'subfamily', 'family', 'dr' ,'meanImputedDR', 'clads', 'meanImputedCLADS', 'bamm')]
	dat <- dat[!is.na(dat$dr), ]

	rownames(dat) <- NULL
	return(dat)

}