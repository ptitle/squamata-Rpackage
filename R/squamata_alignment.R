#' @title Retrieve squamate sequence alignments
#'
#' @description This function will return the multiple sequence alignment used in the inference of the 
#' squamate phylogeny. This alignment can be subset by taxonomy and/or by gene region.
#' 
#' @param species a vector of species names that the alignment will be subset to, Default: NULL
#' @param genus a vector of genus names that the alignment will be subset to, Default: NULL
#' @param family a vector of family names that the alignment will be subset to, Default: NULL
#' @param subclade a vector of major squamate subclades that the tree will be subset to, Default: NULL
#' @param geneRegion names of gene regions to include, see details. All gene regions included if NULL
#' @param concatenate if TRUE (default), then alignment is concatenated, if FALSE, then separate gene-specific alignments are returned as a list. 
#' @param includeSphenodon Should Sphenodon_punctatus be included?, Default: FALSE
#' @param seqClass the alignment can be returned as class \code{DNAbin} or \code{DNAStringSet}, see details. \code{'auto'} means the function will preferentially pick \code{DNAStringSet} if available. 
#' @return if \code{concatenate = TRUE}, then a list with the alignment as the first item, and a table as the second item specifying how the alignment divides into the requested gene regions. If \code{concatenate = FALSE}, then a list of alignments, with one for each requested gene regions.

#' @details
#'	The list of major subclades that can be used for subsetting are listed under \code{\link{squamata}}.
#' 
#'	if \code{geneRegion = NULL}, then all will be returned. Otherwise some combination of the following can be specified:
#'	\code{'ADNP', 'AHR', 'AKAP9', 'AMEL', 'BACH1', 'BDNF', 'BHLHB2', 'BMP2', 'CAND1', 'CARD4', 'CILP', 'cmos', 
#'	'COI', 'CXCR4', 'CYTB', 'DLL1', 'ECEL', 'ENC1', 'FSHR', 'FSTL5', 'GALR1', 'GHSR', 'GPR37', 'HLCS', 'INHIBA', 
#'	'LRRN1', 'LZTSS1', 'MKL1', 'MLL3', 'MSH6', 'ND1', 'ND2', 'ND4', 'NGFB', 'NKTR', 'NTF-3', 'PDC', 'PNN', 'PRLR', 
#'	'PTGER4', 'PTPN', 'R35', 'RAG1_firsthalf', 'RAG1_secondhalf', 'RAG2', 'rRNA_12S', 'rRNA_16S', 'SINCAIP', 
#'	'SLC30A1', 'SLC8A1', 'SLC8A3', 'TRAF6', 'UBN1', 'VCPIP1', 'ZEB2', 'ZFP36L1'}
#'
#'	This function is set up to work with two classes for sequence data: \code{\link[ape]{DNAbin}} class from the 
#' 	R package ape, and \code{seqClass = \link[Biostrings]{DNAStringSet}} from the [bioconductor Biostrings R package](https://bioconductor.org/packages/release/bioc/html/Biostrings.html). 
#' 	If the latter is desired, then the R package Biostrings must be installed (available via bioconductor, not CRAN). 
#' 	Setting \code{seqClass = 'auto'} means the function will preferentially pick \code{DNAStringSet} because it is 
#' 	faster, and will fall back on \code{DNAbin} if the Biostrings package is not installed. 
#' 
#' @examples 
#' \dontrun{
#' 
#' # get the full concatenated alignment
#' squamata_alignment()
#'
#' # get the full concatenated alignment for Gekkota
#' squamata_alignment(subclade = 'Gekkota')
#' 
#' # get the alignment for genus Lerista, and only for gene regions ADNP and ND1
#' squamata_alignment(genus = 'Lerista', geneRegion = c('ADNP', 'ND1'), concatenate = FALSE)
#' 
#' # get the alignment for genus Lerista, and only for gene regions ADNP and ND1
#' ## this time get the concatenated alignment
#' squamata_alignment(genus = 'Lerista', geneRegion = c('ADNP', 'ND1'), concatenate = TRUE)
#' }
#' @seealso \code{\link[ape]{DNAbin}}, \code{\link[Biostrings]{DNAStringSet}}
#' @rdname squamata_alignment
#' @export 





squamata_alignment <- function(species = NULL, genus = NULL, family = NULL, subclade = NULL, geneRegion = NULL, concatenate = TRUE, includeSphenodon = FALSE, seqClass = 'auto') {

	if (sum(c(is.null(species), is.null(genus), is.null(family), is.null(subclade))) < 3) {
		stop("You cannot subset by more than one taxonomic rank.")
	}
	
	subCladeOptions <- c('Acrodonta', 'Alethinophidia', 'Amphisbaenia', 'Anguiformes', 'Caenophidia', 'Colubriformes', 'Colubrinae', 'Colubroidea', 'Colubroides', 'Dipsadinae', 'Episquamata', 'Gekkota', 'Iguania', 'Lacertoidea', 'Pleurodonta', 'Scincoidea', 'Scolecophidia', 'Serpentes', 'Teioidea', 'Toxicofera', 'tropicalDipsadines', 'Unidentata')
	
	if (!is.null(subclade)) {
		subclade <- match.arg(subclade, subCladeOptions, several.ok = TRUE)
	}	

	allgenes <- c('ADNP', 'AHR', 'AKAP9', 'AMEL', 'BACH1', 'BDNF', 'BHLHB2', 'BMP2', 'CAND1', 'CARD4', 'CILP', 'cmos', 'COI', 'CXCR4', 'CYTB', 'DLL1', 'ECEL', 'ENC1', 'FSHR', 'FSTL5', 'GALR1', 'GHSR', 'GPR37', 'HLCS', 'INHIBA', 'LRRN1', 'LZTSS1', 'MKL1', 'MLL3', 'MSH6', 'ND1', 'ND2', 'ND4', 'NGFB', 'NKTR', 'NTF-3', 'PDC', 'PNN', 'PRLR', 'PTGER4', 'PTPN', 'R35', 'RAG1_firsthalf', 'RAG1_secondhalf', 'RAG2', 'rRNA_12S', 'rRNA_16S', 'SINCAIP', 'SLC30A1', 'SLC8A1', 'SLC8A3', 'TRAF6', 'UBN1', 'VCPIP1', 'ZEB2', 'ZFP36L1')
	allgenes <- toupper(allgenes)
	
	if (!is.null(geneRegion)) {
		geneRegion <- toupper(geneRegion)
		if (!all(geneRegion %in% allgenes)) {
			geneNotRecognized <- setdiff(geneRegion, allgenes)
			geneNotRecognized <- paste0('\t', geneNotRecognized)
			stop(paste0('The following gene regions do not exactly match any existing options: \n', paste0(geneNotRecognized, collapse = '\n')))
		}
	}
	
	seqClass <- match.arg(seqClass, choices = c('auto', 'DNAStringSet', 'DNAbin'))
	if (length(seqClass) != 1) stop("seqClass must be one of 'auto', 'DNAbin' or 'DNAStringSet'.")
	
	if (!requireNamespace("Biostrings", quietly = TRUE) & seqClass == 'DNAStringSet') {
		warning('\tSwitching to class DNAbin -- class DNAStringSet is only possible if the Biostrings package is installed.')
		seqClass <- 'DNAbin'
	}

	if (seqClass == 'auto' & requireNamespace("Biostrings", quietly = TRUE)) {
		seqClass <- 'DNAStringSet'
	} else if (seqClass == 'auto' & !requireNamespace("Biostrings", quietly = TRUE)) {
		seqClass <- 'DNAbin'
	}

	# if no subsetting requested, then just get the alignment and return it
	if (is.null(genus) & is.null(family) & is.null(species)) {
		aln <- .getFile(target = 'alignment', seqClass = seqClass)
		if (!includeSphenodon) {
			if (inherits(aln, 'DNAbin')) {
				aln <- aln[setdiff(rownames(aln), 'Sphenodon_punctatus'), ]
			} else {
				aln <- aln[setdiff(names(aln), 'Sphenodon_punctatus')]
			}
		}

	} else {

		# get squamdat table
		tax <- .getFile(target = 'squamdat')
		tax <- tax[tax$inTree == 1, ]
	
		# if genus or family are provided for subsetting, get the associated species
		if (!is.null(genus) | !is.null(family) | !is.null(subclade)) {
			
			if (!is.null(genus)) {
				genus <- gsub('^\\s+|\\s+$', '', genus)
				genusNotFound <- setdiff(genus, tax$genus)
				if (length(genusNotFound) > 0) {
					stop(paste0('The following genera did not match our genus list: \n', paste0(paste0('\t', genusNotFound), collapse = '\n')))
				}
				species <- tax[which(tax$genus %in% genus), 'treename']
			}
			
			if (!is.null(family)) {
				family <- gsub('^\\s+|\\s+$', '', family)
				familyNotFound <- setdiff(family, tax$family)
				if (length(familyNotFound) > 0) {
					stop(paste0('The following families did not match our family list: \n', paste0(paste0('\t', familyNotFound), collapse = '\n')))
				}
				species <- tax[which(tax$family %in% family), 'treename']	
			}

			if (!is.null(subclade)) {
				subclade <- gsub('^\\s+|\\s+$', '', subclade)
				species <- lapply(subclade, function(x) which(tax[, paste0('clade_', x)] == 1))
				species <- unique(unlist(species))
				species <- tax[species, 'treename']
			}			
		}
		
		species <- gsub('^\\s+|\\s+$', '', species)
		species <- gsub('\\s+', '_', species)
	
		if (!all(species %in% tax$treename)) {
			notInTree <- setdiff(species, tax$treename)
			notInTree <- paste0('\t', notInTree)
			stop(paste0('The following taxa did not match any tip labels: \n', paste0(notInTree, collapse = '\n')))
		}
		
		if (includeSphenodon) {
			species <- c(species, 'Sphenodon_punctatus')
		}
		
		# get alignment and subset
		aln <- .getFile(target = 'alignment', seqClass = seqClass)
		
		if (inherits(aln, 'DNAbin')) {
			aln <- aln[species, ]
		} else {
			aln <- aln[species]
		}
	}
	
	geneTable <- .getFile(target = 'geneIndex')
	geneTable[, 1] <- toupper(geneTable[, 1])

	# if gene regions specified, subset to those regions
	if (!is.null(geneRegion)) {
		genes <- geneTable[geneTable[,1] %in% geneRegion, ]
		
		if (inherits(aln, 'DNAbin')) {
			geneInd <- apply(genes, 1, function(x) seq.int(from = x[2], to = x[3]))
			geneSpans <- lengths(geneInd)
			
			if (concatenate) {
				geneInd <- unlist(geneInd, use.names = FALSE)
				aln <- aln[, geneInd]
			} else {
				aln <- lapply(geneInd, function(x) aln[, x])
				names(aln) <- genes[, 1]
			}
			
		} else if (inherits(aln, 'DNAStringSet')) {
			aln <- apply(genes, 1, function(x) {
				Biostrings::subseq(aln, start = as.integer(x[2]), end = as.integer(x[3]))
			})
			names(aln) <- genes[, 1]
			geneSpans <- sapply(aln, function(x) length(x[[1]]))
			
			if (concatenate) {
				alnCat <- do.call(Biostrings::xscat, aln)
				names(alnCat) <- names(aln[[1]])
				aln <- alnCat
				rm(alnCat)
			}
		}
		
		# create table of index ranges for each gene
		geneRanges <- genes
		geneRanges[, 2] <- geneRanges[, 3] <- 0
		geneRanges[1, 2] <- 1
		geneRanges[1, 3] <- geneSpans[1]
		if (length(geneSpans) > 1) {
			for (i in 2:length(geneSpans)) {
				geneRanges[i, 2] <- geneRanges[(i-1), 3] + 1
				geneRanges[i, 3] <- geneRanges[i, 2] + (geneSpans[i] - 1)
			}
		}
		
	} else {
		geneRanges <- geneTable
	}
	
	if (concatenate) {
		return(list('alignment' = aln, 'geneRegions' = geneRanges))
	} else {
		return(aln)
	}
}

