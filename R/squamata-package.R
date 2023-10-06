##' @title squamata R package
##'
##' @description An R package to facilitate access to squamate phylogeny and macroevolution resources.
##' 
##' @author Pascal O. Title, Sonal Singhal, Dan Rabosky and others
##' 
##' @references \url{https://github.com/ptitle/squamata-Rpackage}
##'
##' To cite the squamata package in publications, please cite the data generation paper:
##' Our empirical paper
##' as well as our resources paper:  
##' Citation of R package paper will go here.
##' 
##' @details
##' This R package provides access to a number of different types of phylogenies, as well as a compilation 
##' of species-level ecological, morphological and other attributes. 
##' 
##' \strong{Note that an internet connection is required for use of this package}, as data are downloaded rather than 
##' bundled with the package. When functions from this R package are used, a temporary directory is created, into which
##' data objects are downloaded. Objects are downloaded only once per R session, so repeated runs of the same function
##' will not lead to repeated identical downloads. However, the temporary folder is deleted at the end of the R session.
##' 
##' Major squamate subclades that can be used for subsetting include:
##' \code{'Acrodonta', 'Alethinophidia', 'Amphisbaenia', 'Anguiformes', 'Caenophidia', 'Colubriformes', 'Colubrinae', 
##' 'Colubroidea', 'Colubroides', 'Dipsadinae', 'Episquamata', 'Gekkota', 'Iguania', 'Lacertoidea', 'Pleurodonta', 
##' 'Scincoidea', 'Scolecophidia', 'Serpentes', 'Teioidea', 'Toxicofera', 'tropicalDipsadines', 'Unidentata'}.
##' 
##' @keywords internal
##'
##'
"_PACKAGE"
