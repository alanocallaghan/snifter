#' @importFrom assertthat assert_that
#' @importFrom reticulate py_to_r py_to_r_wrapper import
#' @importFrom basilisk BasiliskEnvironment basiliskStart basiliskRun basiliskStop
NULL

#' @name snifter
#' 
#' @title snifter: fast interpolated t-SNE in R
#' 
#' An R package for running openTSNE's implementation of fast interpolated
#' t-SNE.
#' @references
#'  openTSNE: a modular Python library for t-SNE dimensionality reduction and
#'  embedding
#'  Pavlin G. Poličar, Martin Stražar, Blaž Zupan
#'  bioRxiv (2019) 731877; doi: \url{https://doi.org/10.1101/731877}
#'
#'  Fast interpolation-based t-SNE for improved visualization of single-cell
#'  RNA-seq data
#'  George C. Linderman, Manas Rachh, Jeremy G. Hoskins, Stefan Steinerberger,
#'  and Yuval Kluger
#'  Nature Methods 16, 243–245 (2019)
#'  doi: \url{https://doi.org/10.1038/s41592-018-0308-4}
#' 
#'  Accelerating t-SNE using Tree-Based Algorithms
#'  Laurens van der Maaten
#'  Journal of Machine Learning Research (2014)
#'  \url{http://jmlr.org/papers/v15/vandermaaten14a.html}
#' @docType package
NULL
