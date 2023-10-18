#' The test data set called 'real'
#'
#' @description These are first couple of thousand lines of a vcf file from real world data. It contains genotype information of 16 individuals of two bat species.
#'
#' @format A vcfR object.
#'
#' @examples
#' \dontrun{
#' data("real", package = "GenoPop")
#' head(vcf@gt)
#' }
#' @docType data
#' @usage data("real", package = "GenoPop")
#' @name real
#' @aliases real
#' @keywords datasets
#'
NULL

#' The test data set called 'mys'
#'
#' @description This is the real data set, but seperated by population (mys & dav) with missing data imputed (rf) and allele frequencies calculated.
#'
#' @format A myVcfR object
#'
#' @examples
#' \dontrun{
#' data("mys", package = "GenoPop")
#' head(mys@imp_gt)
#' }
#' @docType data
#' @usage data("mys", package = "GenoPop")
#' @name mys
#' @aliases mys
#' @keywords datasets
#'
NULL

#' The test data set called 'dav'
#'
#' @description This is the real data set, but seperated by population (mys & dav) with missing data imputed (rf) and allele frequencies calculated.
#'
#' @format A myVcfR object
#'
#' @examples
#' \dontrun{
#' data("dav", package = "GenoPop")
#' head(dav@imp_gt)
#' }
#' @docType data
#' @usage data("dav", package = "GenoPop")
#' @name dav
#' @aliases dav
#' @keywords datasets
#'
NULL

#' The test data set called 'sim'
#'
#' @description These are first couple of thousand lines of a vcf file from a simulated data set that was created by msprime. It also contains genotype information of 16 individuals.
#'
#' @format A vcfR object.
#'
#' @examples
#' \dontrun{
#' data("sim", package = "GenoPop")
#' head(vcf@gt)
#' }
#' @docType data
#' @usage data("sim", package = "GenoPop")
#' @name sim
#' @aliases sim
#' @keywords datasets
#'
NULL

