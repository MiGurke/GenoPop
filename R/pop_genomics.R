#' FixedSites
#'
#' Count the number of sites fixed for an alternative allele.
#'
#' @param object An S4 object of type myVcfR. Allele frequencies must be present.
#'
#' @return The number of fixed sites.
#'
#' @examples
#' data("real", package = "GenoPop")
#' vcf <- calculateAlleleFreqs(real, missing_data = "impute", method = "mean")
#' FixedSites(vcf)
#'
#' @export

FixedSites <- function(object) {
  # Extract the allele frequency table from the object
  allele_freq_table <- object@allele_freqs

  # Check if 'allele_freqs' slot in the object has data
  if (is.null(allele_freq_table) || nrow(allele_freq_table) == 0) {
    stop("No allele frequency data available in the object.")
  }

  # The number of fixed sites (excluding those fixed for the reference allele)
  num_fixed_sites <- 0
  # Check each site
  for (i in 1:nrow(allele_freq_table)) {
    site_freqs <- allele_freq_table[i, ]
    # Identify if there's an allele (excluding the reference allele '0') with frequency equal to 1
    if (any(site_freqs[-1] == 1)) {  # -1 to exclude the first column which represents the reference allele
      num_fixed_sites <- num_fixed_sites + 1
    }
  }
  return(num_fixed_sites)
}

#' PolymorphicSites
#'
#' Count the number of polymorphic sites in the data set.
#'
#' @param object An S4 object of type myVcfR. Allele frequencies must be present.
#'
#' @return The number of polymorphic sites.
#'
#' @examples
#' data("real", package = "GenoPop")
#' vcf <- calculateAlleleFreqs(real, missing_data = "impute", method = "mean")
#' PolymorphicSites(vcf)
#'
#' @export

PolymorphicSites <- function(object) {
  # Extract the allele frequency table from the object
  allele_freq_table <- object@allele_freqs

  # Check if 'allele_freqs' slot in the object has data
  if (is.null(allele_freq_table) || nrow(allele_freq_table) == 0) {
    stop("No allele frequency data available in the object.")
  }

  # The number of polymorphic sites
  num_polymorphic_sites <- 0

  # Check each site
  for (i in 1:nrow(allele_freq_table)) {
    site_freqs <- allele_freq_table[i, ]

    # Check if there's more than one allele present at the site
    # i.e., no allele has frequency 1 or 0 (excluding sites fixed for the reference allele)
    if (!any(site_freqs == 1) && !all(site_freqs == 0)) {
      num_polymorphic_sites <- num_polymorphic_sites + 1
    }
  }

  return(num_polymorphic_sites)
}

#' SingeltonSites
#'
#' Count the number of singelton sites in the data set. These are sites where a minor allele occurs only once in the sample.
#'
#' @param object An S4 object of type myVcfR. Allele frequencies must be present.
#'
#' @return The number of singelton sites.
#'
#' @examples
#' data("real", package = "GenoPop")
#' vcf <- calculateAlleleFreqs(real, missing_data = "impute", method = "mean")
#' SingeltonSites(vcf)
#'
#' @export

SingeltonSites <- function(object) {
  # Extract the genotype matrix, if there is an imputed version use that.
  sep <- object@sep_gt
  imp <- object@imp_gt
  if (is.null(imp) || nrow(imp) == 0) {
    genotype_matrix <- sep
  } else {
    genotype_matrix <- imp
  }

  # Check if we actually found some genotype data.
  if (is.null(genotype_matrix) || nrow(genotype_matrix) == 0) {
    stop("No genotype data available in the object.")
  }

  # The number of singleton sites
  num_singleton_sites <- 0

  # Check each site
  for (i in 1:nrow(genotype_matrix)) {
    site_genotypes <- genotype_matrix[i, ]

    # Calculate the allele counts at this site
    allele_counts <- table(site_genotypes)
    # A site is considered to have a singleton if any allele (except the reference allele, assumed to be '0')
    # is present exactly once in the entire sample, considering the ploidy level.
    # This means for diploid organisms, we're looking for alleles that occur twice (as each individual has two sets of chromosomes).
    # For haploid, we look for an allele count of 1, and so on.
    singleton_condition <- (allele_counts == 1) & (names(allele_counts) != "0")  # excluding the reference allele

    if (any(singleton_condition)) {
      num_singleton_sites <- num_singleton_sites + 1
    }
  }

  return(num_singleton_sites)
}
