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

#' SegregatingSites
#'
#' Count the number of segregating sites in the data set.
#'
#' @param object An S4 object of type myVcfR. Allele frequencies must be present.
#'
#' @return The number of segregating sites.
#'
#' @examples
#' data("real", package = "GenoPop")
#' vcf <- calculateAlleleFreqs(real, missing_data = "impute", method = "mean")
#' SegregatingSites(vcf)
#'
#' @export

SegregatingSites <- function(object) {
  # Extract the allele frequency table from the object
  allele_freq_table <- object@allele_freqs

  # Check if 'allele_freqs' slot in the object has data
  if (is.null(allele_freq_table) || nrow(allele_freq_table) == 0) {
    stop("No allele frequency data available in the object.")
  }

  # The number of segregating sites
  num_segregating_sites <- 0

  # Check each site
  for (i in 1:nrow(allele_freq_table)) {
    site_freqs <- allele_freq_table[i, ]

    # Check if there's more than one allele present at the site
    # i.e., no allele has frequency 1 or 0 (excluding sites fixed for the reference allele)
    if (!any(site_freqs == 1) && !all(site_freqs == 0)) {
      num_segregating_sites <- num_segregating_sites + 1
    }
  }

  return(num_segregating_sites)
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
    singleton_condition <- (allele_counts == 1) & (names(allele_counts) != "0")  # excluding the reference allele

    if (any(singleton_condition)) {
      num_singleton_sites <- num_singleton_sites + 1
    }
  }

  return(num_singleton_sites)
}

#' PrivateAlleles
#'
#' Count the number of private alleles in each of several populations. These are alleles found exclusively in one population and not in the others.
#'
#' @param objects A list of S4 objects of type myVcfR, each representing a different population. Allele frequencies must be present.
#'
#' @return A list containing the number of private alleles for each population.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' data("dav", package = "GenoPop")
#' PrivateAlleles(list(mys, dav))
#'
#' @export

PrivateAlleles <- function(objects) {
  # Ensure that objects is indeed a list of multiple populations
  if (length(objects) < 2) {
    stop("Please provide two or more population objects.")
  }

  # Initialize a list to store positions, chromosomes, and genotype matrices from all populations
  all_positions <- list()
  all_chromosomes <- list()
  genotype_matrices <- list()

  # Extract data from each object
  for (i in seq_along(objects)) {
    object <- objects[[i]]

    # Extract and store genotype matrix
    sep <- object@sep_gt
    imp <- object@imp_gt
    genotype_matrices[[i]] <- if (is.null(imp) || nrow(imp) == 0) sep else imp

    # Extract and store positions and chromosomes
    all_positions[[i]] <- object@fix[,2]  # assuming the 2nd column is 'POS'
    all_chromosomes[[i]] <- object@fix[,1]  # assuming the 1st column is 'CHROM'
  }

  # Check if we actually found some genotype data in each object
  for (genotype_matrix in genotype_matrices) {
    if (is.null(genotype_matrix) || nrow(genotype_matrix) == 0) {
      stop("No genotype data available in one or more of the objects.")
    }
  }

  # Combine all chromosomes and positions to identify unique and shared sites
  combined_chromosomes_positions <- list()

  for (i in seq_along(all_positions)) {
    combined_chromosomes_positions[[i]] <- paste(all_chromosomes[[i]], all_positions[[i]], sep = ":")
  }

  unique_combined_sites <- Reduce(union, combined_chromosomes_positions)

  # Initialize a list to keep track of private alleles count for each population
  private_alleles_counts <- rep(0, length(objects))

  # Check each site in the combined set
  for (site in unique_combined_sites) {
    # Split chromosome and position
    site_info <- strsplit(site, ":")[[1]]
    chromosome <- site_info[1]
    position <- as.numeric(site_info[2])

    # Store unique alleles from all populations for the current site
    unique_alleles_all_pops <- lapply(seq_along(all_positions), function(i) {
      # Find indices where chromosome and position match
      indices <- which(all_chromosomes[[i]] == chromosome & all_positions[[i]] == position)
      if (length(indices) > 0) {
        unique(genotype_matrices[[i]][indices, ])
        } else {
          character(0)
        }
    })

    # Count private alleles for each population at this site
    for (i in seq_along(unique_alleles_all_pops)) {
      other_pops <- unique_alleles_all_pops[-i]  # all populations except the current one
      combined_other_pops <- Reduce(union, other_pops)  # combine alleles from other populations

      # Identifying private alleles
      private_alleles <- setdiff(unique_alleles_all_pops[[i]], combined_other_pops)

      # Exclude the reference allele, assumed to be '0'
      private_alleles <- private_alleles[private_alleles != "0"]

      # Update counts
      private_alleles_counts[i] <- private_alleles_counts[i] + length(private_alleles)
    }
  }

  # Prepare the results to be returned, associating each count with the corresponding population
  names(private_alleles_counts) <- paste0('pop', seq_along(private_alleles_counts))
  return(private_alleles_counts)
}

#' Observed Heterozygosity
#'
#' This function calculates the observed heterozygosity in a population.
#'
#' @param object An S4 object of type myVcfR.
#'
#' @return Observed heterozygosity.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' ObservedHeterozygosity(mys)
#'
#' @export

ObservedHeterozygosity <- function(object) {
  # Extract the genotype matrix
  sep <- object@sep_gt
  imp <- object@imp_gt
  genotype_matrix <- if (is.null(imp) || nrow(imp) == 0) sep else imp

  # Check if we actually found some genotype data
  if (is.null(genotype_matrix) || nrow(genotype_matrix) == 0) {
    stop("No genotype data available in the object.")
  }

  # Calculate observed heterozygosity
  num_individuals <- ncol(genotype_matrix) / 2  # assuming diploid organisms
  heterozygotes <- 0

  # Iterate over individuals by stepping 2 columns at a time
  for (indiv_col in seq(1, ncol(genotype_matrix), by = 2)) {
    individual_genotypes <- genotype_matrix[, c(indiv_col, indiv_col + 1)]
    heterozygotes <- heterozygotes + sum(apply(individual_genotypes, 1, function(locus) {
      locus[1] != locus[2]
    }))
  }

  Ho <- heterozygotes / (num_individuals * nrow(genotype_matrix))
  return(Ho)
}

#' Expected Heterozygosity
#'
#' This function calculates the expected heterozygosity of a population.
#'
#' @param object An S4 object of type myVcfR. Allele frequencies must be present.
#'
#' @return Expected heterozygosity.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' ExpectedHeterozygosity(mys)
#'
#' @export

ExpectedHeterozygosity <- function(object) {
  # Extract allele frequencies
  allele_freqs <- object@allele_freqs

  # Check if allele frequencies are present
  if (is.null(allele_freqs) || nrow(allele_freqs) == 0) {
    stop("No allele frequency data available in the object.")
  }

  # Calculate expected heterozygosity per locus
  He_per_locus <- apply(allele_freqs, 1, function(p) {
    1 - sum(p^2)
  })

  # Average over all loci
  He <- mean(He_per_locus)

  return(He)
}

#' Pi
#'
#' Calculate the average number of nucleotide differences per site between two sequences. The formula used for this is equivalent to the one used in vcftools --window-pi (https://vcftools.sourceforge.net/man_latest.html).
#'
#' @param object An S4 object of type myVcfR.
#' @param window_size The size of the window for which Pi is calculated. (Default = 1000)
#' @param step_size The size of the step in between windows. (Default = 0)
#'
#' @return The nucleotide diversity (pi) per window in data frame.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' Pi(mys)
#'
#' @export

Pi <- function(object, seq_length) {
  # Assuming the object has slots for the genotype and other necessary data.
  # Extract the genotype matrix from the object
  sep <- object@sep_gt
  imp <- object@imp_gt
  genotype_matrix <- if (is.null(imp) || nrow(imp) == 0) sep else imp

  N_mismatches_total <- 0
  N_comparisons_total <- 0
  num_chroms <- ncol(genotype_matrix)

  # Check if there are any variants in this window
  if (nrow(genotype_matrix) == 0) {
    # If no variants, the Pi is 0 for this window, considering only monomorphic sites
    return(0)
  }

  for (site_index in seq_len(nrow(genotype_matrix))) {
    site_genotypes <- genotype_matrix[site_index, ]
    site_allele_freqs <- table(site_genotypes)  # Allele frequencies at this site

    # Total alleles at this site (not missing)
    N_non_missing_chr <- sum(site_allele_freqs)

    # Number of actual nucleotide differences (mismatches) for the site
    N_site_mismatches <- 0
    for (allele_count in site_allele_freqs) {
      N_site_mismatches <- N_site_mismatches + (allele_count * (N_non_missing_chr - allele_count))
    }

    N_mismatches_total <- N_mismatches_total + N_site_mismatches
    N_comparisons_total <- N_comparisons_total + (N_non_missing_chr * (N_non_missing_chr - 1))
  }

  # Including monomorphic sites
  N_monomorphic_sites <- seq_length - nrow(genotype_matrix)
  N_pairs_total <- N_comparisons_total + (N_monomorphic_sites * num_chroms * (num_chroms - 1))

  # Pi calculation for the window
  pi_value <- N_mismatches_total / N_pairs_total

  return(pi_value)
}

#' Tajima's D
#'
#' Calculate Tajima's D statistic for a given dataset, a measure for neutrality.
#'
#' @param object An S4 object of type myVcfR. Allele frequencies and genotype matrix must be present.
#'
#' @return Tajima's D value.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' TajimasD(mys)
#'
#' @export

TajimasD <- function(object, seq_length) {
  # Extract the genotype matrix
  sep <- object@sep_gt
  imp <- object@imp_gt
  genotype_matrix <- if (is.null(imp) || nrow(imp) == 0) sep else imp

  # Check if there's actual genotype data
  if (is.null(genotype_matrix) || nrow(genotype_matrix) == 0) {
    stop("No genotype data available in the object.")
  }

  # Calculate the number of segregating sites
  S <- SegregatingSites(object)

  # Calculate number of nucleotide differences (pi)
  pi <- Pi(object, seq_length)

  # Calculate constants based on sample size n
  n <- ncol(genotype_matrix)
  # Step 1: Calculate a1 and a2, which are summations used in the denominator of Tajima's D
  # These summations are over 1/i and 1/i^2, respectively, from i = 1 to (n-1)
  i_values <- 1:(n-1)
  a1 <- sum(1 / i_values)
  a2 <- sum(1 / (i_values^2))
  # Step 2: Calculate b1 and b2, constants used for the estimation of the variance
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- (2 * (n^2 + n + 3)) / (9 * n * (n - 1))
  # Step 3: Calculate c1 and c2, which scale b1 and b2 by a1 and a2 to get the actual variance components
  c1 <- b1 - (1 / a1)
  c2 <- b2 - ((n + 2) / (a1 * n)) + (a2 / (a1^2))
  # Step 4: Calculate e1 and e2, which are used in the denominator of Tajima's D
  # These represent the expectations of the variance and covariance, respectively
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  # Step 5: Calculate the denominator of Tajima's D, which is the square root of the variance of the difference between pi and S/a1
  denominator <- sqrt((e1 * S) + (e2 * S * (S - 1)))
  # Step 6: Calculate Tajima's D
  # This is the difference between pi and the mean number of pairwise differences (S/a1), divided by the standard deviation (denominator)
  D <- (pi - (S / a1)) / denominator

  return(D)
}



