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


#' Heterozygosity
#'
#' This function calculates the observed and expected heterozygosity of a population.
#'
#' @param object An S4 object of type myVcfR. Allele frequencies must be present.
#'
#' @return A list containing observed and expected heterozygosity.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' Heterozygosity(mys)
#'
#' @export

Heterozygosity <- function(object) {
  # Extract the genotype matrix and allele frequencies
  sep <- object@sep_gt
  imp <- object@imp_gt
  genotype_matrix <- if (is.null(imp) || nrow(imp) == 0) sep else imp
  allele_freqs <- object@allele_freqs

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

  # Calculate expected heterozygosity
  He_per_locus <- apply(allele_freqs, 1, function(p) {
    1 - sum(p^2)
  })
  # Average over all loci
  He <- mean(He_per_locus)

  return(list(observed_heterozygosity = Ho, expected_heterozygosity = He))
}

#' WindowedPi
#'
#' Calculate the average number of nucleotide differences per site between two sequences, for windows of specified length. The algorithm used for this is equivalent to the one used in vcftools --window-pi (https://vcftools.sourceforge.net/man_latest.html).
#'
#' @param object An S4 object of type myVcfR.
#' @param window_size The size of the window for which Pi is calculated. (Default = 1000)
#' @param step_size The size of the step in between windows. (Default = 0)
#'
#' @return The nucleotide diversity (pi) per window in data frame.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' WindowedPi(mys)
#'
#' @export

WindowedPi <- function(object, window_size = 1000, step_size = 0) {
  # Extract the positions and the chromosomes
  variant_positions <- as.numeric(object@fix[, 2])  # Assuming the second column contains positions
  chromosomes <- object@fix[, 1]  # Assuming the first column contains chromosome names

  # Extract the genotype matrix and allele frequencies
  sep <- object@sep_gt
  imp <- object@imp_gt
  genotype_matrix <- if (is.null(imp) || nrow(imp) == 0) sep else imp

  # Check if we actually found some genotype data
  if (is.null(genotype_matrix) || nrow(genotype_matrix) == 0) {
    stop("No genotype data available in the object.")
  }

  # Determine the number of chromsomes (each individual has two chroms for diploid organisms)
  num_chroms <- ncol(genotype_matrix)

  # Prepare a list to hold results
  results <- list()

  # We will calculate pi for each chromosome separately
  unique_chromosomes <- unique(chromosomes)
  for (chr in unique_chromosomes) {
    # Get indices of the variants on this chromosome
    chr_indices <- which(chromosomes == chr)
    chr_positions <- variant_positions[chr_indices]

    # Calculate the number of windows based on the step size and the positions of the variants
    num_windows <- ceiling((max(chr_positions) - min(chr_positions)) / (step_size + window_size))
    print(num_windows)
    # Prepare a vector to store pi values for each window
    windowed_pi <- numeric(num_windows)

    # Loop over each window and calculate pi
    for (i in seq_len(num_windows)) {
      start_pos <- (i - 1) * (window_size + step_size)
      end_pos <- start_pos + window_size
      print(c(start_pos,end_pos))

      # Get indices of variants within this window
      window_indices <- which(chr_positions >= start_pos & chr_positions <= end_pos)
      print(length(window_indices))
      # If there are no variants in this window, skip the calculations
      if (length(window_indices) == 0) {
        windowed_pi[i] <- NA  # No data for this window
        next
      }

      # Select the genotype data for these positions
      window_genotypes <- genotype_matrix[chr_indices[window_indices], ]

      # Calculate nucleotide diversity (pi) for this window
      N_mismatches_total <- 0
      N_comparisons_total <- 0

      for (site_index in seq_len(nrow(window_genotypes))) {
        site_genotypes <- window_genotypes[site_index, ]
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

      # Pi calculation for the window, including monomorphic sites
      N_monomorphic_sites <- window_size - length(window_indices)
      N_pairs_total <- N_comparisons_total + (N_monomorphic_sites * num_chroms * (num_chroms - 1))
      windowed_pi[i] <- N_mismatches_total / N_pairs_total
    }

    # Add the results for this chromosome to the results list
    results[[chr]] <- windowed_pi
  }

  return(results)
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

TajimasD <- function(object) {
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

  # Calculate average number of nucleotide differences (pi)
  pi <- # ... [you need to implement the calculation for pi here, based on pairwise differences]

    # Calculate constants based on sample size
  n <- ncol(genotype_matrix) / 2  # assuming diploid organisms
  a1 <- sum(1 / (1:(n - 1)))
  a2 <- sum(1 / (1:(n - 1))^2)
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1 <- b1 - (1 / a1)
  c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)

  # Calculate Tajima's D
  D <- (pi - (S / a1)) / sqrt((e1 * S + e2 * S * (S - 1)) / (a1^2 + a2))

  return(D)
}



