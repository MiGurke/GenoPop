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
#' Function to calculate the number of private alleles in  two populations.
#'
#' @param object A myVcfR object.
#' @param pop_assignments A named vector. Elements are the population names and names are the individual name.
#'
#' @return A list containing the number of private alleles for each population.
#'
#' @examples
#' mys1 <- c("8449", "8128", "8779")
#' mys2 <- c("8816", "8823", "8157")
#'
#' individuals <- c(mys1, mys2)
#' pop_names <- c(rep("mys1", length(mys1)), rep("mys2", length(mys2)))
#' pop_assignments <- setNames(pop_names, individuals)
#'
#' data("mys", package = "GenoPop")
#' PrivateAlleles(mys, pop_assignments)
#'
#' @export

PrivateAlleles <- function(object, pop_assignments) {

  # Use the pop_assignments vector to separate the populations
  separated_pops <- seperateByPopulations(object, pop_assignments, rm_ref_alleles = TRUE)
  pop1 <- separated_pops[[1]]
  pop2 <- separated_pops[[2]]

  # Extract and store positions and chromosomes for both populations
  positions_pop1 <- pop1@fix[,2]  # assuming the 2nd column is 'POS'
  chromosomes_pop1 <- pop1@fix[,1]  # assuming the 1st column is 'CHROM'
  positions_pop2 <- pop2@fix[,2]  # assuming the 2nd column is 'POS'
  chromosomes_pop2 <- pop2@fix[,1]  # assuming the 1st column is 'CHROM'

  # Combine chromosomes and positions into a unique site identifier for both populations
  sites_pop1 <- paste(chromosomes_pop1, positions_pop1, sep=":")
  sites_pop2 <- paste(chromosomes_pop2, positions_pop2, sep=":")

  # Identify private alleles by finding unique sites in each population
  private_sites_pop1 <- setdiff(sites_pop1, sites_pop2)
  private_sites_pop2 <- setdiff(sites_pop2, sites_pop1)

  # Return the count of private alleles as a named vector
  private_alleles_count <- c(pop1 = length(private_sites_pop1), pop2 = length(private_sites_pop2))
  return(private_alleles_count)
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
#' @param seq_length Length of the sequence in number of bases. Must be provided to accurately work with all monomorphic sites, including those monomorphic for the reference, which are generally not included in a vcf.
#'
#' @return The nucleotide diversity (pi) per window in data frame.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' Pi(mys, 265392)
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
#' @param object A S4 object of type myVcfR. Allele frequencies and genotype matrix must be present.
#' @param seq_length Length of the sequence in number of bases. Must be provided to accurately work with all monomorphic sites, including those monomorphic for the reference, which are generally not included in a vcf.
#'
#' @return Tajima's D value.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' TajimasD(mys, 265392)
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
  pi <- Pi(object, seq_length) * seq_length

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

#' WattersonsTheta
#'
#' Calculate Watterson's thea, a measure for neutrality, from an myVcfR object. The metric will be normalized by the sequence length to make it comparable between data sets.
#'
#' @param object A S4 object of type myVcfR. Allele frequencies and genotype matrix must be present.
#' @param seq_length The length of the sequence in the data set.
#'
#' @return Watterson's theta value.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' WattersonsTheta(mys, 265392)
#'
#' @export

WattersonsTheta <- function(object, seq_length) {
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
  # Calculate constants based on sample size n
  n <- ncol(genotype_matrix)
  # The sum of the harmonic series until n-1
  i_values <- 1:(n-1)
  a1 <- sum(1 / i_values)
  WattersonsTheta <- S / a1
  return(WattersonsTheta / seq_length)
}

#' Fst
#'
#' Calculate the mean or weighted (by number of non missing chromsomes) fixiation index (Fst) from two populations in a list of myVcfR objects using the method of Weir and Cockerham (1984).
#'
#' @param object A myVcfR object.
#' @param pop_assignments A named vector. Elements are the population names and names are the individual name.
#' @param weighted Logical, wether weighted Fst or mean Fst is returned. (Default = FALSE (mean Fst is returned))
#'
#' @return Fst value.
#'
#' @examples
#' mys1 <- c("8449", "8128", "8779")
#' mys2 <- c("8816", "8823", "8157")
#'
#' individuals <- c(mys1, mys2)
#' pop_names <- c(rep("mys1", length(mys1)), rep("mys2", length(mys2)))
#' pop_assignments <- setNames(pop_names, individuals)
#'
#' data("mys", package = "GenoPop")
#' Fst(mys, pop_assignments)
#'
#' @export

Fst <- function(object, pop_assignments, weighted = FALSE) {
  # Check for necessary components in the object
  if (!all(c("sep_gt", "allele_freqs") %in% slotNames(object))) {
    stop("The object does not have the necessary components.")
  }

  # Separate the populations
  separated_pops <- seperateByPopulations(object, pop_assignments, rm_ref_alleles = FALSE)
  pop1 <- separated_pops[[1]]
  pop2 <- separated_pops[[2]]

  # Extract the genotype matrices
  sep1 <- pop1@sep_gt
  imp1 <- pop1@imp_gt
  genotype_matrix1 <- if (is.null(imp1) || nrow(imp1) == 0) sep1 else imp1
  sep2 <- pop2@sep_gt
  imp2 <- pop2@imp_gt
  genotype_matrix2 <- if (is.null(imp2) || nrow(imp2) == 0) sep2 else imp2

  # Extract allele frequencies
  allele_freqs1 <- pop1@allele_freqs
  allele_freqs2 <- pop2@allele_freqs

  # Ensure both matrices have the same number of columns by adding zeros if necessary
  # This can happen if one population has f.e. only the reference allele for all positions.
  # Especially when used in small windows for the windowed metric calculation.
  if (ncol(allele_freqs1) != ncol(allele_freqs2)) {
    max_cols <- max(ncol(allele_freqs1), ncol(allele_freqs2))
    if (ncol(allele_freqs1) < max_cols) {
      missing_cols <- setdiff(0:(max_cols - 1), colnames(allele_freqs1))
      for (col in missing_cols) {
        allele_freqs1 <- cbind(allele_freqs1, 0)
        colnames(allele_freqs1)[ncol(allele_freqs1)] <- as.character(col)
      }
    }
    if (ncol(allele_freqs2) < max_cols) {
      missing_cols <- setdiff(0:(max_cols - 1), colnames(allele_freqs2))
      for (col in missing_cols) {
        allele_freqs2 <- cbind(allele_freqs2, 0)
        colnames(allele_freqs2)[ncol(allele_freqs2)] <- as.character(col)
      }
    }
  }

  # Calculate average allele frequencies across populations
  mean_freqs <- (allele_freqs1 + allele_freqs2) / 2

  if (nrow(genotype_matrix1) < 5){
    print(mean_freqs)
  }

  # Determine the number of loci
  num_loci <- nrow(allele_freqs1)

  Fsts <- c()  # Initialize a vector to hold Fst values for each locus
  weights <- c() # Initialize a vector to hold the weights for each locus

  # Loop through each locus and calculate components of genetic variance
  for (i in 1:num_loci) {
    # Get the number of non-missing chromosomes
    n1 <- length(genotype_matrix1[i,genotype_matrix1[i,] != "."])
    n2 <- length(genotype_matrix2[i,genotype_matrix2[i,] != "."])

    # Get all the allele frequencies
    p1 <- as.numeric(allele_freqs1[i, 1])  # frequency of allele 1 in population 1
    p2 <- as.numeric(allele_freqs2[i, 1])  # frequency of allele 1 in population 2
    q1 <- 1 - p1  # frequency of allele 2 in population 1
    q2 <- 1 - p2  # frequency of allele 2 in population 2

    p_bar <- as.numeric(mean_freqs[i, 1])  # average frequency of allele 1
    q_bar <- 1 - p_bar  # average frequency of allele 2

    # Calculate the expected homozygosity within each population for this locus
    h1 <- p1^2 + q1^2
    h2 <- p2^2 + q2^2

    # Calculate the sample variances for this locus
    s1 <- p1 * q1 / (n1 - 1)  # sample variance in population 1
    s2 <- p2 * q2 / (n2 - 1)  # sample variance in population 2

    # Calculate components of genetic variance for this locus
    a <- ((p1 - p_bar)^2 + (p2 - p_bar)^2) / 2 - (s1 + s2) / 2  # among populations
    b <- ((s1 + s2) / 2) - ((1 / (2 * n1)) * (p1 * (1 - p1) + p2 * (1 - p2)) / 2)  # between individuals within populations
    c <- (1 / 2) * ((p1 * (1 - p1) + p2 * (1 - p2)) / 2)  # within individuals

    # Calculate Fst for this locus
    Fst_locus <- a / (a + b + c)
    # Consider negative Fst's to be 0, because they are an artifact of uneven sample sizes.
    if (!is.na(Fst_locus) && Fst_locus < 0) {
      Fst_locus <- 0
    }
    if (weighted) {
      weights <- c(weights, (n1+n2))
    }
    #print(Fst_locus)
    Fsts <- c(Fsts, Fst_locus)
  }

  if (weighted) {
    fst_weighted_products <- Fsts * weights
    sum_weighted_products <- sum(fst_weighted_products, na.rm = TRUE)
    sum_weights <- sum(weights, na.rm = TRUE)
    Fst_final <- sum_weighted_products / sum_weights
  } else {
    Fst_final <- mean(Fsts, na.rm = TRUE)
  }

  return(Fst_final)
}

#' OneDimSFS
#'
#' Calculate a one dimensional site frequency spectrum from an myVcfR object.
#'
#' @param object A S4 object of type myVcfR. Allele frequencies and genotype matrix must be present.
#' @param folded Logical, deciding if folded (TRUE) or unfolded (FALSE) SFS is returned. For the unfolded it is assumed that the genotype "0" represents the ancestral state in the data. (Default is unfolded (FALSE).)
#'
#' @return Site frequency spectrum as a named vector
#'
#' @examples
#' data("mys", package = "GenoPop")
#' OneDimSFS(mys, folded = FALSE)
#'
#' @export

OneDimSFS <- function(object, folded = FALSE) {
  # Extract the genotype matrix
  sep <- object@sep_gt
  imp <- object@imp_gt
  genotype_matrix <- if (is.null(imp) || nrow(imp) == 0) sep else imp

  # Replace '.' with NA for missing data
  genotype_matrix[genotype_matrix == "."] <- NA

  # Convert the entire matrix to numeric in a vectorized manner
  genotype_matrix_numeric <- as.matrix(apply(genotype_matrix, c(1, 2), as.numeric))

  # Number of individuals (assuming diploid, so 2 columns per individual)
  num_individuals <- ncol(genotype_matrix) / 2

  # Initialize a vector to hold the site frequency spectrum
  sfs <- numeric(num_individuals + 1)  # frequencies from 0 to num_individuals

  # Iterate over the sites in the genotype matrix
  for (i in 1:nrow(genotype_matrix_numeric)) {
    site_data <- genotype_matrix_numeric[i, ]

    # Exclude missing data for this site
    valid_data <- site_data[!is.na(site_data)]

    # Count the number of derived alleles (assuming '1' is the derived state)
    derived_count <- sum(valid_data)

    # Calculate the minor allele count for folded SFS
    if (folded) {
      allele_count <- min(derived_count, length(valid_data) - derived_count)
    } else {
      allele_count <- derived_count
    }

    # Adjust the total number of alleles based on missing data
    total_alleles_at_site <- length(valid_data)

    # Skip sites with no valid data
    if (total_alleles_at_site == 0) next

    # Calculate the frequency, adjusting for the varying number of valid alleles
    freq_index <- allele_count * (num_individuals) / total_alleles_at_site

    # Round to the nearest integer to get the discrete frequency category
    freq_category <- round(freq_index) + 1  # R is 1-indexed

    # Update the SFS
    sfs[freq_category] <- sfs[freq_category] + 1
  }

  # Name the vector elements for clearer interpretation
  names(sfs) <- 0:num_individuals
  # For a folded SFS, remove the redundant second half of the vector
  if (folded) {
    # Determine the midpoint of the vector
    midpoint <- ceiling((num_individuals + 1) / 2)
    # Keep only up to the midpoint (inclusive)
    sfs <- sfs[1:midpoint]
  }

  return(sfs)
}

#' TwoDimSFS
#'
#' Calculate a two-dimensional site frequency spectrum from a list of two myVcfR objects.
#'
#' @param objects A list of two S4 objects of type myVcfR. Allele frequencies and genotype matrices must be present.
#' @param folded Logical, deciding if folded (TRUE) or unfolded (FALSE) SFS is returned. (Default is unfolded (FALSE).)
#'
#' @return Two-dimensional site frequency spectrum as a matrix
#'
#' @examples
#' data("mys", package = "GenoPop")
#' data("dav", package = "GenoPop")
#' TwoDimSFS(list(mys, dav), folded = TRUE)
#'
#' @export

TwoDimSFS <- function(objects, folded = FALSE) {
  if (length(objects) != 2) {
    stop("Please provide a list of two myVcfR objects.")
  }

  # Extract the variant information
  variant_info1 <- objects[[1]]@fix
  variant_info2 <- objects[[2]]@fix

  # Extract the genotype matrices
  sep1 <- objects[[1]]@sep_gt
  imp1 <- objects[[1]]@imp_gt
  genotype_matrix1 <- if (is.null(imp1) || nrow(imp1) == 0) sep1 else imp1
  sep2 <- objects[[2]]@sep_gt
  imp2 <- objects[[2]]@imp_gt
  genotype_matrix2 <- if (is.null(imp2) || nrow(imp2) == 0) sep2 else imp2

  # Replace '.' with NA for missing data and convert to numeric
  genotype_matrix1[genotype_matrix1 == "."] <- NA
  genotype_matrix1 <- as.matrix(apply(genotype_matrix1, c(1, 2), as.numeric))

  genotype_matrix2[genotype_matrix2 == "."] <- NA
  genotype_matrix2 <- as.matrix(apply(genotype_matrix2, c(1, 2), as.numeric))

  # Create a unified set of variants based on chromosome and position
  common_variants <- merge(variant_info1, variant_info2, by = c("CHROM", "POS"), all = TRUE)

  # Number of individuals (assuming diploid, so 2 columns per individual)
  num_individuals1 <- ncol(genotype_matrix1) / 2
  num_individuals2 <- ncol(genotype_matrix2) / 2

  # Initialize a matrix to hold the 2d site frequency spectrum
  sfs_2d <- matrix(0, nrow = num_individuals1 + 1, ncol = num_individuals2 + 1)

  # Process each common variant
  for (i in 1:nrow(common_variants)) {
    chrom <- common_variants[i, "CHROM"]
    pos <- common_variants[i, "POS"]

    # Find the index of this variant in each dataset
    idx1 <- which(variant_info1[,"CHROM"] == chrom & variant_info1[,"POS"] == pos)
    idx2 <- which(variant_info2[,"CHROM"] == chrom & variant_info2[,"POS"] == pos)

    # If the variant is missing in a dataset, we assume it's monomorphic there
    site_data1 <- if (length(idx1) == 1) genotype_matrix1[idx1, ] else rep(0, 2 * num_individuals1)
    site_data2 <- if (length(idx2) == 1) genotype_matrix2[idx2, ] else rep(0, 2 * num_individuals2)

    # Exclude missing data for this site
    valid_data1 <- site_data1[!is.na(site_data1)]
    valid_data2 <- site_data2[!is.na(site_data2)]

    # Count the number of derived alleles (assuming '1' is the derived state)
    derived_count1 <- sum(valid_data1)
    derived_count2 <- sum(valid_data2)

    # Calculate the minor allele count for folded SFS
    if (folded) {
      allele_count1 <- min(derived_count1, length(valid_data1) - derived_count1)
      allele_count2 <- min(derived_count2, length(valid_data2) - derived_count2)
    } else {
      allele_count1 <- derived_count1
      allele_count2 <- derived_count2
    }
    total_alleles_at_site1 <- length(valid_data1)
    total_alleles_at_site2 <- length(valid_data2)
    # Calculate the frequency, adjusting for the varying number of valid alleles
    freq_index1 <- allele_count1 * (num_individuals1) / total_alleles_at_site1
    freq_index2 <- allele_count2 * (num_individuals2) / total_alleles_at_site2

    # Round to the nearest integer to get the discrete frequency category
    freq_category1 <- round(freq_index1) + 1  # R is 1-indexed
    freq_category2 <- round(freq_index2) + 1

    # Update the 2dSFS
    sfs_2d[freq_category1, freq_category2] <- sfs_2d[freq_category1, freq_category2] + 1
  }

  # If the SFS is folded, remove the empty categories.
  if (folded) {
    sfs_2d <- sfs_2d[rowSums(sfs_2d[,-1]) != 0,]
    sfs_2d <- sfs_2d[,colSums(sfs_2d[-1,]) != 0]
  }

  return(sfs_2d)
}
