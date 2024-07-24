#' FixedSites
#'
#' This function counts the number of sites fixed for the alternative allele ("1") in a VCF file.
#' It processes the file in two modes: the entire file at once or in specified windows across the genome.
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach but tailored to process specific genomic windows (`process_vcf_in_windows`).
#'
#' @param vcf_path Path to the VCF file.
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param batch_size The number of variants to be processed in each batch
#'                  (used in batch mode only, default of 10,000 should be suitable for most use cases).
#' @param window_size Size of the window for windowed analysis in base pairs (optional).
#'                   When specified, `skip_size` must also be provided.
#' @param skip_size Number of base pairs to skip between windows (optional).
#'                  Used in conjunction with `window_size` for windowed analysis.
#' @param exclude_ind Optional vector of individual IDs to exclude from the analysis.
#'        If provided, the function will remove these individuals from the genotype matrix
#'        before applying the custom function. Default is NULL, meaning no individuals are excluded.
#'
#' @return
#' In batch mode (no window_size or skip_size provided): A single integer representing the total number of fixed sites for the alternative allele across the entire VCF file.
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'FixedSites', representing the count of fixed sites within each window.
#'
#' @details
#' The function has two modes of operation:
#' 1. Batch Mode: Processes the entire VCF file in batches to count the total number of fixed sites for the alternative allele. Suitable for a general overview of the entire dataset.
#' 2. Window Mode: Processes the VCF file in windows of a specified size and skip distance. This mode is useful for identifying regions with high numbers of fixed sites, which could indicate selective sweeps or regions of low recombination.
#'
#' @examples
#' # Batch mode example
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' num_fixed_sites <- FixedSites(vcf_file)
#'
#' # Window mode example
#' fixed_sites_df <- FixedSites(vcf_file, window_size = 100000, skip_size = 50000)
#'
#' @export


FixedSites <- function(vcf_path, threads = 1, write_log = FALSE, logfile = "log.txt", batch_size = 10000, window_size = NULL, skip_size = NULL, exclude_ind = NULL) {
  if (is.null(window_size) | is.null(skip_size)) {
    # Validate inputs for batch mode
    if (!is.null(window_size) || !is.null(skip_size)) {
      stop("Both 'window_size' and 'skip_size' must be provided for window mode.")
    }

    # Use process_vcf_in_batches to process the file
    batch_results <- process_vcf_in_batches(vcf_path,
                                            batch_size = batch_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            exclude_ind = exclude_ind,
                                            custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL, ploidy = 2) {
                                              allele_freqs <- calculateAlleleFreqs(sep_gt)
                                              # Count fixed sites for the alternative allele in this batch
                                              return(sum(apply(allele_freqs, 1, function(x) any(x[2] == 1))))
                                            })
    # Sum up the counts from all batches
    total_fixed_sites <- sum(do.call("rbind", batch_results))
    return(total_fixed_sites)
  } else {
    # Use process_vcf_in_batches to process the file
    window_results <- process_vcf_in_windows(vcf_path,
                                            window_size = window_size,
                                            skip_size = skip_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            exclude_ind = exclude_ind,
                                            custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals = NULL, pop2_individuals = NULL) {
                                              allele_freqs <- calculateAlleleFreqs(sep_gt)
                                              # Count fixed sites for the alternative allele in this batch
                                              return(c(chrom, start_pos, end_pos, sum(apply(allele_freqs, 1, function(x) any(x[2] == 1)))))
                                            })
    # Bind results per window into a data frame
    fixed_sites_df <- as.data.frame(do.call("rbind", window_results))
    #colnames(fixed_sites_df) <- c("Chromosome", "Start", "End", "FixedSites")
    return(fixed_sites_df)
  }
}


#' SegregatingSites
#'
#' This function counts the number of polymorphic or segregating sites (sites not fixed for the alternative allele)
#' in a VCF file. It processes the file in batches or specified windows across the genome.
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows (`process_vcf_in_windows`).
#'
#' @param vcf_path Path to the VCF file.
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param batch_size The number of variants to be processed in each batch
#'                  (used in batch mode only, default of 10,000 should be suitable for most use cases).
#' @param window_size Size of the window for windowed analysis in base pairs (optional).
#'                   When specified, `skip_size` must also be provided.
#' @param skip_size Number of base pairs to skip between windows (optional).
#'                  Used in conjunction with `window_size` for windowed analysis.
#' @param exclude_ind Optional vector of individual IDs to exclude from the analysis.
#'        If provided, the function will remove these individuals from the genotype matrix
#'        before applying the custom function. Default is NULL, meaning no individuals are excluded.
#'
#' @return
#' In batch mode (no window_size or skip_size provided): A single integer representing the total number of polymorphic sites across the entire VCF file.
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'PolymorphicSites', representing the count of polymorphic sites within each window.
#'
#' @examples
#' # Batch mode example
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' num_polymorphic_sites <- SegregatingSites(vcf_file)
#'
#' # Window mode example
#' polymorphic_sites_df <- SegregatingSites(vcf_file, window_size = 100000, skip_size = 50000)
#'
#' @export

SegregatingSites <- function(vcf_path, threads = 1, write_log = FALSE, logfile = "log.txt", batch_size = 10000, window_size = NULL, skip_size = NULL, exclude_ind = NULL) {
  if (is.null(window_size) || is.null(skip_size)) {
    # Validate inputs for batch mode
    if (!is.null(window_size) || !is.null(skip_size)) {
      stop("Both 'window_size' and 'skip_size' must be provided for window mode.")
    }

    # Batch mode processing
    batch_results <- process_vcf_in_batches(vcf_path,
                                            batch_size = batch_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            exclude_ind = exclude_ind,
                                            custom_function = function(index, fix, sep_gt,pop1_individuals = NULL, pop2_individuals = NULL, ploidy = 2) {
                                              allele_freqs <- calculateAlleleFreqs(sep_gt)
                                              # Count polymorphic sites in this batch
                                              return(sum(apply(allele_freqs, 1, function(x) !any(x == 1) && !all(x == 0))))
                                            })

    # Sum up the counts from all batches
    total_polymorphic_sites <- sum(do.call("rbind", batch_results))
    return(total_polymorphic_sites)
  } else {
    # Window mode processing
    window_results <- process_vcf_in_windows(vcf_path,
                                             window_size = window_size,
                                             skip_size = skip_size,
                                             threads = threads,
                                             write_log = write_log,
                                             logfile = logfile,
                                             exclude_ind = exclude_ind,
                                             custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals = NULL, pop2_individuals = NULL) {
                                               allele_freqs <- calculateAlleleFreqs(sep_gt)
                                               # Count polymorphic sites in this window
                                               return(c(chrom, start_pos, end_pos, sum(apply(allele_freqs, 1, function(x) !any(x == 1) && !all(x == 0)))))
                                             })

    # Bind results per window into a data frame
    polymorphic_sites_df <- as.data.frame(do.call("rbind", window_results))
    colnames(polymorphic_sites_df) <- c("Chromosome", "Start", "End", "PolymorphicSites")
    return(polymorphic_sites_df)
  }
}


#' SingletonSites
#'
#' This function counts the number of singleton sites (sites where a minor allele occurs only once in the sample)
#' in a VCF file. It processes the file in batches or specified windows across the genome.
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows (`process_vcf_in_windows`).
#'
#' @param vcf_path Path to the VCF file.
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param batch_size The number of variants to be processed in each batch
#'                  (used in batch mode only, default of 10,000 should be suitable for most use cases).
#' @param window_size Size of the window for windowed analysis in base pairs (optional).
#'                   When specified, `skip_size` must also be provided.
#' @param skip_size Number of base pairs to skip between windows (optional).
#'                  Used in conjunction with `window_size` for windowed analysis.
#' @param exclude_ind Optional vector of individual IDs to exclude from the analysis.
#'        If provided, the function will remove these individuals from the genotype matrix
#'        before applying the custom function. Default is NULL, meaning no individuals are excluded.
#'
#' @return
#' In batch mode (no window_size or skip_size provided): A single integer representing the total number of singleton sites across the entire VCF file.
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'SingletonSites', representing the count of singleton sites within each window.
#'
#' @examples
#' # Batch mode example
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' num_singleton_sites <- SingletonSites(vcf_file)
#'
#' # Window mode example
#' vcf_path <- "path/to/vcf/file"
#' singleton_sites_df <- SingletonSites(vcf_file, window_size = 100000, skip_size = 50000)
#'
#' @export

SingletonSites <- function(vcf_path, threads = 1, write_log = FALSE, logfile = "log.txt", batch_size = 10000, window_size = NULL, skip_size = NULL, exclude_ind = NULL) {
  if (is.null(window_size) || is.null(skip_size)) {
    # Validate inputs for batch mode
    if (!is.null(window_size) || !is.null(skip_size)) {
      stop("Both 'window_size' and 'skip_size' must be provided for window mode.")
    }

    # Batch mode processing
    batch_results <- process_vcf_in_batches(vcf_path,
                                            batch_size = batch_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            exclude_ind = exclude_ind,
                                            custom_function = function(index, fix, sep_gt,pop1_individuals = NULL, pop2_individuals = NULL, ploidy = 2) {
                                              allele_freqs <- calculateAlleleFreqs(sep_gt)
                                              # Count singleton sites in this batch
                                              return(sum(apply(allele_freqs, 1, function(x) any((x == 1/length(x)) & (names(x) != "0")))))
                                            })

    # Sum up the counts from all batches
    total_singleton_sites <- sum(do.call("rbind", batch_results))
    return(total_singleton_sites)
  } else {
    # Window mode processing
    window_results <- process_vcf_in_windows(vcf_path,
                                             window_size = window_size,
                                             skip_size = skip_size,
                                             threads = threads,
                                             write_log = write_log,
                                             logfile = logfile,
                                             exclude_ind = exclude_ind,
                                             custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals = NULL, pop2_individuals = NULL) {
                                               allele_freqs <- calculateAlleleFreqs(sep_gt)
                                               # Count singleton sites in this window
                                               return(c(chrom, start_pos, end_pos, sum(apply(allele_freqs, 1, function(x) any((x == 1/length(x)) & (names(x) != "0"))))))
                                             })

    # Bind results per window into a data frame
    singleton_sites_df <- as.data.frame(do.call("rbind", window_results))
    colnames(singleton_sites_df) <- c("Chromosome", "Start", "End", "SingletonSites")
    return(singleton_sites_df)
  }
}


#' PrivateAlleles
#'
#' This function calculates the number of private alleles in two populations from a VCF file. (Alleles which are not present in the other popualtion.)
#' It processes the file in batches or specified windows across the genome.
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows (`process_vcf_in_windows`).
#'
#' @param vcf_path Path to the VCF file.
#' @param pop1_individuals Vector of individual names belonging to the first population.
#' @param pop2_individuals Vector of individual names belonging to the second population.
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param batch_size The number of variants to be processed in each batch
#'                  (used in batch mode only, default of 10,000 should be suitable for most use cases).
#' @param window_size Size of the window for windowed analysis in base pairs (optional).
#'                   When specified, `skip_size` must also be provided.
#' @param skip_size Number of base pairs to skip between windows (optional).
#'                  Used in conjunction with `window_size` for windowed analysis.
#'
#' @return
#' In batch mode (no window_size or skip_size provided): A list containing the number of private alleles for each population.
#' In window mode (window_size and skip_size provided): A list of data frames, each with columns 'Chromosome', 'Start', 'End', 'PrivateAllelesPop1', and 'PrivateAllelesPop2', representing the count of private alleles within each window for each population.
#'
#' @examples
#' # Batch mode example
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' pop1_individuals <- c("tsk_0", "tsk_1", "tsk_2")
#' pop2_individuals <- c("tsk_3", "tsk_4", "tsk_5")
#' private_alleles <- PrivateAlleles(vcf_file, pop1_individuals, pop2_individuals)
#'
#' # Window mode example
#' private_alleles_windows <- PrivateAlleles(vcf_file, pop1_individuals, pop2_individuals,
#'                                           window_size = 100000, skip_size = 50000)
#'
#' @export

PrivateAlleles <- function(vcf_path, pop1_individuals, pop2_individuals, threads = 1, write_log = FALSE, logfile = "log.txt", batch_size = 10000, window_size = NULL, skip_size = NULL) {

  if (is.null(window_size) || is.null(skip_size)) {
    # Validate inputs for batch mode
    if (!is.null(window_size) || !is.null(skip_size)) {
      stop("Both 'window_size' and 'skip_size' must be provided for window mode.")
    }
    # Batch mode processing
    batch_results <- process_vcf_in_batches(vcf_path,
                                            batch_size = batch_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            pop1_individuals = pop1_individuals,
                                            pop2_individuals = pop2_individuals,
                                            custom_function = function(index, fix, sep_gt, pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals, ploidy = 2) {
                                              # Separate populations
                                              sep <- separateByPopulations(sep_gt, pop1_names = pop1_individuals, pop2_names = pop2_individuals, rm_ref_alleles = FALSE)
                                              sep_pop1 <- sep$pop1
                                              sep_pop2 <- sep$pop2

                                              # Calculate allele frequencies for each population
                                              allele_freqs_pop1 <- calculateAlleleFreqs(sep_pop1)
                                              allele_freqs_pop2 <- calculateAlleleFreqs(sep_pop2)

                                              private_alleles_pop1 <- 0
                                              private_alleles_pop2 <- 0
                                              # Identify private alleles
                                              for (i in 1:nrow(allele_freqs_pop1)) {
                                                if (allele_freqs_pop1[i,1] == 1 && allele_freqs_pop2[i,2] > 0){
                                                  private_alleles_pop2 <- private_alleles_pop2 + 1
                                                }
                                                if (allele_freqs_pop2[i,1] == 1 && allele_freqs_pop1[i,2] > 0){
                                                  private_alleles_pop1 <- private_alleles_pop1 + 1
                                                }
                                              }
                                              return(c(private_alleles_pop1, private_alleles_pop2))
                                            })

    # Combine the counts from all batches
    private_alleles <- do.call("rbind", batch_results)
    total_private_alleles <- c(sum(private_alleles[,1]), sum(private_alleles[,2]))
    return(total_private_alleles)
  } else {
    # Window mode processing
    window_results <- process_vcf_in_windows(vcf_path,
                                             window_size = window_size,
                                             skip_size = skip_size,
                                             threads = threads,
                                             write_log = write_log,
                                             pop1_individuals = pop1_individuals,
                                             pop2_individuals = pop2_individuals,
                                             logfile = logfile,
                                             custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals) {
                                               # Separate populations
                                               sep <- separateByPopulations(sep_gt, pop1_names = pop1_individuals, pop2_names = pop2_individuals, rm_ref_alleles = FALSE)
                                               sep_pop1 <- sep$pop1
                                               sep_pop2 <- sep$pop2

                                               # Calculate allele frequencies for each population
                                               allele_freqs_pop1 <- calculateAlleleFreqs(sep_pop1)
                                               allele_freqs_pop2 <- calculateAlleleFreqs(sep_pop2)

                                               private_alleles_pop1 <- 0
                                               private_alleles_pop2 <- 0
                                               # Identify private alleles
                                               for (i in 1:nrow(allele_freqs_pop1)) {
                                                 if (allele_freqs_pop1[i,1] == 1 & allele_freqs_pop2[2] > 0){
                                                   private_alleles_pop2 <- private_alleles_pop2 + 1
                                                 }
                                                 if (allele_freqs_pop2[1] == 1 & allele_freqs_pop1[2] > 0){
                                                   private_alleles_pop1 <- private_alleles_pop1 + 1
                                                 }
                                               }
                                               return(c(chrom, start_pos, end_pos, private_alleles_pop1, private_alleles_pop2))
                                             })

    # Bind results per window into a list of data frames
    private_alleles_windows <- as.data.frame(do.call("rbind", window_results))
    colnames(private_alleles_windows) <- c("Chromosome", "Start", "End", "PrivateAllelesPop1", "PrivateAllelesPop2")
    return(private_alleles_windows)
  }
}


#' Heterozygosity Rate
#'
#' This function calculates the rate of heterozygosity for samples in a VCF file. (The proportion of heterozygote genotypes.)
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows (`process_vcf_in_windows`).
#'
#' @param vcf_path Path to the VCF file.
#' @param batch_size The number of variants to be processed in each batch
#'                  (used in batch mode only, default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param window_size Size of the window for windowed analysis in base pairs (optional).
#'                   When specified, `skip_size` must also be provided.
#' @param skip_size Number of base pairs to skip between windows (optional).
#'                  Used in conjunction with `window_size` for windowed analysis.
#' @param exclude_ind Optional vector of individual IDs to exclude from the analysis.
#'        If provided, the function will remove these individuals from the genotype matrix
#'        before applying the custom function. Default is NULL, meaning no individuals are excluded.
#'
#' @return
#' In batch mode (no window_size or skip_size provided): Observed heterozygosity rate averaged over all loci.
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'Ho', representing the observed heterozygosity rate within each window.
#'
#' @examples
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' # Batch mode example
#' Ho <- Heterozygosity(vcf_file)
#' # Window mode example
#' Ho_windows <- Heterozygosity(vcf_file, window_size = 100000, skip_size = 50000)
#'
#' @export

Heterozygosity <- function(vcf_path, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt", window_size = NULL, skip_size = NULL, exclude_ind = NULL) {
  if (is.null(window_size) || is.null(skip_size)) {
    # Validate inputs for batch mode
    if (!is.null(window_size) || !is.null(skip_size)) {
      stop("Both 'window_size' and 'skip_size' must be provided for window mode.")
    }
    # Batch mode processing
    batch_results <- process_vcf_in_batches(vcf_path,
                                            batch_size = batch_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            exclude_ind = exclude_ind,
                                            custom_function = function(index, fix, sep_gt,pop1_individuals = NULL, pop2_individuals = NULL, ploidy = 2) {
                                              # Replace '.' with NA for missing data
                                              sep_gt[sep_gt == "."] <- NA
                                              num_individuals <- ncol(sep_gt) / 2  # assuming diploid organisms
                                              heterozygotes <- 0

                                              # Iterate over individuals by stepping 2 columns at a time
                                              for (indiv_col in seq(1, ncol(sep_gt), by = 2)) {
                                                individual_genotypes <- sep_gt[, c(indiv_col, indiv_col + 1)]
                                                heterozygotes <- heterozygotes + sum(apply(individual_genotypes, 1, function(locus) {
                                                  # Only count as a heterozygote if neither allele is NA and they are different
                                                  !is.na(locus[1]) && !is.na(locus[2]) && locus[1] != locus[2]
                                                }))
                                              }

                                              # Return the count of heterozygotes and the number of loci
                                              return(c(heterozygotes, num_individuals * nrow(sep_gt)))
                                            })

    # Combine the counts from all batches and calculate the mean
    all_heterozygotes <- sum(sapply(batch_results, function(x) x[1]))
    all_valid_loci <- sum(sapply(batch_results, function(x) x[2]))
    Ho <- all_heterozygotes / all_valid_loci

    return(Ho)
  } else {
    # Window mode processing
    window_results <- process_vcf_in_windows(vcf_path,
                                             window_size = window_size,
                                             skip_size = skip_size,
                                             threads = threads,
                                             write_log = write_log,
                                             logfile = logfile,
                                             exclude_ind = exclude_ind,
                                             custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals = NULL, pop2_individuals = NULL) {
                                               sep_gt[sep_gt == "."] <- NA
                                               num_individuals <- ncol(sep_gt) / 2  # assuming diploid organisms
                                               heterozygotes <- 0

                                               for (indiv_col in seq(1, ncol(sep_gt), by = 2)) {
                                                 individual_genotypes <- sep_gt[, c(indiv_col, indiv_col + 1)]
                                                 heterozygotes <- heterozygotes + sum(apply(individual_genotypes, 1, function(locus) {
                                                   !is.na(locus[1]) && !is.na(locus[2]) && locus[1] != locus[2]
                                                 }))
                                               }

                                               total_loci <- num_individuals * nrow(sep_gt)
                                               Ho <- if (total_loci > 0) heterozygotes / total_loci else NA
                                               return(c(chrom, start_pos, end_pos, Ho))
                                             })

    # Bind results per window into a data frame
    Ho_windows <- as.data.frame(do.call("rbind", window_results))
    colnames(Ho_windows) <- c("Chromosome", "Start", "End", "Ho")
    return(Ho_windows)
  }
}


#' Pi
#'
#' This function calculates the nucleotide diversity (Pi) for a sample in a VCF file as defined by Nei & Li, 1979 (https://doi.org/10.1073/pnas.76.10.5269).
#' The formula used for this is equivalent to the one used in vcftools --window-pi (https://vcftools.sourceforge.net/man_latest.html).
#' Handling missing alleles at one site is equivalent to Korunes & Samuk, 2021 ( https://doi.org/10.1111/1755-0998.13326).
#' The function calculates the number of monomorphic sites using the sequence length and the number of variants in the VCF file. This assumes, that all sites not present in the VCF file are invariant sites, which will underestimate Pi, because of commonly done (and necessary) variant filtering. However, otherwise this calculation would only work with VCF files that include all monomorphic sites, which is quite unpractical for common use cases and will increase computational demands significantly.
#' If you happen to know the number of filtered our sites vs the number of monomorphic sites, please use the number of monomorphic + the number of polymorphic (number of variants in your VCF) sites as the sequence length to get the most accurate estimation of Pi. (This does not work for the window mode of this function, which assumes the sequence length to be the window size.)
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows (`process_vcf_in_windows`).
#'
#' @param vcf_path Path to the VCF file.
#' @param seq_length Total length of the sequence in number of bases (used in batch mode only).
#' @param batch_size The number of variants to be processed in each batch
#'                  (used in batch mode only, default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param window_size Size of the window for windowed analysis in base pairs (optional).
#'                   When specified, `skip_size` must also be provided.
#' @param skip_size Number of base pairs to skip between windows (optional).
#'                  Used in conjunction with `window_size` for windowed analysis.
#' @param exclude_ind Optional vector of individual IDs to exclude from the analysis.
#'        If provided, the function will remove these individuals from the genotype matrix
#'        before applying the custom function. Default is NULL, meaning no individuals are excluded.
#'
#' @return
#' In batch mode (no window_size or skip_size provided): Nucleotide diversity (Pi) across the sequence.
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'Pi', representing the nucleotide diversity within each window.
#'
#' @examples
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' total_sequence_length <- 999299  # Total length of the sequence in vcf
#' # Batch mode example
#' pi_value <- Pi(vcf_file, total_sequence_length)
#' # Window mode example
#' pi_windows <- Pi(vcf_file, seq_length = total_sequence_length,
#'                  window_size = 100000, skip_size = 50000)
#'
#' @export

Pi <- function(vcf_path, seq_length, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt", window_size = NULL, skip_size = NULL, exclude_ind = NULL) {
  if (is.null(window_size) || is.null(skip_size)) {
    # Validate inputs for batch mode
    if (!is.null(window_size) || !is.null(skip_size)) {
      stop("Both 'window_size' and 'skip_size' must be provided for window mode.")
    }
    # Batch mode processing
    batch_results <- process_vcf_in_batches(vcf_path,
                                            batch_size = batch_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            exclude_ind = exclude_ind,
                                            custom_function = function(index, fix, sep_gt,pop1_individuals = NULL, pop2_individuals = NULL, ploidy = 2) {
                                              sep_gt[sep_gt == "."] <- NA  # Replace '.' with NA for missing data
                                              N_mismatches_batch <- as.numeric(0)
                                              N_comparisons_batch <- as.numeric(0)
                                              num_chroms <- ncol(sep_gt)

                                              for (site_index in seq_len(nrow(sep_gt))) {
                                                site_genotypes <- sep_gt[site_index, !is.na(sep_gt[site_index, ])]
                                                site_allele_freqs <- table(site_genotypes)  # Allele frequencies at this site

                                                # Total alleles at this site (not missing)
                                                N_non_missing_chr <- sum(site_allele_freqs)

                                                # Number of actual nucleotide differences (mismatches) for the site
                                                N_site_mismatches <- 0
                                                for (allele_count in site_allele_freqs) {
                                                  N_site_mismatches <- N_site_mismatches + (allele_count * (N_non_missing_chr - allele_count))
                                                }

                                                N_mismatches_batch <- N_mismatches_batch + N_site_mismatches
                                                N_comparisons_batch <- N_comparisons_batch + (N_non_missing_chr * (N_non_missing_chr - 1))
                                              }

                                              return(c(N_mismatches_batch, N_comparisons_batch, nrow(sep_gt), num_chroms))
                                            })

    # Combine the counts from all batches
    all_N_mismatches <- sum(sapply(batch_results, function(x) x[1]))
    all_N_comparisons <- sum(sapply(batch_results, function(x) x[2]))
    all_variants_counted <- sum(sapply(batch_results, function(x) x[3]))
    num_chroms <- batch_results[[1]][4]
    # Including monomorphic sites
    N_monomorphic_sites <- seq_length - all_variants_counted
    N_pairs_total <- as.numeric(all_N_comparisons) + (as.numeric(N_monomorphic_sites) * as.numeric(num_chroms) * (as.numeric(num_chroms) - 1))

    # Pi calculation for the sequence
    pi_value <- all_N_mismatches / N_pairs_total

    return(pi_value)
  } else {
    # Window mode processing
    window_results <- process_vcf_in_windows(vcf_path,
                                             window_size = window_size,
                                             skip_size = skip_size,
                                             threads = threads,
                                             write_log = write_log,
                                             logfile = logfile,
                                             exclude_ind = exclude_ind,
                                             custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals = NULL, pop2_individuals = NULL) {
                                               sep_gt[sep_gt == "."] <- NA  # Replace '.' with NA for missing data
                                               N_mismatches_batch <- as.numeric(0)
                                               N_comparisons_batch <- as.numeric(0)
                                               num_chroms <- ncol(sep_gt)

                                               for (site_index in seq_len(nrow(sep_gt))) {
                                                 site_genotypes <- sep_gt[site_index, !is.na(sep_gt[site_index, ])]
                                                 site_allele_freqs <- table(site_genotypes)  # Allele frequencies at this site

                                                 # Total alleles at this site (not missing)
                                                 N_non_missing_chr <- sum(site_allele_freqs)

                                                 # Number of actual nucleotide differences (mismatches) for the site
                                                 N_site_mismatches <- 0
                                                 for (allele_count in site_allele_freqs) {
                                                   N_site_mismatches <- N_site_mismatches + (allele_count * (N_non_missing_chr - allele_count))
                                                 }

                                                 N_mismatches_batch <- N_mismatches_batch + N_site_mismatches
                                                 N_comparisons_batch <- N_comparisons_batch + (N_non_missing_chr * (N_non_missing_chr - 1))
                                               }

                                               window_length <- end_pos - start_pos
                                               N_monomorphic_sites <- window_length - nrow(sep_gt)
                                               N_pairs_batch <- N_comparisons_batch + (N_monomorphic_sites * num_chroms * (num_chroms - 1))
                                               pi_window <- N_mismatches_batch / N_pairs_batch

                                               return(c(chrom, start_pos, end_pos, pi_window))
                                             })

    # Bind results per window into a data frame
    pi_windows <- as.data.frame(do.call("rbind", window_results))
    colnames(pi_windows) <- c("Chromosome", "Start", "End", "Pi")
    return(pi_windows)
  }
}


#' TajimasD
#'
#' This function calculates Tajima's D statistic for a given dataset (Tajima, 1989 (10.1093/genetics/123.3.585)).
#' The formula used for this is equivalent to the one used in vcftools --TajimaD (https://vcftools.sourceforge.net/man_latest.html).
#' The function calculates the number of monomorphic sites using the sequence length and the number of variants in the VCF file. This assumes, that all sites not present in the VCF file are invariant sites, which will underestimate Pi, because of commonly done (and necessary) variant filtering. However, otherwise this calculation would only work with VCF files that include all monomorphic sites, which is quite unpractical for common use cases and will increase computational demands significantly.
#' If you happen to know the number of filtered our sites vs the number of monomorphic sites, please use the number of monomorphic + the number of polymorphic (number of variants in your VCF) sites as the sequence length to get the most accurate estimation of Pi. (This does not work for the window mode of this function, which assumes the sequence length to be the window size.)
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows (`process_vcf_in_windows`).
#'
#' @param vcf_path Path to the VCF file.
#' @param seq_length Total length of the sequence in number of bases (used in batch mode only).
#' @param batch_size The number of variants to be processed in each batch
#'                  (used in batch mode only, default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param window_size Size of the window for windowed analysis in base pairs (optional).
#'                   When specified, `skip_size` must also be provided.
#' @param skip_size Number of base pairs to skip between windows (optional).
#'                  Used in conjunction with `window_size` for windowed analysis.
#' @param exclude_ind Optional vector of individual IDs to exclude from the analysis.
#'        If provided, the function will remove these individuals from the genotype matrix
#'        before applying the custom function. Default is NULL, meaning no individuals are excluded.
#'
#' @return
#' In batch mode (no window_size or skip_size provided): Tajima's D value.
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'TajimasD', representing Tajima's D within each window.
#'
#' @examples
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' total_sequence_length <- 999299  # Total length of the sequence
#' # Batch mode example
#' tajimas_d <- TajimasD(vcf_file, total_sequence_length)
#' # Window mode example
#' tajimas_d_windows <- TajimasD(vcf_file, seq_length = total_sequence_length,
#'                               window_size = 100000, skip_size = 50000)
#'
#' @importFrom stats na.omit
#' @export

TajimasD <- function(vcf_path, seq_length, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt", window_size = NULL, skip_size = NULL, exclude_ind = NULL) {
  if (is.null(window_size) || is.null(skip_size)) {
    # Validate inputs for batch mode
    if (!is.null(window_size) || !is.null(skip_size)) {
      stop("Both 'window_size' and 'skip_size' must be provided for window mode.")
    }
    # Batch mode processing
    batch_results <- process_vcf_in_batches(vcf_path,
                                            batch_size = batch_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            exclude_ind = exclude_ind,
                                            custom_function = function(index, fix, sep_gt,pop1_individuals = NULL, pop2_individuals = NULL, ploidy = 2) {
                                              sep_gt[sep_gt == "."] <- NA  # Replace '.' with NA for missing data
                                              num_chroms <- ncol(sep_gt)
                                              n <- num_chroms  # Number of chromosomes

                                              # Calculate segregating sites (S) for the batch
                                              S_batch <- sum(apply(sep_gt, 1, function(row) {
                                                length(unique(na.omit(row))) > 1
                                              }))

                                              # Calculate Pi for the batch
                                              N_mismatches_batch <- as.numeric(0)
                                              N_comparisons_batch <- as.numeric(0)
                                              for (site_index in seq_len(nrow(sep_gt))) {
                                                site_genotypes <- sep_gt[site_index, !is.na(sep_gt[site_index, ])]
                                                site_allele_freqs <- table(site_genotypes)  # Allele frequencies at this site

                                                # Total alleles at this site (not missing)
                                                N_non_missing_chr <- sum(site_allele_freqs)

                                                # Number of actual nucleotide differences (mismatches) for the site
                                                N_site_mismatches <- 0
                                                for (allele_count in site_allele_freqs) {
                                                  N_site_mismatches <- N_site_mismatches + (allele_count * (N_non_missing_chr - allele_count))
                                                }

                                                N_mismatches_batch <- N_mismatches_batch + N_site_mismatches
                                                N_comparisons_batch <- N_comparisons_batch + (N_non_missing_chr * (N_non_missing_chr - 1))
                                              }

                                              return(c(S_batch, N_mismatches_batch, N_comparisons_batch, nrow(sep_gt), num_chroms))
                                            })

    # Aggregate the results from all batches
    total_S <- sum(sapply(batch_results, function(x) x[1]))
    all_N_mismatches <- sum(sapply(batch_results, function(x) x[2]))
    all_N_comparisons <- sum(sapply(batch_results, function(x) x[3]))
    all_variants_counted <- sum(sapply(batch_results, function(x) x[4]))
    num_chroms <- batch_results[[1]][5]  # Assuming number of chromosomes is consistent across all batches

    # Including monomorphic sites
    N_monomorphic_sites <- seq_length - all_variants_counted
    N_pairs_total <- as.numeric(all_N_comparisons) + (as.numeric(N_monomorphic_sites) * as.numeric(num_chroms) * (as.numeric(num_chroms) - 1))

    # Pi calculation for the sequence
    total_pi <- all_N_mismatches / N_pairs_total

    # Constants for Tajima's D calculation
    n <- num_chroms
    i_values <- 1:(n-1)
    a1 <- sum(1 / i_values)
    a2 <- sum(1 / (i_values^2))
    b1 <- (n + 1) / (3 * (n - 1))
    b2 <- (2 * (n^2 + n + 3)) / (9 * n * (n - 1))
    c1 <- b1 - (1 / a1)
    c2 <- b2 - ((n + 2) / (a1 * n)) + (a2 / (a1^2))
    e1 <- c1 / a1
    e2 <- c2 / (a1^2 + a2)

    pi <- total_pi * seq_length  # Adjust total pi by sequence length to get average pi per site

    # Tajima's D calculation
    denominator <- sqrt((e1 * total_S) + (e2 * total_S * (total_S - 1)))
    D <- (pi - (total_S / a1)) / denominator

    return(D)
  } else {
    # Window mode processing
    window_results <- process_vcf_in_windows(vcf_path,
                                             window_size = window_size,
                                             skip_size = skip_size,
                                             threads = threads,
                                             write_log = write_log,
                                             logfile = logfile,
                                             exclude_ind = exclude_ind,
                                             custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals = NULL, pop2_individuals = NULL) {
                                               sep_gt[sep_gt == "."] <- NA  # Replace '.' with NA for missing data
                                               num_chroms <- ncol(sep_gt)
                                               n <- num_chroms  # Number of chromosomes

                                               # Calculate segregating sites (S) for the batch
                                               S_window <- sum(apply(sep_gt, 1, function(row) {
                                                 length(unique(na.omit(row))) > 1
                                               }))

                                               # Calculate Pi for the batch
                                               N_mismatches_batch <- as.numeric(0)
                                               N_comparisons_batch <- as.numeric(0)
                                               for (site_index in seq_len(nrow(sep_gt))) {
                                                 site_genotypes <- sep_gt[site_index, !is.na(sep_gt[site_index, ])]
                                                 site_allele_freqs <- table(site_genotypes)  # Allele frequencies at this site

                                                 # Total alleles at this site (not missing)
                                                 N_non_missing_chr <- sum(site_allele_freqs)

                                                 # Number of actual nucleotide differences (mismatches) for the site
                                                 N_site_mismatches <- 0
                                                 for (allele_count in site_allele_freqs) {
                                                   N_site_mismatches <- N_site_mismatches + (allele_count * (N_non_missing_chr - allele_count))
                                                 }

                                                 N_mismatches_batch <- N_mismatches_batch + N_site_mismatches
                                                 N_comparisons_batch <- N_comparisons_batch + (N_non_missing_chr * (N_non_missing_chr - 1))
                                               }

                                               window_length <- end_pos - start_pos
                                               N_monomorphic_sites <- window_length - nrow(sep_gt)
                                               pi_window <- (as.numeric(N_mismatches_batch) / (as.numeric(N_comparisons_batch) + (as.numeric(N_monomorphic_sites) * as.numeric(num_chroms) * ((as.numeric(num_chroms) - 1))))) * window_length


                                               # Constants for Tajima's D calculation
                                               n <- num_chroms
                                               i_values <- 1:(n-1)
                                               a1 <- sum(1 / i_values)
                                               a2 <- sum(1 / (i_values^2))
                                               b1 <- (n + 1) / (3 * (n - 1))
                                               b2 <- (2 * (n^2 + n + 3)) / (9 * n * (n - 1))
                                               c1 <- b1 - (1 / a1)
                                               c2 <- b2 - ((n + 2) / (a1 * n)) + (a2 / (a1^2))
                                               e1 <- c1 / a1
                                               e2 <- c2 / (a1^2 + a2)

                                               # Tajima's D calculation for the window
                                               denominator <- sqrt((e1 * S_window) + (e2 * S_window * (S_window - 1)))
                                               D_window <- (pi_window - (S_window / a1)) / denominator

                                               return(c(chrom, start_pos, end_pos, D_window))
                                             })

    # Bind results per window into a data frame
    tajimas_d_windows <- as.data.frame(do.call("rbind", window_results))
    colnames(tajimas_d_windows) <- c("Chromosome", "Start", "End", "TajimasD")
    return(tajimas_d_windows)
  }
}

#' WattersonsTheta
#'
#' This function calculates Watterson's Theta, a measure for neutrality, from a VCF file (Watterson, 1975 (https://doi.org/10.1016/0040-5809(75)90020-9)).
#' The function calculates the number of monomorphic sites using the sequence length and the number of variants in the VCF file. This assumes, that all sites not present in the VCF file are invariant sites, which will underestimate Pi, because of commonly done (and necessary) variant filtering. However, otherwise this calculation would only work with VCF files that include all monomorphic sites, which is quite unpractical for common use cases and will increase computational demands significantly.
#' If you happen to know the number of filtered our sites vs the number of monomorphic sites, please use the number of monomorphic + the number of polymorphic (number of variants in your VCF) sites as the sequence length to get the most accurate estimation of Pi. (This does not work for the window mode of this function, which assumes the sequence length to be the window size.)
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows (`process_vcf_in_windows`).
#'
#' @param vcf_path Path to the VCF file.
#' @param seq_length The length of the sequence in the data set (used in batch mode only).
#' @param batch_size The number of variants to be processed in each batch
#'                  (used in batch mode only, default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param window_size Size of the window for windowed analysis in base pairs (optional).
#'                   When specified, `skip_size` must also be provided.
#' @param skip_size Number of base pairs to skip between windows (optional).
#'                  Used in conjunction with `window_size` for windowed analysis.
#' @param exclude_ind Optional vector of individual IDs to exclude from the analysis.
#'        If provided, the function will remove these individuals from the genotype matrix
#'        before applying the custom function. Default is NULL, meaning no individuals are excluded.
#'
#' @return
#' In batch mode (no window_size or skip_size provided): Watterson's theta value normalized by the sequence length.
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'WattersonsTheta', representing Watterson's theta within each window normalized by the window length.
#'
#' @examples
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' total_sequence_length <- 999299  # Total length of the sequence
#' # Batch mode example
#' wattersons_theta <- WattersonsTheta(vcf_file, total_sequence_length)
#' # Window mode example
#' wattersons_theta_windows <- WattersonsTheta(vcf_file, seq_length = total_sequence_length,
#'                                             window_size = 100000, skip_size = 50000)
#'
#' @importFrom stats na.omit
#' @export

WattersonsTheta <- function(vcf_path, seq_length, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt", window_size = NULL, skip_size = NULL, exclude_ind = NULL) {
  if (is.null(window_size) || is.null(skip_size)) {
    # Validate inputs for batch mode
    if (!is.null(window_size) || !is.null(skip_size)) {
      stop("Both 'window_size' and 'skip_size' must be provided for window mode.")
    }
    # Batch mode processing
    batch_results <- process_vcf_in_batches(vcf_path,
                                            batch_size = batch_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            exclude_ind = exclude_ind,
                                            custom_function = function(index, fix, sep_gt,pop1_individuals = NULL, pop2_individuals = NULL, ploidy = 2) {
                                              sep_gt[sep_gt == "."] <- NA  # Replace '.' with NA for missing data

                                              # Calculate segregating sites (S) for the batch
                                              S_batch <- sum(apply(sep_gt, 1, function(row) {
                                                length(unique(na.omit(row))) > 1
                                              }))

                                              return(c(S_batch, ncol(sep_gt)))  # Return S and number of chromosomes
                                            })

    # Aggregate the results from all batches
    total_S <- sum(sapply(batch_results, function(x) x[1]))
    num_chroms <- batch_results[[1]][2]  # Assuming number of chromosomes is consistent across all batches

    # Constants for Watterson's Theta calculation
    n <- num_chroms  # Sample size
    i_values <- 1:(n-1)
    a1 <- sum(1 / i_values)

    # Watterson's Theta calculation
    wattersons_theta <- total_S / a1

    # Normalize by the sequence length
    normalized_wattersons_theta <- wattersons_theta / seq_length

    return(normalized_wattersons_theta)
  } else {
    # Window mode processing
    window_results <- process_vcf_in_windows(vcf_path,
                                             window_size = window_size,
                                             skip_size = skip_size,
                                             threads = threads,
                                             write_log = write_log,
                                             logfile = logfile,
                                             exclude_ind = exclude_ind,
                                             custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals = NULL, pop2_individuals = NULL) {
                                               sep_gt[sep_gt == "."] <- NA  # Replace '.' with NA for missing data

                                               # Calculate segregating sites (S) for the window
                                               S_window <- sum(apply(sep_gt, 1, function(row) {
                                                 length(unique(na.omit(row))) > 1
                                               }))

                                               num_chroms <- ncol(sep_gt)
                                               n <- num_chroms  # Sample size for the window
                                               i_values <- 1:(n-1)
                                               a1 <- sum(1 / i_values)

                                               # Watterson's Theta calculation for the window
                                               wattersons_theta_window <- S_window / a1

                                               # Normalize by the window length
                                               window_length <- end_pos - start_pos + 1
                                               normalized_wattersons_theta_window <- wattersons_theta_window / window_length

                                               return(c(chrom, start_pos, end_pos, normalized_wattersons_theta_window))
                                             })

    # Bind results per window into a data frame
    wattersons_theta_windows <- as.data.frame(do.call("rbind", window_results))
    colnames(wattersons_theta_windows) <- c("Chromosome", "Start", "End", "WattersonsTheta")
    return(wattersons_theta_windows)
  }
}


#' Dxy
#'
#' This function calculates the average number of nucleotide differences per site (Dxy) between two populations from a VCF file (Nei & Li, 1979 (https://doi.org/10.1073/pnas.76.10.5269)).
#' Handling missing alleles at one site is equivalent to Korunes & Samuk, 2021 ( https://doi.org/10.1111/1755-0998.13326).
#' The function calculates the number of monomorphic sites using the sequence length and the number of variants in the VCF file. This assumes, that all sites not present in the VCF file are invariant sites, which will underestimate Pi, because of commonly done (and necessary) variant filtering. However, otherwise this calculation would only work with VCF files that include all monomorphic sites, which is quite unpractical for common use cases and will increase computational demands significantly.
#' If you happen to know the number of filtered our sites vs the number of monomorphic sites, please use the number of monomorphic + the number of polymorphic (number of variants in your VCF) sites as the sequence length to get the most accurate estimation of Pi. (This does not work for the window mode of this function, which assumes the sequence length to be the window size.)
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows (`process_vcf_in_windows`).
#'
#' @param vcf_path Path to the VCF file.
#' @param pop1_individuals Vector of individual names belonging to the first population.
#' @param pop2_individuals Vector of individual names belonging to the second population.
#' @param seq_length Length of the sequence in number of bases, including monomorphic sites (used in batch mode only).
#' @param batch_size The number of variants to be processed in each batch
#'                  (used in batch mode only, default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param window_size Size of the window for windowed analysis in base pairs (optional).
#'                   When specified, `skip_size` must also be provided.
#' @param skip_size Number of base pairs to skip between windows (optional).
#'                  Used in conjunction with `window_size` for windowed analysis.
#'
#' @return
#' In batch mode (no window_size or skip_size provided): The average number of nucleotide substitutions per site between the individuals of two populations (Dxy).
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'Dxy', representing the average nucleotide differences within each window.
#'
#' @examples
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' pop1_individuals <- c("tsk_0", "tsk_1", "tsk_2")
#' pop2_individuals <- c("tsk_3", "tsk_4", "tsk_5")
#' total_sequence_length <- 999299  # Total length of the sequence
#' # Batch mode example
#' dxy_value <- Dxy(vcf_file, pop1_individuals, pop2_individuals, total_sequence_length)
#' # Window mode example
#' dxy_windows <- Dxy(vcf_file, pop1_individuals, pop2_individuals, seq_length = total_sequence_length,
#'                    window_size = 100000, skip_size = 50000)
#'
#' @export

Dxy <- function(vcf_path, pop1_individuals, pop2_individuals, seq_length, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt", window_size = NULL, skip_size = NULL) {
  if (is.null(window_size) || is.null(skip_size)) {
    # Validate inputs for batch mode
    if (!is.null(window_size) || !is.null(skip_size)) {
      stop("Both 'window_size' and 'skip_size' must be provided for window mode.")
    }
    # Batch mode processing
    batch_results <- process_vcf_in_batches(vcf_path,
                                            batch_size = batch_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            pop1_individuals = pop1_individuals,
                                            pop2_individuals = pop2_individuals,
                                            custom_function = function(index, fix, sep_gt,pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals, ploidy = 2) {
                                              # Separate populations
                                              sep <- separateByPopulations(sep_gt, pop1_names = pop1_individuals, pop2_names = pop2_individuals, rm_ref_alleles = FALSE)
                                              pop1_genotypes <- sep$pop1
                                              pop2_genotypes <- sep$pop2

                                              # Initialize counters for nucleotide differences and comparisons
                                              diffs_batch <- as.numeric(0)
                                              comparisons_batch <- as.numeric(0)

                                              # Iterate over polymorphic sites
                                              for (site_index in seq_len(nrow(pop1_genotypes))) {
                                                site_genotypes1 <- pop1_genotypes[site_index, ]
                                                site_genotypes2 <- pop2_genotypes[site_index, ]

                                                # Calculate differences for each allele combination between populations
                                                for (i in seq_along(site_genotypes1)) {
                                                  for (j in seq_along(site_genotypes2)) {
                                                    if (site_genotypes1[i] != "." && site_genotypes2[j] != ".") {
                                                      diffs_batch <- diffs_batch + as.numeric(site_genotypes1[i] != site_genotypes2[j])
                                                      comparisons_batch <- comparisons_batch + 1
                                                    }
                                                  }
                                                }
                                              }

                                              return(c(diffs_batch, comparisons_batch, nrow(sep_gt)))
                                            })

    # Aggregate the results from all batches
    total_diffs <- sum(sapply(batch_results, function(x) x[1]))
    total_comparisons <- sum(sapply(batch_results, function(x) x[2]))
    all_variants_counted <- sum(sapply(batch_results, function(x) x[3]))

    # Including monomorphic sites in total comparisons
    monomorphic_sites <- seq_length - all_variants_counted
    total_comparisons <- total_comparisons + (monomorphic_sites * (length(pop1_individuals) * 2) * (length(pop2_individuals) * 2))

    # Calculate Dxy
    dxy_value <- total_diffs / total_comparisons

    return(dxy_value)

  } else {
    # Window mode processing
    window_results <- process_vcf_in_windows(vcf_path,
                                             window_size = window_size,
                                             skip_size = skip_size,
                                             threads = threads,
                                             write_log = write_log,
                                             logfile = logfile,
                                             pop1_individuals = pop1_individuals,
                                             pop2_individuals = pop2_individuals,
                                             custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals, pop2_individuals) {
                                               sep <- separateByPopulations(sep_gt, pop1_names = pop1_individuals, pop2_names = pop2_individuals, rm_ref_alleles = FALSE)
                                               pop1_genotypes <- sep$pop1
                                               pop2_genotypes <- sep$pop2

                                               diffs_window <- as.numeric(0)
                                               comparisons_window <- as.numeric(0)

                                               for (site_index in seq_len(nrow(pop1_genotypes))) {
                                                 site_genotypes1 <- pop1_genotypes[site_index, ]
                                                 site_genotypes2 <- pop2_genotypes[site_index, ]

                                                 for (i in seq_along(site_genotypes1)) {
                                                   for (j in seq_along(site_genotypes2)) {
                                                     if (site_genotypes1[i] != "." && site_genotypes2[j] != ".") {
                                                       diffs_window <- diffs_window + as.numeric(site_genotypes1[i] != site_genotypes2[j])
                                                       comparisons_window <- comparisons_window + 1
                                                     }
                                                   }
                                                 }
                                               }

                                               window_length <- end_pos - start_pos
                                               monomorphic_sites <- window_length - nrow(sep_gt)
                                               comparisons_window <- comparisons_window + (monomorphic_sites * (length(pop1_individuals) * 2) * (length(pop2_individuals) * 2))
                                               dxy_window <- diffs_window / comparisons_window

                                               return(c(chrom, start_pos, end_pos, dxy_window))
                                             })

    # Bind results per window into a data frame
    dxy_windows <- as.data.frame(do.call("rbind", window_results))
    colnames(dxy_windows) <- c("Chromosome", "Start", "End", "Dxy")
    return(dxy_windows)
  }
}


#' Fst
#'
#' This function calculates the fixation index (Fst) between two populations from a VCF file using the method of Weir and Cockerham (1984).
#' The formula used for this is equivalent to the one used in vcftools --weir-fst-pop (https://vcftools.sourceforge.net/man_latest.html).
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows (`process_vcf_in_windows`).
#'
#' @param vcf_path Path to the VCF file.
#' @param pop1_individuals Vector of individual names belonging to the first population.
#' @param pop2_individuals Vector of individual names belonging to the second population.
#' @param weighted Logical, whether weighted Fst or mean Fst is returned (Default = FALSE (mean Fst is returned)).
#' @param batch_size The number of variants to be processed in each batch
#'                  (used in batch mode only, default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param window_size Size of the window for windowed analysis in base pairs (optional).
#'                   When specified, `skip_size` must also be provided.
#' @param skip_size Number of base pairs to skip between windows (optional).
#'                  Used in conjunction with `window_size` for windowed analysis.
#'
#' @return
#' In batch mode (no window_size or skip_size provided): Fst value (either mean or weighted).
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'Fst', representing the fixation index within each window.
#'
#' @examples
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' pop1_individuals <- c("tsk_0", "tsk_1", "tsk_2")
#' pop2_individuals <- c("tsk_3", "tsk_4", "tsk_5")
#' # Batch mode example
#' fst_value <- Fst(vcf_file, pop1_individuals, pop2_individuals, weighted = TRUE)
#' # Window mode example
#' fst_windows <- Fst(vcf_file, pop1_individuals, pop2_individuals, weighted = TRUE,
#'                    window_size = 100000, skip_size = 50000)
#'
#' @export

Fst <- function(vcf_path, pop1_individuals, pop2_individuals, weighted = FALSE, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt", window_size = NULL, skip_size = NULL) {
  if (is.null(window_size) || is.null(skip_size)) {
    # Validate inputs for batch mode
    if (!is.null(window_size) || !is.null(skip_size)) {
      stop("Both 'window_size' and 'skip_size' must be provided for window mode.")
    }
    # Batch mode processing
    batch_results <- process_vcf_in_batches(vcf_path,
                                            batch_size = batch_size,
                                            threads = threads,
                                            write_log = write_log,
                                            logfile = logfile,
                                            pop1_individuals = pop1_individuals,
                                            pop2_individuals = pop2_individuals,
                                            custom_function = function(index, fix, sep_gt,pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals, ploidy = 2) {
                                              # Separate populations
                                              sep <- separateByPopulations(sep_gt, pop1_names = pop1_individuals, pop2_names = pop2_individuals, rm_ref_alleles = FALSE)
                                              genotype_matrix1 <- sep$pop1
                                              genotype_matrix2 <- sep$pop2
                                              allele_freqs1 <- calculateAlleleFreqs(genotype_matrix1)
                                              allele_freqs2 <- calculateAlleleFreqs(genotype_matrix2)

                                              # Preparing variables for
                                              sum1 <- 0
                                              sum2 <- 0
                                              sum3 <- 0
                                              count <- 0

                                              # Iterate over polymorphic sites
                                              for (i in seq_len(nrow(genotype_matrix1))) {
                                                ## Sample size values ##
                                                # Get the number of non-missing individuals
                                                n_p1 <- length(genotype_matrix1[i,genotype_matrix1[i,] != "."]) / 2
                                                n_p2 <- length(genotype_matrix2[i,genotype_matrix2[i,] != "."]) / 2

                                                # Total sum of individuals across populations
                                                n_sum <- n_p1 + n_p2

                                                # Average number of individuals per population
                                                nbar <- n_sum / 2

                                                # Sum of squares of the number of individuals
                                                sum_nsqr <- n_p1^2 + n_p2^2

                                                # Effective Sample Size per Population
                                                nc <- (n_sum - (sum_nsqr / n_sum))

                                                ## Allele frequency values ##
                                                # Variant allele frequencies
                                                p1_a0 <- allele_freqs1[i,1]
                                                p1_a1 <- allele_freqs1[i,2]
                                                p2_a0 <- allele_freqs2[i,1]
                                                p2_a1 <- allele_freqs2[i,2]

                                                # Total allele counts
                                                pc_a0 <- (p1_a0 * (n_p1*2)) + (p2_a0 * (n_p2*2)) # times 2 cause diploidy
                                                pc_a1 <- (p1_a1 * (n_p1*2)) + (p2_a1 * (n_p2*2))

                                                # Average allele frequencies
                                                pbar_a0 <- pc_a0 / (n_sum * 2) # times 2 cause diploidy
                                                pbar_a1 <- pc_a1 / (n_sum * 2)

                                                # Heterozygosity values
                                                # Counts of heterozygotes
                                                N_het_p1 <- 0
                                                for (indiv_col in seq(1, ncol(genotype_matrix1), by = 2)) {
                                                  individual_genotypes <- genotype_matrix1[i, c(indiv_col, indiv_col + 1)]
                                                  is_heterozygote <- !is.na(individual_genotypes[1]) && !is.na(individual_genotypes[2]) && individual_genotypes[1] != individual_genotypes[2]
                                                  if (is_heterozygote) {
                                                    N_het_p1 <- N_het_p1 + 1
                                                  }
                                                }
                                                N_het_p2 <- 0
                                                for (indiv_col in seq(1, ncol(genotype_matrix2), by = 2)) {
                                                  individual_genotypes <- genotype_matrix2[i, c(indiv_col, indiv_col + 1)]
                                                  is_heterozygote <- !is.na(individual_genotypes[1]) && !is.na(individual_genotypes[2]) && individual_genotypes[1] != individual_genotypes[2]
                                                  if (is_heterozygote) {
                                                    N_het_p2 <- N_het_p2 + 1
                                                  }
                                                }

                                                # Average heterozygosity
                                                hbar <- (N_het_p1 + N_het_p2) / n_sum

                                                # Squared Deviations of Allele Frequencies
                                                ssqr <- ((n_p1 * (p1_a0 - pbar_a0)^2) + (n_p2 * (p2_a0 - pbar_a0)^2)) / nbar

                                                ## Components a, b, and c
                                                # Component a
                                                a <- (ssqr - (pbar_a0 * (1 - pbar_a0) - (ssqr / 2) - (hbar / 4)) / (nbar - 1)) * nbar / nc

                                                # Component b
                                                b <- (pbar_a0 * (1 - pbar_a0) - (ssqr / 2) - hbar * (((2 * nbar) - 1) / (4 * nbar))) * nbar / (nbar - 1)

                                                # Component c
                                                c <- hbar / 2

                                                ## Fst ##
                                                fst <- a / (a + b +c)

                                                ## Sums for return ##
                                                if (!is.na(a) && !is.na(b) && !is.na(c)) {
                                                  sum1 <- sum1 + a
                                                  sum2 <- sum2 + (a + b + c)
                                                }
                                                if (!is.na(fst)) {
                                                  sum3 <- sum3 + fst
                                                }
                                                count <- count + 1

                                              }
                                              return(c(sum1, sum2, sum3, count))
                                            })

    if (weighted) {
      sum1 <- sum(sapply(batch_results, function(x) x[1]))
      sum2 <- sum(sapply(batch_results, function(x) x[2]))
      final_fst <- sum1 / sum2
    } else {
      sum3 <- sum(sapply(batch_results, function(x) x[3]))
      count <- sum(sapply(batch_results, function(x) x[4]))
      final_fst <- sum3 / count
    }
    return(final_fst)
  } else {
    # Window mode processing
    window_results <- process_vcf_in_windows(vcf_path,
                                             window_size = window_size,
                                             skip_size = skip_size,
                                             threads = threads,
                                             write_log = write_log,
                                             logfile = logfile,
                                             pop1_individuals = pop1_individuals,
                                             pop2_individuals = pop2_individuals,
                                             custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals, pop2_individuals) {
                                               # Separate populations
                                               sep <- separateByPopulations(sep_gt, pop1_names = pop1_individuals, pop2_names = pop2_individuals, rm_ref_alleles = FALSE)
                                               genotype_matrix1 <- sep$pop1
                                               genotype_matrix2 <- sep$pop2
                                               allele_freqs1 <- calculateAlleleFreqs(genotype_matrix1)
                                               allele_freqs2 <- calculateAlleleFreqs(genotype_matrix2)

                                               # Preparing variables for
                                               sum1 <- 0
                                               sum2 <- 0
                                               sum3 <- 0
                                               count <- 0

                                               # Iterate over polymorphic sites
                                               for (i in seq_len(nrow(genotype_matrix1))) {
                                                 ## Sample size values ##
                                                 # Get the number of non-missing individuals
                                                 n_p1 <- length(genotype_matrix1[i,genotype_matrix1[i,] != "."]) / 2
                                                 n_p2 <- length(genotype_matrix2[i,genotype_matrix2[i,] != "."]) / 2

                                                 # Total sum of individuals across populations
                                                 n_sum <- n_p1 + n_p2

                                                 # Average number of individuals per population
                                                 nbar <- n_sum / 2

                                                 # Sum of squares of the number of individuals
                                                 sum_nsqr <- n_p1^2 + n_p2^2

                                                 # Effective Sample Size per Population
                                                 nc <- (n_sum - (sum_nsqr / n_sum))

                                                 ## Allele frequency values ##
                                                 # Variant allele frequencies
                                                 p1_a0 <- allele_freqs1[i,1]
                                                 p1_a1 <- allele_freqs1[i,2]
                                                 p2_a0 <- allele_freqs2[i,1]
                                                 p2_a1 <- allele_freqs2[i,2]

                                                 # Total allele counts
                                                 pc_a0 <- (p1_a0 * (n_p1*2)) + (p2_a0 * (n_p2*2)) # times 2 cause diploidy
                                                 pc_a1 <- (p1_a1 * (n_p1*2)) + (p2_a1 * (n_p2*2))

                                                 # Average allele frequencies
                                                 pbar_a0 <- pc_a0 / (n_sum * 2) # times 2 cause diploidy
                                                 pbar_a1 <- pc_a1 / (n_sum * 2)

                                                 # Heterozygosity values
                                                 # Counts of heterozygotes
                                                 N_het_p1 <- 0
                                                 for (indiv_col in seq(1, ncol(genotype_matrix1), by = 2)) {
                                                   individual_genotypes <- genotype_matrix1[i, c(indiv_col, indiv_col + 1)]
                                                   is_heterozygote <- !is.na(individual_genotypes[1]) && !is.na(individual_genotypes[2]) && individual_genotypes[1] != individual_genotypes[2]
                                                   if (is_heterozygote) {
                                                     N_het_p1 <- N_het_p1 + 1
                                                   }
                                                 }
                                                 N_het_p2 <- 0
                                                 for (indiv_col in seq(1, ncol(genotype_matrix2), by = 2)) {
                                                   individual_genotypes <- genotype_matrix2[i, c(indiv_col, indiv_col + 1)]
                                                   is_heterozygote <- !is.na(individual_genotypes[1]) && !is.na(individual_genotypes[2]) && individual_genotypes[1] != individual_genotypes[2]
                                                   if (is_heterozygote) {
                                                     N_het_p2 <- N_het_p2 + 1
                                                   }
                                                 }

                                                 # Average heterozygosity
                                                 hbar <- (N_het_p1 + N_het_p2) / n_sum

                                                 # Squared Deviations of Allele Frequencies
                                                 ssqr <- ((n_p1 * (p1_a0 - pbar_a0)^2) + (n_p2 * (p2_a0 - pbar_a0)^2)) / nbar

                                                 ## Components a, b, and c
                                                 # Component a
                                                 a <- (ssqr - (pbar_a0 * (1 - pbar_a0) - (ssqr / 2) - (hbar / 4)) / (nbar - 1)) * nbar / nc

                                                 # Component b
                                                 b <- (pbar_a0 * (1 - pbar_a0) - (ssqr / 2) - hbar * (((2 * nbar) - 1) / (4 * nbar))) * nbar / (nbar - 1)

                                                 # Component c
                                                 c <- hbar / 2

                                                 ## Fst ##
                                                 fst <- a / (a + b +c)

                                                 ## Sums for return ##
                                                 if (!is.na(a) && !is.na(b) && !is.na(c)) {
                                                   sum1 <- sum1 + a
                                                   sum2 <- sum2 + (a + b + c)
                                                 }
                                                 if (!is.na(fst)) {
                                                   sum3 <- sum3 + fst
                                                 }
                                                 count <- count + 1
                                               }
                                               if (weighted) {
                                                 final_fst <- sum1 / sum2
                                               } else {
                                                 final_fst <- sum3 / count
                                               }
                                               return(c(chrom, start_pos, end_pos, final_fst))
                                             })

    # Bind results per window into a data frame
    fst_windows <- as.data.frame(do.call("rbind", window_results))
    colnames(fst_windows) <- c("Chromosome", "Start", "End", "Fst")
    return(fst_windows)
  }
}


#' OneDimSFS
#'
#' This function calculates a one-dimensional site frequency spectrum from a VCF file. It processes the file in batches for efficient memory usage.
#' The user can decide between a folded or unfolded spectrum.
#'
#' @param vcf_path Path to the VCF file.
#' @param folded Logical, deciding if folded (TRUE) or unfolded (FALSE) SFS is returned.
#' @param batch_size The number of variants to be processed in each batch
#'                  (default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param exclude_ind Optional vector of individual IDs to exclude from the analysis.
#'        If provided, the function will remove these individuals from the genotype matrix
#'        before applying the custom function. Default is NULL, meaning no individuals are excluded.
#'
#' @return Site frequency spectrum as a named vector
#'
#' @examples
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' sfs <- OneDimSFS(vcf_file, folded = FALSE)
#'
#' @export

OneDimSFS <- function(vcf_path, folded = FALSE, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt", exclude_ind = NULL) {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          exclude_ind = exclude_ind,
                                          custom_function = function(index, fix, sep_gt,pop1_individuals = NULL, pop2_individuals = NULL, ploidy = 2) {
                                            sep_gt[sep_gt == "."] <- NA  # Replace '.' with NA for missing data
                                            num_individuals <- ncol(sep_gt)

                                            # Initialize a vector to hold the site frequency spectrum for this batch
                                            sfs_batch <- numeric(num_individuals + 1)

                                            # Iterate over the sites in the genotype matrix
                                            for (i in 1:nrow(sep_gt)) {
                                              site_data <- sep_gt[i, ]
                                              valid_data <- site_data[!is.na(site_data)]  # Exclude missing data for this site
                                              derived_count <- sum(as.numeric(valid_data))  # Count the number of derived alleles

                                              # Calculate the minor allele count for folded SFS
                                              allele_count <- if (folded) {
                                                min(derived_count, length(valid_data) - derived_count)
                                              } else {
                                                derived_count
                                              }

                                              # Update the SFS for this batch
                                              sfs_batch[allele_count] <- sfs_batch[allele_count] + 1
                                            }

                                            return(sfs_batch)
                                          })

  # Aggregate SFS values from all batches
  total_sfs <- Reduce("+", lapply(batch_results, function(x) x))

  # Name the vector elements for clearer interpretation
  names(total_sfs) <- 0:(length(total_sfs) - 1)

  # For a folded SFS, remove the redundant second half of the vector
  if (folded) {
    # Determine the midpoint of the vector
    midpoint <- ceiling(length(total_sfs) / 2)
    # Keep only up to the midpoint (inclusive)
    total_sfs <- total_sfs[1:midpoint]
  }

  return(total_sfs)
}


#' TwoDimSFS
#'
#' This function calculates a two-dimensional site frequency spectrum from a VCF file for two populations. It processes the file in batches for efficient memory usage.
#' The user can decide between a folded or unfolded spectrum.
#'
#' @param vcf_path Path to the VCF file.
#' @param pop1_individuals Vector of individual names belonging to the first population.
#' @param pop2_individuals Vector of individual names belonging to the second population.
#' @param folded Logical, deciding if folded (TRUE) or unfolded (FALSE) SFS is returned.
#' @param batch_size The number of variants to be processed in each batch
#'                  (default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#' @param exclude_ind Optional vector of individual IDs to exclude from the analysis.
#'        If provided, the function will remove these individuals from the genotype matrix
#'        before applying the custom function. Default is NULL, meaning no individuals are excluded.
#'
#' @return Two-dimensional site frequency spectrum as a matrix.
#'
#' @examples
#' vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
#' index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")
#' pop1_individuals <- c("tsk_0", "tsk_1", "tsk_2")
#' pop2_individuals <- c("tsk_3", "tsk_4", "tsk_5")
#' sfs_2d <- TwoDimSFS(vcf_file, pop1_individuals, pop2_individuals, folded = TRUE)
#'
#' @export

TwoDimSFS <- function(vcf_path, pop1_individuals, pop2_individuals, folded = FALSE, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt", exclude_ind = NULL) {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          exclude_ind = exclude_ind,
                                          pop1_individuals = pop1_individuals,
                                          pop2_individuals = pop2_individuals,
                                          custom_function = function(index, fix, sep_gt,pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals, ploidy = 2) {
                                            # Separate populations
                                            sep <- separateByPopulations(sep_gt, pop1_names = pop1_individuals, pop2_names = pop2_individuals, rm_ref_alleles = FALSE)
                                            genotype_matrix1 <- sep$pop1
                                            genotype_matrix2 <- sep$pop2

                                            genotype_matrix1[genotype_matrix1 == "."] <- NA
                                            genotype_matrix2[genotype_matrix2 == "."] <- NA

                                            # Initialize a matrix to hold the 2d site frequency spectrum for this batch
                                            sfs_2d_batch <- matrix(0, nrow = length(pop1_individuals) * 2 + 1, ncol = length(pop2_individuals) * 2 + 1)
                                            x <- c()
                                            # Iterate over the sites in the genotype matrix
                                            for (i in 1:nrow(genotype_matrix1)) {
                                              site_data1 <- genotype_matrix1[i,]
                                              site_data2 <- genotype_matrix2[i,]

                                              # Exclude missing data for this site
                                              valid_data1 <- site_data1[!is.na(site_data1)]
                                              valid_data2 <- site_data2[!is.na(site_data2)]

                                              # Count the number of derived alleles (assuming '1' is the derived state)
                                              derived_count1 <- sum(as.numeric(valid_data1))
                                              derived_count2 <- sum(as.numeric(valid_data2))

                                              # Calculate the minor allele count for folded SFS
                                              if (folded) {
                                                allele_count1 <- min(derived_count1, length(valid_data1) - derived_count1)
                                                allele_count2 <- min(derived_count2, length(valid_data2) - derived_count2)
                                              } else {
                                                allele_count1 <- derived_count1
                                                allele_count2 <- derived_count2
                                              }
                                              # Update the 2dSFS
                                              sfs_2d_batch[allele_count1, allele_count2] <- sfs_2d_batch[allele_count1, allele_count2] + 1
                                            }
                                            return(sfs_2d_batch)
                                          })

  # Aggregate 2dSFS values from all batches
  total_sfs_2d <- Reduce("+", batch_results)

  # If the SFS is folded, remove the empty categories
  if (folded) {
    total_sfs_2d <- total_sfs_2d[rowSums(total_sfs_2d[,-1]) != 0,]
    total_sfs_2d <- total_sfs_2d[,colSums(total_sfs_2d[-1,]) != 0]
  }

  return(total_sfs_2d)
}

