### ADD WINDOW MODE TO EACH OF THEM ###

#' Count Fixed Sites for Alternative Allele in VCF File
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
#' vcf_path <- "path/to/vcf/file"
#' num_fixed_sites <- FixedSites(vcf_path, threads = 4, write_log = TRUE, logfile = "fixed_sites_log.txt")
#'
#' # Window mode example
#' vcf_path <- "path/to/vcf/file"
#' fixed_sites_df <- FixedSites(vcf_path, threads = 4, write_log = TRUE, logfile = "windowed_fixed_sites_log.txt", window_size = 100000, skip_size = 50000)
#'
#' @export


FixedSites <- function(vcf_path, threads = 1, write_log = FALSE, logfile = "log.txt", batch_size = 10000, window_size = NULL, skip_size = NULL) {
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
                                            custom_function = function(index, fix, sep_gt,pop1_individuals = NULL, pop2_individuals = NULL) {
                                              allele_freqs <- calculateAlleleFreqs(sep_gt)
                                              # Count fixed sites for the alternative allele in this batch
                                              return(sum(apply(allele_freqs, 1, function(x) any(x[2] == 1))))  # -1 to exclude the reference allele
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
                                            custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos,pop1_individuals = NULL, pop2_individuals = NULL) {
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


#' Count Segregating Sites in VCF File
#'
#' This function counts the number of polymorphic sites (sites not fixed for the alternative allele)
#' in a VCF file. It processes the file in batches or specified windows across the genome.
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows.
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
#'
#' @return
#' In batch mode (no window_size or skip_size provided): A single integer representing the total number of polymorphic sites across the entire VCF file.
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'PolymorphicSites', representing the count of polymorphic sites within each window.
#'
#' @examples
#' # Batch mode example
#' vcf_path <- "path/to/vcf/file"
#' num_polymorphic_sites <- SegregatingSites(vcf_path, threads = 4, write_log = TRUE, logfile = "polymorphic_sites_log.txt")
#'
#' # Window mode example
#' vcf_path <- "path/to/vcf/file"
#' polymorphic_sites_df <- SegregatingSites(vcf_path, threads = 4, write_log = TRUE, logfile = "windowed_polymorphic_sites_log.txt", window_size = 100000, skip_size = 50000)
#'
#' @export

SegregatingSites <- function(vcf_path, threads = 1, write_log = FALSE, logfile = "log.txt", batch_size = 10000, window_size = NULL, skip_size = NULL) {
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
                                            custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
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
                                             custom_function = function(fix, sep_gt, chrom, start_pos, end_pos,pop1_individuals = NULL, pop2_individuals = NULL) {
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


#' Count Singleton Sites in VCF File
#'
#' This function counts the number of singleton sites (sites where a minor allele occurs only once in the sample)
#' in a VCF file. It processes the file in batches or specified windows across the genome.
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows.
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
#'
#' @return
#' In batch mode (no window_size or skip_size provided): A single integer representing the total number of singleton sites across the entire VCF file.
#' In window mode (window_size and skip_size provided): A data frame with columns 'Chromosome', 'Start', 'End', and 'SingletonSites', representing the count of singleton sites within each window.
#'
#' @examples
#' # Batch mode example
#' vcf_path <- "path/to/vcf/file"
#' num_singleton_sites <- SingletonSites(vcf_path, threads = 4, write_log = TRUE, logfile = "singleton_sites_log.txt")
#'
#' # Window mode example
#' vcf_path <- "path/to/vcf/file"
#' singleton_sites_df <- SingletonSites(vcf_path, threads = 4, write_log = TRUE, logfile = "windowed_singleton_sites_log.txt", window_size = 100000, skip_size = 50000)
#'
#' @export

SingletonSites <- function(vcf_path, threads = 1, write_log = FALSE, logfile = "log.txt", batch_size = 10000, window_size = NULL, skip_size = NULL) {
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
                                            custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
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
                                             custom_function = function(fix, sep_gt, chrom, start_pos, end_pos,pop1_individuals = NULL, pop2_individuals = NULL) {
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


#' Count Private Alleles in VCF File
#'
#' This function calculates the number of private alleles in two populations from a VCF file.
#' It processes the file in batches or specified windows across the genome.
#' For batch processing, it uses `process_vcf_in_batches`. For windowed analysis, it uses a similar
#' approach tailored to process specific genomic windows.
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
#' vcf_path <- "path/to/vcf/file"
#' pop1_individuals <- c("8449", "8128", "8779")
#' pop2_individuals <- c("8816", "8823", "8157")
#' private_alleles <- PrivateAlleles(vcf_path, pop1_individuals, pop2_individuals, threads = 4, write_log = TRUE, logfile = "private_alleles_log.txt")
#'
#' # Window mode example
#' private_alleles_windows <- PrivateAlleles(vcf_path, pop1_individuals, pop2_individuals, threads = 4, write_log = TRUE, logfile = "windowed_private_alleles_log.txt", window_size = 100000, skip_size = 50000)
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
                                            custom_function = function(index, fix, sep_gt,pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals) {
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
                                             logfile = logfile,
                                             custom_function = function(fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals, pop2_individuals) {
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
    private_alleles_windows <-as.data.frame(do.call("rbind", window_results))
    colnames(df) <- c("Chromosome", "Start", "End", "PrivateAllelesPop1", "PrivateAllelesPop2")
    return(private_alleles_windows)
  }
}


#' Calculate Observed Heterozygosity from VCF File
#'
#' This function calculates the observed heterozygosity (Ho) for each variant in a VCF file.
#' It processes the file in batches for efficient memory usage.
#'
#' @param vcf_path Path to the VCF file.
#' @param batch_size The number of variants to be processed in each batch
#'                  (default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#'
#' @return Observed heterozygosity averaged over all loci.
#'
#' @examples
#' vcf_path <- "path/to/your/vcf/file"
#' Ho <- ObservedHeterozygosity(vcf_path)
#'
#' @export

ObservedHeterozygosity <- function(vcf_path, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt") {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
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

                                            # Return the count of heterozygotes and the number of valid loci per individual
                                            return(c(heterozygotes, num_individuals * nrow(sep_gt)))
                                          })

  # Combine the counts from all batches and calculate the mean
  all_heterozygotes <- sum(sapply(batch_results, function(x) x[1]))
  all_valid_loci <- sum(sapply(batch_results, function(x) x[2]))
  Ho <- all_heterozygotes / all_valid_loci

  return(Ho)
}


#' Calculate Expected Heterozygosity from VCF File
#'
#' This function calculates the expected heterozygosity (He) for each variant in a VCF file.
#' It processes the file in batches for efficient memory usage.
#'
#' @param vcf_path Path to the VCF file.
#' @param batch_size The number of variants to be processed in each batch
#'                  (default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#'
#' @return Expected heterozygosity averaged over all loci.
#'
#' @examples
#' vcf_path <- "path/to/your/vcf/file"
#' He <- ExpectedHeterozygosity(vcf_path)
#'
#' @export

ExpectedHeterozygosity <- function(vcf_path, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt") {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
                                            # Calculate allele frequencies for this batch
                                            allele_freqs <- calculateAlleleFreqs(sep_gt)
                                            # Calculate expected heterozygosity per locus for this batch
                                            He_per_locus <- apply(allele_freqs, 1, function(p) {
                                              1 - sum(p^2)
                                            })
                                            return(He_per_locus)
                                          })

  # Combine the He per locus from all batches and calculate the mean
  all_He_per_locus <- unlist(batch_results)
  He <- mean(all_He_per_locus, na.rm = TRUE)

  return(He)
}


#' Calculate Nucleotide Diversity (Pi) from VCF File
#'
#' This function calculates the nucleotide diversity (Pi) for a VCF file. Nei & Li, 1979 (https://doi.org/10.1073/pnas.76.10.5269).
#' It processes the file in batches for efficient memory usage and requires the overall sequence length.
#' The formula used for this is equivalent to the one used in vcftools --window-pi (https://vcftools.sourceforge.net/man_latest.html).
#' Handling missing alleles at one site is equivalent to Korunes & Samuk, 2021 ( https://doi.org/10.1111/1755-0998.13326), but for simplicity assuming that completely missing sites are invariant sites, which will underestimate Pi.
#' Otherwise this would only function with VCF files that include all monomorphic sites, which may be unpractical given common data sets.
#' If you happen to know the number of missing sites vs the number of monomorphic sites, please use the number of monomorphic + the number of polymorphic sites as the sequence length to get the most accurate estimation of Pi.
#'
#' @param vcf_path Path to the VCF file.
#' @param seq_length Total length of the sequence in number of bases.
#' @param batch_size The number of variants to be processed in each batch
#'                  (default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#'
#' @return Nucleotide diversity (Pi) across the sequence.
#'
#' @examples
#' vcf_path <- "path/to/your/vcf/file"
#' total_sequence_length <- 265392  # Total length of the sequence
#' pi_value <- Pi(vcf_path, total_sequence_length)
#'
#' @export

Pi <- function(vcf_path, seq_length, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt") {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
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
}


#' Calculate Tajima's D from VCF File
#'
#' This function calculates Tajima's D statistic for a given dataset. It processes the file in batches for efficient memory usage.
#' Tajima's D is a measure for neutrality based on the difference between the number of segregating sites and the average number of nucleotide differences.
#'
#' @param vcf_path Path to the VCF file.
#' @param seq_length Total length of the sequence in number of bases.
#' @param batch_size The number of variants to be processed in each batch
#'                  (default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#'
#' @return Tajima's D value.
#'
#' @examples
#' vcf_path <- "path/to/your/vcf/file"
#' total_sequence_length <- 265392  # Total length of the sequence
#' tajimas_d <- TajimasD(vcf_path, total_sequence_length)
#'
#' @export

TajimasD <- function(vcf_path, seq_length, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt") {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
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

  # Including monomorphic sites
  N_monomorphic_sites <- seq_length - all_variants_counted
  pi <- total_pi * seq_length  # Adjust total pi by sequence length to get average pi per site

  # Tajima's D calculation
  denominator <- sqrt((e1 * total_S) + (e2 * total_S * (total_S - 1)))
  D <- (pi - (total_S / a1)) / denominator

  return(D)
}


#' Calculate Watterson's Theta from VCF File
#'
#' This function calculates Watterson's Theta, a measure for neutrality, from a VCF file. It processes the file in batches for efficient memory usage.
#' The metric will be normalized by the sequence length to make it comparable between datasets.
#'
#' @param vcf_path Path to the VCF file.
#' @param seq_length The length of the sequence in the data set.
#' @param batch_size The number of variants to be processed in each batch
#'                  (default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#'
#' @return Watterson's theta value normalized by the sequence length.
#'
#' @examples
#' vcf_path <- "path/to/your/vcf/file"
#' total_sequence_length <- 265392  # Total length of the sequence
#' wattersons_theta <- WattersonsTheta(vcf_path, total_sequence_length)
#'
#' @export

WattersonsTheta <- function(vcf_path, seq_length, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt") {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
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
}


#' Calculate Average Nucleotide Differences (Dxy) from VCF File
#'
#' This function calculates the average number of nucleotide differences per site (Dxy) between two populations from a VCF file.  Nei & Li, 1979 (https://doi.org/10.1073/pnas.76.10.52699).
#' Handling missing alleles at one site is equivalent to Korunes & Samuk, 2021 ( https://doi.org/10.1111/1755-0998.13326), but for simplicity assuming that completely missing sites are invariant sites, which will underestimate Dxy.
#' Otherwise this would only function with VCF files that include all monomorphic sites, which may be unpractical given common data sets.
#' If you happen to know the number of missing sites vs the number of monomorphic sites, please just use the number of monomorphic + the number of polymorphic sites as the sequence length to the the most accurate estimation of Dxy.
#'
#' @param vcf_path Path to the VCF file.
#' @param pop1_individuals Vector of individual names belonging to the first population.
#' @param pop2_individuals Vector of individual names belonging to the second population.
#' @param seq_length Length of the sequence in number of bases, including monomorphic sites.
#' @param batch_size The number of variants to be processed in each batch
#'                  (default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#'
#' @return The average number of nucleotide substitutions per site between the individuals of two populations (Dxy).
#'
#' @examples
#' vcf_path <- "path/to/your/vcf/file"
#' pop1_individuals <- c("8449", "8128", "8779")
#' pop2_individuals <- c("8816", "8823", "8157")
#' total_sequence_length <- 265392  # Total length of the sequence
#' dxy_value <- Dxy(vcf_path, pop1_individuals, pop2_individuals, total_sequence_length)
#'
#' @export

Dxy <- function(vcf_path, pop1_individuals, pop2_individuals, seq_length, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt") {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          pop1_individuals = pop1_individuals,
                                          pop2_individuals = pop2_individuals,
                                          custom_function = function(index, fix, sep_gt,pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals) {
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
}


#' Calculate Fst from VCF File
#'
#' This function calculates the fixation index (Fst) between two populations from a VCF file using the method of Weir and Cockerham (1984). It processes the file in batches for efficient memory usage.
#' The result can be either the mean or weighted Fst based on the number of non-missing chromosomes.
#'
#' @param vcf_path Path to the VCF file.
#' @param pop1_individuals Vector of individual names belonging to the first population.
#' @param pop2_individuals Vector of individual names belonging to the second population.
#' @param weighted Logical, whether weighted Fst or mean Fst is returned (Default = FALSE (mean Fst is returned)).
#' @param batch_size The number of variants to be processed in each batch
#'                  (default of 10,000 should be suitable for most use cases).
#' @param threads Number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile Path to the log file where progress will be logged.
#'
#' @return Fst value (either mean or weighted).
#'
#' @examples
#' vcf_path <- "path/to/your/vcf/file"
#' pop1_individuals <- c("8449", "8128", "8779")
#' pop2_individuals <- c("8816", "8823", "8157")
#' fst_value <- Fst(vcf_path, pop1_individuals, pop2_individuals, weighted = TRUE)
#'
#' @export

Fst <- function(vcf_path, pop1_individuals, pop2_individuals, weighted = FALSE, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt") {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          pop1_individuals = pop1_individuals,
                                          pop2_individuals = pop2_individuals,
                                          custom_function = function(index, fix, sep_gt,pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals) {
                                            # Separate populations
                                            sep <- separateByPopulations(sep_gt, pop1_names = pop1_individuals, pop2_names = pop2_individuals, rm_ref_alleles = FALSE)
                                            genotype_matrix1 <- sep$pop1
                                            genotype_matrix2 <- sep$pop2
                                            allele_freqs1 <- calculateAlleleFreqs(genotype_matrix1)
                                            allele_freqs2 <- calculateAlleleFreqs(genotype_matrix2)

                                            # Calculate average allele frequencies across populations
                                            mean_freqs <- (allele_freqs1 + allele_freqs2) / 2

                                            # Initialize variables for Fst calculation
                                            Fsts <- c()
                                            weights <- c()

                                            # Iterate over polymorphic sites
                                            for (i in seq_len(nrow(genotype_matrix1))) {
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
                                              weights <- c(weights, (n1+n2))
                                              Fsts <- c(Fsts, Fst_locus)
                                            }
                                            return(list(Fsts = Fsts, weights = weights))
                                          })

  # Aggregate Fst values from all batches
  all_Fsts <- unlist(lapply(batch_results, function(x) x$Fsts))
  all_weights <- unlist(lapply(batch_results, function(x) x$weights))

  # Calculate final Fst value
  if (weighted) {
    fst_weighted_products <- all_Fsts * all_weights
    sum_weighted_products <- sum(fst_weighted_products, na.rm = TRUE)
    sum_weights <- sum(all_weights, na.rm = TRUE)
    Fst_final <- sum_weighted_products / sum_weights
  } else {
    Fst_final <- mean(all_Fsts, na.rm = TRUE)
  }

  return(Fst_final)
}


#' Calculate One Dimensional Site Frequency Spectrum from VCF File
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
#'
#' @return Site frequency spectrum as a named vector
#'
#' @examples
#' vcf_path <- "path/to/your/vcf/file"
#' sfs <- OneDimSFS(vcf_path, folded = FALSE)
#'
#' @export

OneDimSFS <- function(vcf_path, folded = FALSE, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt") {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
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


#' Calculate Two-Dimensional Site Frequency Spectrum from VCF File
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
#'
#' @return Two-dimensional site frequency spectrum as a matrix.
#'
#' @examples
#' vcf_path <- "path/to/your/vcf/file"
#' pop1_individuals <- c("8449", "8128", "8779")
#' pop2_individuals <- c("8816", "8823", "8157")
#' sfs_2d <- TwoDimSFS(vcf_path, pop1_individuals, pop2_individuals, folded = TRUE)
#'
#' @export

TwoDimSFS <- function(vcf_path, pop1_individuals, pop2_individuals, folded = FALSE, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt") {
  # Batch mode processing
  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          pop1_individuals = pop1_individuals,
                                          pop2_individuals = pop2_individuals,
                                          custom_function = function(index, fix, sep_gt,pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals) {
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

