#' Internal Batch Processing of VCF Data
#'
#' This function is designed for efficient internal processing of large VCF files.
#' It reads the VCF file in batches, processes each batch in parallel, and applies a custom function
#' to the processed data. It's optimized for performance with large genomic datasets and is not
#' intended to be used directly by users.
#'
#' @param vcf_path The path to the VCF file.
#' @param batch_size The number of variants to be processed in each batch.
#' @param custom_function A custom function that takes two arguments: a data frame similar to
#'        the `@fix` slot of a `vcfR` object and a genotype matrix similar to the `@sep_gt` slot.
#'        This function is applied to each batch of data.
#' @param threads The number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile The path to the log file.
#' @param exclude_ind Optional vector of individual IDs to exclude from the analysis.
#'        If provided, the function will remove these individuals from the genotype matrix
#'        before applying the custom function. Default is NULL, meaning no individuals are excluded.
#'
#' @details The function divides the VCF file into batches based on the specified `batch_size`.
#'          Each batch is processed to extract fixed information and genotype data, filter for
#'          biallelic SNPs, and then apply the `custom_function`. The function uses parallel
#'          processing to improve efficiency and can handle large genomic datasets.
#'
#'          This function is part of the internal workings of the package and is not intended
#'          to be called directly by the user. It is documented for the sake of completeness
#'          and to assist in maintenance and understanding of the package's internal mechanics.
#'
#' @return A data frame that is the combined result of applying the `custom_function` to each batch.
#'         The structure of this data frame depends on the `custom_function` used.
#'
#' @keywords internal
#'
#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom GenomicRanges GRanges
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom missForest missForest
#' @importFrom IRanges IRanges
#'
#' @noRd


process_vcf_in_batches <- function(vcf_path, batch_size, custom_function, threads = 1, write_log = FALSE, logfile = "logfile.txt", pop1_individuals = NULL, pop2_individuals = NULL, add_packages = NULL, exclude_ind = NULL) {
  tbx <- open(TabixFile(vcf_path, yieldSize=batch_size))
  # Initialize variables
  batch_coordinates <- list()
  index <- 1
  batch_coord <- NULL
  while(length(res <- scanTabix(tbx)[[1]])) {
    # Get all chromosomes in the current batch
    chroms <- sapply(res, function(x) strsplit(x, "\t")[[1]][1])
    unique_chroms <- unique(chroms)

    # Iterate over each unique chromosome and find start and end positions
    for (chrom in unique_chroms) {
      chrom_indices <- which(chroms == chrom)
      first_index <- min(chrom_indices)
      last_index <- max(chrom_indices)

      first_variant_info <- strsplit(res[first_index], "\t")[[1]]
      last_variant_info <- strsplit(res[last_index], "\t")[[1]]

      start_pos <- first_variant_info[2]
      end_pos <- last_variant_info[2]

      # Construct and store the genomic region for each chromosome segment
      batch_coordinates <- c(batch_coordinates, paste0(index, "\t", chrom, "\t", start_pos, "\t", end_pos))
      index <- index + 1
    }
  }
  # Close the Tabix file
  close(tbx)

  # Get individual names:
  # Open the VCF file
  con <- gzfile(vcf_path, "r")
  on.exit(close(con))

  # Read lines until the header line starting with #CHROM
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (startsWith(line, "#CHROM")) {
      break
    }
  }

  # Split the line and extract individual names
  line <- strsplit(line, "\t")
  individual_names <- line[[1]][10:length(line[[1]])]

  # Determine the number of cores
  if (is.null(threads)) {
    num_cores <- detectCores() - 1  # Reserve one core for the system
  } else {
    num_cores <- threads
  }

  # Set up the parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("calculateAlleleFreqs", "separateByPopulations"))

  # Prepare log file if write_log is true
  #Create log file and prepare progress tracking if write log is true
  if (write_log) {
    message(paste0("Preparing log file ", logfile))
    # Function to safely write to a log file
    log_progress <- function(message, log_file) {
      # Obtain a lock on the file to avoid write conflicts
      fileConn <- file(log_file, open = "a")
      tryCatch({
        # Write the progress message
        writeLines(message, con = fileConn)
      }, finally = {
        # Release the file lock
        close(fileConn)
      })
    }
    file.create(logfile)
  }

  if (!is.null(add_packages)) {
    packages = c("Rsamtools", "GenomicRanges", add_packages)
  } else {
    packages = c("Rsamtools", "GenomicRanges")
  }

  # Perform the calculations in parallel
  results <- foreach(batch_coord = batch_coordinates, .packages = packages) %dopar% {
    tryCatch({
      # Extract chromosome and positions from the batch coordinates
      coords <- strsplit(batch_coord, "\t")[[1]]
      index <- as.numeric(coords[1])
      chrom <- coords[2]
      start_pos <- as.numeric(coords[3])
      end_pos <- as.numeric(coords[4])

      # Read the specific batch from the VCF file
      tbx <- TabixFile(vcf_path)
      batch_data <- scanTabix(tbx, param = GRanges(chrom, IRanges(start = start_pos, end = end_pos)))

      # Split each line of the batch data by tabs and create a matrix
      data_matrix <- do.call(rbind, strsplit(batch_data[[1]], "\t"))
      rm(batch_data)

      # Create a data frame from the matrix, selecting the relevant columns
      fix <- data.frame(
        CHROM = data_matrix[, 1],
        POS = as.numeric(data_matrix[, 2]),
        ID = data_matrix[, 3],
        REF = data_matrix[, 4],
        ALT = data_matrix[, 5],
        QUAL = data_matrix[, 6],
        FILTER = data_matrix[, 7],
        INFO = data_matrix[, 8],
        stringsAsFactors = FALSE
      )

      # Remove SNP's that are not biallelic or complex
      biallelic_indices <- which(nchar(as.character(fix$ALT)) == 1)
      fix <- fix[biallelic_indices, ]

      # Assuming genotype data is in the 10th column onwards in VCF format
      gt_matrix <- data_matrix[biallelic_indices, 10:ncol(data_matrix)]
      rm(data_matrix)

      # Check if exclude_ind is provided and not NULL
      ex_ind_names <- individual_names
      if (!is.null(exclude_ind)) {
        # Check for errors in individual names
        if (length(which(!(exclude_ind %in% individual_names))) > 0) {
          wrong <- which(!(exclude_ind %in% individual_names))
          e <- simpleError(paste0("Individual name '", exclude_ind[wrong], "' not found in VCF file."))
          stop(e)
        }
        # Find columns (individuals) to exclude
        cols_to_exclude <- which(individual_names %in% exclude_ind)
        # Exclude the specified individuals from the genotype matrix
        gt_matrix <- gt_matrix[, -cols_to_exclude, drop = FALSE]
        ex_ind_names <- individual_names[-cols_to_exclude]
      }

      # If pop1 and pop2 individuals are given, also check them for errors
      if (!is.null(pop1_individuals)) {
        # Check for errors in individual names
        if (length(which(!(pop1_individuals %in% individual_names))) > 0) {
          wrong <- which(!(pop1_individuals %in% individual_names))
          e <- simpleError(paste0("Individual name '", pop1_individuals[wrong], "' not found in VCF file."))
          stop(e)
        }
      }
      if (!is.null(pop2_individuals)) {
        # Check for errors in individual names
        if (length(which(!(pop2_individuals %in% individual_names))) > 0) {
          wrong <- which(!(pop2_individuals %in% individual_names))
          e <- simpleError(paste0("Individual name '", pop2_individuals[wrong], "' not found in VCF file."))
          stop(e)
        }
      }

      # Detect separators for alleles (commonly '/' or '|')
      allele_separators <- unique(gsub("[^/|]", "", gt_matrix))
      separator <- allele_separators[1]  # Assuming consistent use of a single separator

      # Escape the pipe character if it's the separator
      if (separator == "|") {
        separator <- "\\|"
      }

      # Estimate ploidy level from the genotype data
      example_gt <- gt_matrix[1, 1]  # Using the first variant as an example
      ploidy <- length(unlist(strsplit(example_gt, separator)))

      # Separate alleles into different columns
      sep_gt <- matrix(NA, nrow = nrow(gt_matrix), ncol = ploidy * ncol(gt_matrix))

      for (i in seq_len(nrow(gt_matrix))) {
        # Check if there are multiple fields in the genotype data
        if (any(grepl(":", gt_matrix[i, ]))) {
          # Extract the GT part (assumes it's the first field)
          gt_matrix[i, ] <- sapply(strsplit(gt_matrix[i, ], ":"), `[`, 1)
        }
        # Separate alleles for each variant and assign to the matrix
        sep_gt[i, ] <- unlist(strsplit(gt_matrix[i, ], separator))
      }

      # Assign column names (e.g., "Sample1_1", "Sample1_2" for diploid)
      colnames(sep_gt) <- paste(rep(ex_ind_names, each = ploidy),
                                rep(1:ploidy, times = length(ex_ind_names)),
                                sep = "_")

      # Removing rows now with only missing data or monomorphic
      rows_to_remove <- apply(gt_matrix, 1, function(x) all(x == "./." | x == ".|."))
      gt_matrix <- gt_matrix[!rows_to_remove, ]
      sep_gt <- sep_gt[!rows_to_remove, ]
      fix <- fix[!rows_to_remove, ]

      # Apply the custom function to the batch data
      process_result <- custom_function(index, fix, sep_gt, pop1_individuals, pop2_individuals, ploidy)

      if (write_log) {
        log_progress(paste0(Sys.time(), " Completed batch ", index, ": ", chrom, ":", start_pos, "-", end_pos, "\n"), logfile)
      }

      return(process_result)

    }, error = function(e) {
      # Return or log the error information
      if (write_log) {
        log_progress(paste0(Sys.time(), " Error in batch ", index, ": ", chrom, ":", start_pos, "-", end_pos," Error: ", e$message, "\n"), logfile)
      }
      return(NULL)  # Return NULL or some error indication
    })
  }

  # Stop the cluster
  stopCluster(cl)

  return(results)
}

#' Internal Window-Based Processing of VCF Data
#'
#' This function is designed for efficient internal processing of large VCF files
#' based on specified window sizes and genomic skips. It reads the VCF file in windows,
#' processes each window in parallel, and applies a custom function to the processed data.
#' It's optimized for performance with large genomic datasets and is not intended
#' to be used directly by end users.
#'
#' @param vcf_path The path to the VCF file.
#' @param window_size The genomic length of each window in base pairs.
#' @param skip_size The size of the genomic region to skip between windows in base pairs.
#' @param custom_function A custom function that takes two arguments: a data frame similar to
#'        the `@fix` slot of a `vcfR` object and a genotype matrix similar to the `@sep_gt` slot.
#'        This function is applied to each window of data.
#' @param threads The number of threads to use for parallel processing.
#' @param write_log Logical, indicating whether to write progress logs.
#' @param logfile The path to the log file.
#'
#' @details The function divides the VCF file into windows based on the specified `window_size`
#'          and `skip_size`. Each window is processed to extract fixed information and genotype data,
#'          filter for biallelic SNPs, and then apply the `custom_function`. The function uses parallel
#'          processing to improve efficiency and can handle large genomic datasets.
#'
#'          This function is part of the internal workings of the package and is not intended
#'          to be called directly by the user. It is documented for the sake of completeness
#'          and to assist in maintenance and understanding of the package's internal mechanics.
#'
#' @return A list that is the combined result of applying the `custom_function` to each window.
#'         The structure of this list depends on the `custom_function` used.
#'
#' @keywords internal
#'
#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom GenomicRanges GRanges
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom IRanges IRanges
#'
#' @noRd

process_vcf_in_windows <- function(vcf_path, window_size, skip_size, custom_function, threads = 1, write_log = FALSE, logfile = "logfile.txt", pop1_individuals = NULL, pop2_individuals = NULL, add_packages = NULL, exclude_ind = NULL) {
  # Read the VCF header to extract chromosome information and individual names
  con <- gzfile(vcf_path, "r")
  on.exit(close(con))
  chrom_info <- list()
  individual_names <- NULL
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (startsWith(line, "##contig=<ID=")) {
      chrom_id <- gsub(".*ID=([^,]+),.*", "\\1", line)
      length <- as.numeric(gsub(".*length=([0-9]+).*", "\\1", line))
      chrom_info[[chrom_id]] <- length
    } else if (startsWith(line, "#CHROM")) {
      line <- strsplit(line, "\t")
      individual_names <- line[[1]][10:length(line[[1]])]
      break
    }
  }
  # Initialize variables for window coordinates
  window_coordinates <- list()
  index <- 1
  window_coord <- NULL
  # Generate window coordinates for each chromosome
  for (chrom in names(chrom_info)) {
    chrom_length <- chrom_info[[chrom]]
    start_pos <- 1
    while (start_pos < chrom_length) {
      end_pos <- min(start_pos + window_size - 1, chrom_length)
      window_coordinates <- c(window_coordinates, paste0(index, "\t", chrom, "\t", start_pos, "\t", end_pos))
      start_pos <- start_pos + window_size + skip_size
      index <- index + 1
    }
  }

  # Determine the number of cores
  if (is.null(threads)) {
    num_cores <- detectCores() - 1  # Reserve one core for the system
  } else {
    num_cores <- threads
  }
  # Set up the parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("calculateAlleleFreqs", "separateByPopulations"))

  # Prepare log file if write_log is true
  # Create log file and prepare progress tracking if write log is true
  if (write_log) {
    message(paste0("Preparing log file ", logfile))
    # Function to safely write to a log file
    log_progress <- function(message, log_file) {
      # Obtain a lock on the file to avoid write conflicts
      fileConn <- file(log_file, open = "a")
      tryCatch({
        # Write the progress message
        writeLines(message, con = fileConn)
      }, finally = {
        # Release the file lock
        close(fileConn)
      })
    }
    file.create(logfile)
  }

  if (!is.null(add_packages)) {
    packages = c("Rsamtools", "GenomicRanges", add_packages)
  } else {
    packages = c("Rsamtools", "GenomicRanges")
  }

  # Process each window in parallel
  results <- foreach(window_coord = window_coordinates, .packages = packages) %dopar% {
    tryCatch({
      # Extract chromosome and positions from the batch coordinates
      coords <- strsplit(window_coord, "\t")[[1]]
      index <- as.numeric(coords[1])
      chrom <- coords[2]
      start_pos <- as.numeric(coords[3])
      end_pos <- as.numeric(coords[4])

      # Read the specific batch from the VCF file
      tbx <- TabixFile(vcf_path)
      batch_data <- scanTabix(tbx, param = GRanges(chrom, IRanges(start = start_pos, end = end_pos)))

      # Check if the window is empty
      if (length(batch_data[[1]]) == 0) {
        # Return an empty result or NULL, depending on how you want to handle it
        process_result <- NULL
      } else {

        # Split each line of the batch data by tabs and create a matrix
        data_matrix <- do.call(rbind, strsplit(batch_data[[1]], "\t"))
        rm(batch_data)

        # Create a data frame from the matrix, selecting the relevant columns
        fix <- data.frame(
          CHROM = data_matrix[, 1],
          POS = as.numeric(data_matrix[, 2]),
          ID = data_matrix[, 3],
          REF = data_matrix[, 4],
          ALT = data_matrix[, 5],
          QUAL = data_matrix[, 6],
          FILTER = data_matrix[, 7],
          INFO = data_matrix[, 8],
          stringsAsFactors = FALSE
        )

        # Remove SNP's that are not biallelic
        biallelic_indices <- which(nchar(as.character(fix$ALT)) == 1)
        fix <- fix[biallelic_indices, ]

        # Assuming genotype data is in the 10th column onwards in VCF format
        gt_matrix <- data_matrix[biallelic_indices, 10:ncol(data_matrix)]
        rm(data_matrix)

        # Check if exclude_ind is provided and not NULL
        ex_ind_names <- individual_names
        if (!is.null(exclude_ind)) {
          # Check for errors in individual names
          if (length(which(!(exclude_ind %in% individual_names))) > 0) {
            wrong <- which(!(exclude_ind %in% individual_names))
            e <- simpleError(paste0("Individuals names '", exclude_ind[wrong], "' not found in VCF file."))
            stop(e)
          }
        }
        # If pop1 and pop2 individuals are given, also check them for errors
        if (!is.null(pop1_individuals)) {
          # Check for errors in individual names
          if (length(which(!(pop1_individuals %in% individual_names))) > 0) {
            wrong <- which(!(pop1_individuals %in% individual_names))
            e <- simpleError(paste0("Individual name '", pop1_individuals[wrong], "' not found in VCF file."))
            stop(e)
          }
        }
        if (!is.null(pop2_individuals)) {
          # Check for errors in individual names
          if (length(which(!(pop2_individuals %in% individual_names))) > 0) {
            wrong <- which(!(pop2_individuals %in% individual_names))
            e <- simpleError(paste0("Individual name '", pop2_individuals[wrong], "' not found in VCF file."))
            stop(e)
          }
        }
        # Detect separators for alleles (commonly '/' or '|')
        allele_separators <- unique(gsub("[^/|]", "", gt_matrix))
        separator <- allele_separators[1]  # Assuming consistent use of a single separator

        # Escape the pipe character if it's the separator
        if (separator == "|") {
          separator <- "\\|"
        }

        # Estimate ploidy level from the genotype data
        example_gt <- gt_matrix[1, 1]  # Using the first variant as an example
        ploidy <- length(unlist(strsplit(example_gt, separator)))

        # Separate alleles into different columns
        sep_gt <- matrix(NA, nrow = nrow(gt_matrix), ncol = ploidy * ncol(gt_matrix))

        for (i in seq_len(nrow(gt_matrix))) {
          # Check if there are multiple fields in the genotype data
          if (any(grepl(":", gt_matrix[i, ]))) {
            # Extract the GT part (assumes it's the first field)
            gt_matrix[i, ] <- sapply(strsplit(gt_matrix[i, ], ":"), `[`, 1)
          }
          # Separate alleles for each variant and assign to the matrix
          sep_gt[i, ] <- unlist(strsplit(gt_matrix[i, ], separator))
        }

        # Assign column names (e.g., "Sample1_1", "Sample1_2" for diploid)
        colnames(sep_gt) <- paste(rep(ex_ind_names, each = ploidy),
                                  rep(1:ploidy, times = length(ex_ind_names)),
                                  sep = "_")

        # Removing rows now with only missing data or monomorphic
        rows_to_remove <- apply(gt_matrix, 1, function(x) all(x == "./." | x == ".|."))
        gt_matrix <- gt_matrix[!rows_to_remove, ]
        sep_gt <- sep_gt[!rows_to_remove, ]
        fix <- fix[!rows_to_remove, ]

        # Apply the custom function to the batch data
        process_result <- custom_function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals, pop2_individuals)

        if (write_log) {
          log_progress(paste0(Sys.time(), " Completed window ", index, ": ", chrom, ":", start_pos, "-", end_pos, "\n"), logfile)
        }
        return(process_result)

      }
    }, error = function(e) {
      # Return or log the error information
      if (write_log) {
        log_progress(paste0(Sys.time(), " Error in window ", index, ": ", chrom, ":", start_pos, "-", end_pos," Error: ", e$message, "\n"), logfile)
      }
      return(NULL)  # Return NULL or some error indication
    })
  }

  # Clean up and return results
  stopCluster(cl)
  return(results)
}


#' Separate Genotype Matrix by Populations
#'
#' This function separates a genotype matrix into two data frames based on population assignments.
#' It's designed to work with the batches or windows processed by `process_vcf_in_batches` and `process_vcf_in_windows`.
#'
#' @param sep_gt A genotype matrix similar to the `@sep_gt` slot of a `vcfR` object.
#' @param pop1_names A character vector of individual names for the first population.
#' @param pop2_names A character vector of individual names for the second population.
#' @param rm_ref_alleles Logical, whether variants that only have the reference allele
#'        should be removed from the respective subpopulations data frame. (Default = TRUE)
#'
#' @return A list containing two data frames, one for each population.
#'
#' @keywords internal
#'
#' @export

separateByPopulations <- function(sep_gt, pop1_names, pop2_names, ploidy = 2, rm_ref_alleles = TRUE) {

  colnames1 <- paste(rep(pop1_names, each = ploidy),
                     rep(1:ploidy, times = length(pop1_names)),
                     sep = "_")
  colnames2 <- paste(rep(pop2_names, each = ploidy),
                     rep(1:ploidy, times = length(pop2_names)),
                     sep = "_")

  # Create data frames for each population
  pop1_gt <- sep_gt[, colnames1, drop = FALSE]
  pop2_gt <- sep_gt[, colnames2, drop = FALSE]

  if (rm_ref_alleles) {
    non_ref_rows1 <- apply(pop1_gt, 1, function(x) !all(x %in% c("0", ".")))
    pop1_gt <- pop1_gt[non_ref_rows1, , drop = FALSE]
    non_ref_rows2 <- apply(pop2_gt, 1, function(x) !all(x %in% c("0", ".")))
    pop2_gt <- pop2_gt[non_ref_rows2, , drop = FALSE]
  }

  return(list(pop1 = pop1_gt, pop2 = pop2_gt))
}


#' Calculate Allele Frequencies from Genotype Matrix
#'
#' This function calculates allele frequencies from a genotype matrix (sep_gt) for each variant.
#' It is designed to be used within the batch or window processing framework of `process_vcf_in_batches` and `process_vcf_in_windows`.
#'
#' @param sep_gt Genotype matrix similar to the `@sep_gt` slot of a `vcfR` object.
#'
#' @return A data frame containing allele frequencies for each variant.
#'
#' @keywords internal
#'
#' @export

calculateAlleleFreqs <- function(sep_gt) {
  allele_frequencies_per_site <- vector("list", length = nrow(sep_gt))
  unique_alleles <- c("0", "1")

  # Calculate allele frequencies for each variant
  for (i in seq_len(nrow(sep_gt))) {
    alleles <- sep_gt[i, ]
    # Removing missing genotypes from the calculation all over
    alleles <- alleles[alleles != "."]
    allele_counts <- table(factor(alleles, levels = unique_alleles))
    allele_frequencies <- allele_counts / sum(allele_counts)
    allele_frequencies_per_site[[i]] <- allele_frequencies
  }

  # Combine the frequencies into a data frame
  allele_frequencies_df <- do.call(rbind, allele_frequencies_per_site)
  colnames(allele_frequencies_df) <- unique_alleles

  return(allele_frequencies_df)
}
