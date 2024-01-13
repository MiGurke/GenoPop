#' meanImputation
#'
#' Imputes missing data in VCF files using mean genotypes, either by variant or by individual. This function is designed for large-scale genomic datasets and utilizes efficient batch processing and parallel computing to handle VCF files directly. The function reads a VCF file in batches, performs mean imputation, and writes the imputed data to a new VCF file.
#'
#' @param vcf_path Path to the input VCF file.
#' @param output_vcf Path for the output VCF file with imputed data.
#' @param batch_size Number of variants to process per batch (default: 10000).
#' @param threads Number of threads to use for parallel processing (default: 1).
#' @param write_log If TRUE, writes a log of the imputation process (default: FALSE).
#' @param logfile Path to the log file, used if `write_log` is TRUE.
#' @param mode Method for calculating means, either "variant" or "individual" (default: "variant").
#'
#' @details The function operates in two modes: "variant" and "individual". In "variant" mode, it calculates the mean genotype for each variant across all individuals and replaces missing values with this mean. In "individual" mode, it calculates the mean genotype for each individual across all variants and replaces missing values with this mean. The means are rounded to maintain integer genotypes. The imputed data is written to a new VCF file, maintaining the original format.
#'
#' @return Path to the output VCF file with imputed data.
#'
#' @examples
#' meanImputation("/path/to/input.vcf",
#'                "/path/to/output.vcf",
#'                batch_size = 1000,
#'                threads = 5,
#'                write_log = TRUE,
#'                logfile = "/path/to/logfile.txt",
#'                mode = "variant")
#'
#' @export

meanImputation <- function(vcf_path, output_vcf, batch_size = 10000, threads = 1, write_log = FALSE, logfile = "log.txt", mode = "variant") {
  # Other parameters and initial setup

  # Open the original VCF file to read the header
  con <- file(vcf_path, "r")
  on.exit(close(con))

  header_lines <- c()
  chrom_line <- NULL
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (startsWith(line, "##")) {
      header_lines <- c(header_lines, line)
    } else if (startsWith(line, "#CHROM")) {
      chrom_line <- line  # Save the #CHROM line separately
      break
    }
  }

  # Create the comment line about imputation
  imputation_comment <- paste("##GenoPop_MeanImputation=mode=", mode, ";date=", Sys.Date(), sep="")

  # Insert the imputation comment before the #CHROM line
  modified_header <- c(header_lines, imputation_comment, chrom_line)

  # Write the modified header to the new VCF file
  write_lines <- function(lines, path) {
    con <- file(path, "w")
    on.exit(close(con))
    writeLines(lines, con)
  }
  write_lines(modified_header, output_vcf)

  # Define the directory to store temporary files
  temp_dir <- dirname(output_vcf)

  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
                                            if (mode == "variant") {
                                              # Convert missing values ("." entries) to NA
                                              sep_gt[sep_gt == "."] <- NA

                                              # Convert the character matrix to a numeric matrix
                                              genotype_matrix <- matrix(as.numeric(sep_gt), nrow = nrow(sep_gt))

                                              # Calculate the mean for each variant (row) excluding NAs
                                              variant_means <- apply(genotype_matrix, 1, mean, na.rm = TRUE)

                                              # Create a matrix of means with the same dimensions as the original data
                                              mean_matrix <- matrix(rep(variant_means, each = ncol(genotype_matrix)),
                                                                    nrow = nrow(genotype_matrix),
                                                                    ncol = ncol(genotype_matrix),
                                                                    byrow = TRUE)

                                              # Replace NAs in the original data with corresponding means
                                              imputed <- ifelse(is.na(genotype_matrix), mean_matrix, genotype_matrix)
                                              imputed <- round(imputed)
                                              colnames(imputed) <- colnames(sep_gt)
                                            } else if (mode == "individual") {
                                              # Convert missing values ("." entries) to NA
                                              sep_gt[sep_gt == "."] <- NA

                                              # Convert the character matrix to a numeric matrix, coercing NA where appropriate
                                              genotype_matrix <- matrix(as.numeric(sep_gt), nrow = nrow(sep_gt))

                                              # Calculate the mean for each individual (column) excluding NAs
                                              indiv_means <- apply(genotype_matrix, 2, mean, na.rm = TRUE)

                                              # Create a matrix of means with the same dimensions as the original data
                                              mean_matrix <- matrix(rep(indiv_means, each = ncol(genotype_matrix)),
                                                                    nrow = nrow(genotype_matrix),
                                                                    ncol = ncol(genotype_matrix),
                                                                    byrow = TRUE)

                                              # Replace NAs in the original data with corresponding means
                                              imputed <- ifelse(is.na(genotype_matrix), mean_matrix, genotype_matrix)
                                              colnames(imputed) <- colnames(sep_gt)
                                            }

                                            # Prepare matrix for reformatting the imputed genotypes
                                            vcf_formatted_gt <- matrix(NA,
                                                                       ncol = (ncol(imputed) / 2) + 1,
                                                                       nrow = nrow(imputed))
                                            # Write the data rows
                                            for (i in seq_len(nrow(imputed))) {
                                              # Combine genotypes into the VCF genotype format
                                              gt <- apply(matrix(imputed[i, ], ncol = 2, byrow = TRUE), 1, function(g) {
                                                paste(g, collapse = "/")
                                              })
                                              vcf_formatted_gt[i, ] <- c("GT", gt)
                                            }

                                            # Combine the fix information with the imputed genotypes to get full VCF lines
                                            full_vcf_lines <- cbind(fix, vcf_formatted_gt)

                                            # Create a temporary file path in the specified directory with gzip compression
                                            temp_file <- tempfile(pattern = ".imputed_batch_", tmpdir = temp_dir, fileext = ".vcf.gz")

                                            # Write and compress the full VCF lines to the temporary file
                                            write.table(full_vcf_lines, gzfile(temp_file), col.names = FALSE, row.names = FALSE, quote = FALSE)

                                            # Return the path to the temporary file with an identifier (e.g., first position in the batch)
                                            return(list(file = temp_file, index = index))
                                          })
  # Removing NULL elements from the batch results list (NULL is returned if an error occurred in the process, error message is written to log file.)
  batch_results <- Filter(Negate(is.null), batch_results)
  # Assuming batch_results is a list of lists with 'file' and 'first_pos'
  # Sort the temporary file paths based on the first position in each batch to ensure correct order
  ordered_temp_files <- batch_results[order(sapply(batch_results, `[[`, "index"))]
  i <- 1
  for (temp_info in ordered_temp_files) {
    temp_file <- temp_info$file

    if (file.exists(temp_file)) {
      # Read the compressed imputed data from each temporary file
      imputed_data <- read.table(gzfile(temp_file))

      # Remove scientific notation from vcf position field to avoid issues with index building
      imputed_data[, 2] <- lapply(imputed_data[, 2, drop = FALSE], function(column) {
        format(column, scientific = FALSE)
      })

      # Append the imputed data to the final VCF file
      write.table(imputed_data, output_vcf, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

      # Delete the temporary file
      file.remove(temp_file)
    } else {
      message(paste0("Batch file number ", i, " does not exist, more information in the log file. Skipping."))
    }
    i <- i + 1
  }
  zipped <- bgzip(output_vcf, overwrite = TRUE)
  file.remove(output_vcf)
  indexTabix(zipped, format = "vcf")
  message(paste("Imputation completed. Imputed VCF written to:", zipped))
  return(zipped)
}

#' knn_imputeR
#'
#' This is just the wrapper function for the kNN imputation algorithm that is written in C++. Please use the function kNNImputation().
#'
#' @param data NumericMatrix, the data matrix.
#' @param k Integer, the number of neighbors.
#'
#' @return NumericMatrix, the imputed data matrix.
#' @export
#'
#' @examples
#' example_matrix <- matrix(c(0, 1, NA, 1, 1, NA, 0, 0, NA, 1, 1, NA, 0, 0, 0, 0, NA, 1, 0, 0, NA, 1, 1, NA, 0), nrow = 5, byrow = TRUE)
#' knn_imputeR(example_matrix, k = 3)

knn_imputeR <- function(data, k) {
  knn_impute(data, k)
}


#' kNNImputation
#'
#' Performs imputation of missing genomic data using a k-nearest-neighbor algorithm, optimized for large VCF files. This function reads VCF files in batches, applies kNN imputation, and writes the results to a new VCF file. It's designed for efficient processing through parallel computing and chunk-based handling of SNP data. The function uses the 'annoy' library for efficient neighbor detection. Note that neighboring SNPs lacking data will be excluded from the imputation process, and positions with large proportions of missing data may remain unimputed. The choice of batch size, number of neighbors (k), and the number of iterations (maxiter) is critical for balancing accuracy and computational demand.
#'
#' @param vcf_path Path to the input VCF file.
#' @param output_vcf Path for the output VCF file with imputed data.
#' @param k Number of nearest neighbors used for imputation (default: 10).
#' @param maxiter Maximum number of iterations for the kNN algorithm (default: 3).
#' @param batch_size Number of variants to process per batch (default: 1000).
#' @param threads Number of threads used for computation (default: 1).
#' @param write_log If TRUE, writes a log file of the process (advised for large datasets).
#' @param logfile Path to the log file, used if `write_log` is TRUE.
#'
#' @return Path to the output VCF file with imputed data.
#'
#' @examples
#' kNNImputation("/path/to/input.vcf",
#'               "/path/to/output.vcf",
#'               k = 10,
#'               maxiter = 3,
#'               batch_size = 1000,
#'               threads = 5,
#'               write_log = TRUE,
#'               logfile = "/path/to/logfile.txt")
#'
#' @import Rcpp
#' @export

kNNImputation <- function(vcf_path, output_vcf, k = 10, maxiter = 3, batch_size = 1000, threads = 1, write_log = FALSE, logfile = "logfile.txt") {
  # Other parameters and initial setup

  # Open the original VCF file to read the header
  con <- file(vcf_path, "r")
  on.exit(close(con))

  header_lines <- c()
  chrom_line <- NULL
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (startsWith(line, "##")) {
      header_lines <- c(header_lines, line)
    } else if (startsWith(line, "#CHROM")) {
      chrom_line <- line  # Save the #CHROM line separately
      break
    }
  }

  # Create the comment line about imputation
  imputation_comment <- paste("##GenoPop_kNNImputation=k=", k, ";maxiter=", maxiter,";date=", Sys.Date(), sep="")

  # Insert the imputation comment before the #CHROM line
  modified_header <- c(header_lines, imputation_comment, chrom_line)

  # Write the modified header to the new VCF file
  write_lines <- function(lines, path) {
    con <- file(path, "w")
    on.exit(close(con))
    writeLines(lines, con)
  }
  write_lines(modified_header, output_vcf)

  # Define the directory to store temporary files
  temp_dir <- dirname(output_vcf)

  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          add_packages = "Rcpp",
                                          custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
                                            # Convert missing values (".") to NA and prepare the matrix
                                            sep_gt[sep_gt == "."] <- NA
                                            numeric_gt <- matrix(as.numeric(sep_gt), nrow = nrow(sep_gt), ncol = ncol(sep_gt))

                                            # kNN imputation via the c++ code
                                            imputed <- knn_impute(numeric_gt, k, maxiter)

                                            # Change NA back to "."
                                            imputed[is.na(imputed)] <- "."

                                            # Prepare matrix for reformatting the imputed genotypes
                                            vcf_formatted_gt <- matrix(NA,
                                                                       ncol = (ncol(imputed) / 2) + 1,
                                                                       nrow = nrow(imputed))

                                            # Write the data rows
                                            for (i in seq_len(nrow(imputed))) {
                                              # Combine genotypes into the VCF genotype format
                                              gt <- apply(matrix(imputed[i, ], ncol = 2, byrow = TRUE), 1, function(g) {
                                                paste(g, collapse = "/")
                                              })
                                              vcf_formatted_gt[i, ] <- c("GT", gt)
                                            }

                                            # Combine the fix information with the imputed genotypes to get full VCF lines
                                            full_vcf_lines <- cbind(fix, vcf_formatted_gt)

                                            # Create a temporary file path in the specified directory with gzip compression
                                            temp_file <- tempfile(pattern = ".imputed_batch_", tmpdir = temp_dir, fileext = ".vcf.gz")

                                            # Write and compress the full VCF lines to the temporary file
                                            write.table(full_vcf_lines, gzfile(temp_file), col.names = FALSE, row.names = FALSE, quote = FALSE)

                                            # Return the path to the temporary file with an identifier (e.g., first position in the batch)
                                            return(list(file = temp_file, index = index))
                                          })
  # Removing NULL elements from the batch results list (NULL is returned if an error occurred in the process, error message is written to log file.)
  batch_results <- Filter(Negate(is.null), batch_results)
  # Assuming batch_results is a list of lists with 'file' and 'first_pos'
  # Sort the temporary file paths based on the first position in each batch to ensure correct order
  ordered_temp_files <- batch_results[order(sapply(batch_results, `[[`, "index"))]
  i <- 1
  for (temp_info in ordered_temp_files) {
    temp_file <- temp_info$file

    if (file.exists(temp_file)) {
      # Read the compressed imputed data from each temporary file
      imputed_data <- read.table(gzfile(temp_file))

      # Remove scientific notation from vcf position field to avoid issues with index building
      imputed_data[, 2] <- lapply(imputed_data[, 2, drop = FALSE], function(column) {
        format(column, scientific = FALSE)
      })

      # Append the imputed data to the final VCF file
      write.table(imputed_data, output_vcf, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

      # Delete the temporary file
      file.remove(temp_file)
    } else {
      message(paste0("Batch file number ", i, " does not exist, more information in the log file. Skipping."))
    }
    i <- i + 1
  }
  zipped <- bgzip(output_vcf, overwrite = TRUE)
  file.remove(output_vcf)
  indexTabix(zipped, format = "vcf")
  message(paste("Imputation completed. Imputed VCF written to:", zipped))
  return(zipped)
}

kNNImputation <- function(vcf_path, output_vcf, k = 10, maxiter = 3, batch_size = 1000, threads = 1, write_log = FALSE, logfile = "logfile.txt") {
  # Other parameters and initial setup

  # Open the original VCF file to read the header
  con <- file(vcf_path, "r")
  on.exit(close(con))

  header_lines <- c()
  chrom_line <- NULL
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (startsWith(line, "##")) {
      header_lines <- c(header_lines, line)
    } else if (startsWith(line, "#CHROM")) {
      chrom_line <- line  # Save the #CHROM line separately
      break
    }
  }

  # Create the comment line about imputation
  imputation_comment <- paste("##GenoPop_kNNImputation=k=", k, ";maxiter=", maxiter,";date=", Sys.Date(), sep="")

  # Insert the imputation comment before the #CHROM line
  modified_header <- c(header_lines, imputation_comment, chrom_line)

  # Write the modified header to the new VCF file
  write_lines <- function(lines, path) {
    con <- file(path, "w")
    on.exit(close(con))
    writeLines(lines, con)
  }
  write_lines(modified_header, output_vcf)

  # Define the directory to store temporary files
  temp_dir <- dirname(output_vcf)

  batch_results <- process_vcf_in_batches(vcf_path,
                                          batch_size = batch_size,
                                          threads = threads,
                                          write_log = write_log,
                                          logfile = logfile,
                                          add_packages = "Rcpp",
                                          custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
                                            # Convert missing values (".") to NA and prepare the matrix
                                            sep_gt[sep_gt == "."] <- NA
                                            numeric_gt <- matrix(as.numeric(sep_gt), nrow = nrow(sep_gt), ncol = ncol(sep_gt))

                                            # kNN imputation via the c++ code
                                            imputed <- knn_impute(numeric_gt, k, maxiter)

                                            # Change NA back to "."
                                            imputed[is.na(imputed)] <- "."

                                            # Prepare matrix for reformatting the imputed genotypes
                                            vcf_formatted_gt <- matrix(NA,
                                                                       ncol = (ncol(imputed) / 2) + 1,
                                                                       nrow = nrow(imputed))

                                            # Write the data rows
                                            for (i in seq_len(nrow(imputed))) {
                                              # Combine genotypes into the VCF genotype format
                                              gt <- apply(matrix(imputed[i, ], ncol = 2, byrow = TRUE), 1, function(g) {
                                                paste(g, collapse = "/")
                                              })
                                              vcf_formatted_gt[i, ] <- c("GT", gt)
                                            }

                                            # Combine the fix information with the imputed genotypes to get full VCF lines
                                            full_vcf_lines <- cbind(fix, vcf_formatted_gt)

                                            # Create a temporary file path in the specified directory with gzip compression
                                            temp_file <- tempfile(pattern = ".imputed_batch_", tmpdir = temp_dir, fileext = ".vcf.gz")

                                            # Write and compress the full VCF lines to the temporary file
                                            write.table(full_vcf_lines, gzfile(temp_file), col.names = FALSE, row.names = FALSE, quote = FALSE)

                                            # Return the path to the temporary file with an identifier (e.g., first position in the batch)
                                            return(list(file = temp_file, index = index))
                                          })
  # Removing NULL elements from the batch results list (NULL is returned if an error occurred in the process, error message is written to log file.)
  batch_results <- Filter(Negate(is.null), batch_results)
  # Assuming batch_results is a list of lists with 'file' and 'first_pos'
  # Sort the temporary file paths based on the first position in each batch to ensure correct order
  ordered_temp_files <- batch_results[order(sapply(batch_results, `[[`, "index"))]
  i <- 1
  for (temp_info in ordered_temp_files) {
    temp_file <- temp_info$file

    if (file.exists(temp_file)) {
      # Read the compressed imputed data from each temporary file
      imputed_data <- read.table(gzfile(temp_file))

      # Remove scientific notation from vcf position field to avoid issues with index building
      imputed_data[, 2] <- lapply(imputed_data[, 2, drop = FALSE], function(column) {
        format(column, scientific = FALSE)
      })

      # Append the imputed data to the final VCF file
      write.table(imputed_data, output_vcf, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

      # Delete the temporary file
      file.remove(temp_file)
    } else {
      message(paste0("Batch file number ", i, " does not exist, more information in the log file. Skipping."))
    }
    i <- i + 1
  }
  zipped <- bgzip(output_vcf, overwrite = TRUE)
  file.remove(output_vcf)
  indexTabix(zipped, format = "vcf")
  message(paste("Imputation completed. Imputed VCF written to:", zipped))
  return(zipped)
}

#' rfImputation
#'
#' Performs imputation of missing genomic data using a random forest algorithm, optimized for large-scale VCF files. This function reads VCF files in batches, applies random forest imputation via the missForest package, and writes the results to a new VCF file. The choice of batch size, number of improvement iterations, and the number of decision trees in the random forest is critical for balancing accuracy and computational demand.
#'
#' @param vcf_path Path to the input VCF file.
#' @param output_vcf Path for the output VCF file with imputed data.
#' @param batch_size Number of variants to process per batch (default: 1000).
#' @param maxiter Number of improvement iterations for the random forest algorithm (default: 10).
#' @param ntree Number of decision trees in the random forest (default: 100).
#' @param threads Number of threads used for computation (default: 1, or one less than available on the system).
#' @param write_log If TRUE, writes a log file of the process (advised for large datasets).
#' @param logfile Path to the log file, used if `write_log` is TRUE.
#'
#' @return Path to the output VCF file with imputed data.
#'
#' @examples
#' rfImputation("/path/to/input.vcf",
#'              "/path/to/output.vcf",
#'              batch_size = 1000,
#'              maxiter = 10,
#'              ntree = 100,
#'              threads = 5,
#'              write_log = TRUE,
#'              logfile = "/path/to/logfile.txt")
#'
#' @importFrom missForest missForest
#' @importFrom Rsamtools bgzip indexTabix
#' @export

rfImputation <- function(vcf_path, output_vcf, batch_size = 1000, maxiter = 10, ntree = 100, threads = 1, write_log = FALSE, logfile = "log.txt") {
  # Other parameters and initial setup

  # Open the original VCF file to read the header
  con <- file(vcf_path, "r")
  on.exit(close(con))

  header_lines <- c()
  chrom_line <- NULL
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (startsWith(line, "##")) {
      header_lines <- c(header_lines, line)
    } else if (startsWith(line, "#CHROM")) {
      chrom_line <- line  # Save the #CHROM line separately
      break
    }
  }

  # Create the comment line about imputation
  imputation_comment <- paste("##GenoPop_rfImputation=maxiter=", maxiter, ";ntree=", ntree, ";date=", Sys.Date(), sep="")

  # Insert the imputation comment before the #CHROM line
  modified_header <- c(header_lines, imputation_comment, chrom_line)

  # Write the modified header to the new VCF file
  write_lines <- function(lines, path) {
    con <- file(path, "w")
    on.exit(close(con))
    writeLines(lines, con)
  }
  write_lines(modified_header, output_vcf)

  # Define the directory to store temporary files
  temp_dir <- dirname(output_vcf)

  batch_results <- process_vcf_in_batches(vcf_path,
                         batch_size = batch_size,
                         threads = threads,
                         write_log = write_log,
                         logfile = logfile,
                         add_packages = "missForest",
                         custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL) {
                           # Convert missing values (".") to NA and prepare the matrix
                           sep_gt[sep_gt == "."] <- NA
                           numeric_gt <- matrix(as.numeric(sep_gt), nrow = nrow(sep_gt), ncol = ncol(sep_gt))

                           # Perform missForest imputation on the batch
                           imputed <- missForest(numeric_gt, maxiter = maxiter, ntree = ntree)$ximp

                           # Ensure imputed values are integers
                           imputed <- matrix(as.character(round(as.numeric(imputed))),
                                                   ncol = ncol(imputed),
                                                   nrow = nrow(imputed))

                           if (ncol(imputed) != ncol(sep_gt)) {
                             e <- simpleError("One or more individuals have only missing Genotypes, batch not imputed.")
                             stop(e)
                           }

                           # Prepare matrix for reformatting the imputed genotypes
                           vcf_formatted_gt <- matrix(NA,
                                                      ncol = (ncol(imputed) / 2) + 1,
                                                      nrow = nrow(imputed))
                           # Write the data rows
                           for (i in seq_len(nrow(imputed))) {
                             # Combine genotypes into the VCF genotype format
                             gt <- apply(matrix(imputed[i, ], ncol = 2, byrow = TRUE), 1, function(g) {
                               paste(g, collapse = "/")
                             })
                             vcf_formatted_gt[i, ] <- c("GT", gt)
                           }

                           # Combine the fix information with the imputed genotypes to get full VCF lines
                           full_vcf_lines <- cbind(fix, vcf_formatted_gt)

                           # Create a temporary file path in the specified directory with gzip compression
                           temp_file <- tempfile(pattern = ".imputed_batch_", tmpdir = temp_dir, fileext = ".vcf.gz")

                           # Write and compress the full VCF lines to the temporary file
                           write.table(full_vcf_lines, gzfile(temp_file), col.names = FALSE, row.names = FALSE, quote = FALSE)

                           # Return the path to the temporary file with an identifier (e.g., first position in the batch)
                           return(list(file = temp_file, index = index))
                         })
  # Removing NULL elements from the batch results list (NULL is returned if an error occurred in the process, error message is written to log file.)
  batch_results <- Filter(Negate(is.null), batch_results)
  # Assuming batch_results is a list of lists with 'file' and 'first_pos'
  # Sort the temporary file paths based on the first position in each batch to ensure correct order
  ordered_temp_files <- batch_results[order(sapply(batch_results, `[[`, "index"))]
  i <- 1
  for (temp_info in ordered_temp_files) {
    temp_file <- temp_info$file

    if (file.exists(temp_file)) {
      # Read the compressed imputed data from each temporary file
      imputed_data <- read.table(gzfile(temp_file))

      # Remove scientific notation from vcf position field to avoid issues with index building
      imputed_data[, 2] <- lapply(imputed_data[, 2, drop = FALSE], function(column) {
        format(column, scientific = FALSE)
      })

      # Append the imputed data to the final VCF file
      write.table(imputed_data, output_vcf, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

      # Delete the temporary file
      file.remove(temp_file)
    } else {
      message(paste0("Batch file number ", i, " does not exist, more information in the log file. Skipping."))
    }
    i <- i + 1
  }
  zipped <- bgzip(output_vcf, overwrite = TRUE)
  file.remove(output_vcf)
  indexTabix(zipped, format = "vcf")
  message(paste("Imputation completed. Imputed VCF written to:", zipped))
  return(zipped)
}
