#' GenoPop-Impute
#'
#' Performs imputation of missing genomic data in batches using the missForest (Stekhoven & Bühlmanm, 2012) algorithm. This function reads VCF files, divides it into batches of a fixed number of SNPs, applies the missForest algorithm to each batch, and writes the results to a new VCF file, which will be returned bgzipped and tabix indexed. The choice of the batch size is critical for balancing accuracy and computational demand. We found that a batch size of 500 SNPs is the most accurate for recombination rates typical of mammalians. For on average higher recombination rates (> 5 cM/Mb) we recommend a batch size of 100 SNPs.
#'
#' @param vcf_path Path to the input VCF file.
#' @param output_vcf Path for the output VCF file with imputed data.
#' @param batch_size Number of SNPs to process per batch (default: 500).
#' @param maxiter Number of improvement iterations for the random forest algorithm (default: 10).
#' @param ntree Number of decision trees in the random forest (default: 100).
#' @param threads Number of threads used for computation (default: 1).
#' @param write_log If TRUE, writes a log file of the process (advised for large datasets).
#' @param logfile Path to the log file, used if `write_log` is TRUE.
#'
#' @return Path to the output VCF file with imputed data.
#'
#' @importFrom missForest missForest
#' @importFrom Rsamtools bgzip indexTabix
#' @importFrom utils read.table write.table
#' @export
#'
#' @examples
#'  vcf_file <- system.file("tests/testthat/sim_miss.vcf.gz", package = "GenoPop")
#'  index_file <- system.file("tests/testthat/sim_miss.vcf.gz.tbi", package = "GenoPop")
#'  output_file <- tempfile(fileext = ".vcf")
#'  GenoPop_Impute(vcf_file, output_vcf = output_file, batch_size = 500)




GenoPop_Impute <- function(vcf_path, output_vcf, batch_size = 1000, maxiter = 10, ntree = 100, threads = 1, write_log = FALSE, logfile = "log.txt") {
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
                         custom_function = function(index, fix, sep_gt, pop1_individuals = NULL, pop2_individuals = NULL, ploidy = 2) {
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
                                                      ncol = (ncol(imputed) / ploidy) + 1,
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
