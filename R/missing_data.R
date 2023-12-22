### Chnage to work with process_vcf_in_batches and write new vcf ###
#' rmMissingData
#'
#' Remove variants from vcfR object with to much missing data.
#'
#' @param object A S4 object of class vcfR.
#' @param threshold Fraction of missing individuals per variant that is still accepted. Default: 0.1
#'
#' @return A S4 object of the same class, but without variants that did not meet the missingness threshold. In addition, the slot "missing_data" will be added to the object. It contains a list with information about the removed variants and the missingness per variant and individual.
#'
#' @examples
#' data("real", package = "GenoPop")
#' vcf <- rmMissingData(real, 0.1)
#'
#' @export

rmMissingData <- function(object, threshold = 0.1) {
  # Check if needed slots are available, if not create them
  if (!inherits(object, "GPvcfR") | nrow(object@sep_gt) == 0) {
    message("Calculating ploidy and separating genotype data...")
    object <- calculatePloidyAndSepGT(object)
  }

  separated_gt_data <- object@sep_gt
  gt_data <- object@gt
  fix_data <- object@fix
  ploidy <- object@ploidy
  num_inds <- ncol(separated_gt_data)
  # Manage Missing Data
  # Identify missing data
  missing_data_indices <- which(separated_gt_data == ".", arr.ind = TRUE)

  if (threshold < 1) {
    # If missing data is present, manage it
    if (length(missing_data_indices) > 0) {
      # Count missing data per variant (row)
      missing_data_count <- table(missing_data_indices[, 1])

      # Identify variants (rows) exceeding the missing data threshold
      high_missing_data_variants <- as.numeric(names(missing_data_count[missing_data_count / (ploidy * num_inds) > threshold]))

      # Exclude variants (rows) with too much missing data
      separated_gt_data <- separated_gt_data[-high_missing_data_variants, ]
      gt_data <- gt_data[-high_missing_data_variants, ]
      fix_data <- fix_data[-high_missing_data_variants, ]

    } else {
      missing_data_count <- NULL
      high_missing_data_variants <- NULL
    }

    # Calculate and output missing data per individual
    individual_indices <- ceiling(missing_data_indices[, 2] / ploidy)
    missing_data_per_individual <- table(individual_indices) / 2
    names(missing_data_per_individual) <- colnames(gt_data)[as.numeric(names(missing_data_per_individual))]

    missing_data = list(
      count_per_variant = missing_data_count,
      fraction_per_variant = missing_data_count / ncol(gt_data),
      excluded_variants = high_missing_data_variants,
      num_of_excluded_variants = length(high_missing_data_variants),
      count_per_individual = missing_data_per_individual,
      fraction_per_individual = missing_data_per_individual / nrow(gt_data)
    )

    object@missing_data <- missing_data
    object@sep_gt <- separated_gt_data
    object@gt <- gt_data
    object@fix <- fix_data

    message(paste0(as.character(length(high_missing_data_variants)), " SNPs removed because more then a fraction of ", as.character(threshold), " data was missing."))
  } else {
    # If missing data is present, manage it
    if (length(missing_data_indices) > 0) {
      # Count missing data per variant (row)
      missing_data_count <- table(missing_data_indices[, 1])

    } else {
      missing_data_count <- NULL
      high_missing_data_variants <- NULL
    }

    # Calculate and output missing data per individual
    individual_indices <- ceiling(missing_data_indices[, 2] / ploidy)
    missing_data_per_individual <- table(individual_indices) / 2
    names(missing_data_per_individual) <- colnames(gt_data)[as.numeric(names(missing_data_per_individual))]

    missing_data = list(
      count_per_variant = missing_data_count,
      fraction_per_variant = missing_data_count / ncol(gt_data),
      count_per_individual = missing_data_per_individual,
      fraction_per_individual = missing_data_per_individual / nrow(gt_data)
    )

    object@missing_data <- missing_data

  }

   return(object)
}

### Change to work with process_vcf_in_batches and write new vcf ###
#' meanImputation
#'
#' Use the mean over each variant or individual to replace missing data with that mean. Mean is rounded to keep genotype integers, so this just corresponds to the most often occuring genotyp. If you want to use this algorithm on your data, please use the \code{\link{imputeMissingData}} function which will do the operation on GPvcfR object.
#'
#' @param sep_gt separated genotype matrix from myvcfR object.
#' @param mode Means are calculated either yb "variant" or by "individual".
#'
#' @return Matrix with imputed missing data.
#'
#' @examples
#' example_matrix <- matrix(c("0", "1", ".", "1", "1", ".", "0", "0", ".", "1", "1", ".", "0", "0", "0", "0", ".", "1", "0", "0", ".", "1", "1", ".", "0"), nrow = 5, byrow = TRUE)
#' meanImputation(example_matrix, mode ="variant")
#' meanImputation(example_matrix, mode ="individual")
#'
#' @export

meanImputation <- function(sep_gt, mode = "variant") {

  if (mode == "variant") {
    # Convert missing values ("." entries) to NA
    sep_gt[sep_gt == "."] <- NA

    # Count NAs before imputation
    na_count_before <- sum(is.na(sep_gt))

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
    imputed_matrix <- ifelse(is.na(genotype_matrix), mean_matrix, genotype_matrix)
    colnames(imputed_matrix) <- colnames(sep_gt)
  } else if (mode == "individual") {
    # Convert missing values ("." entries) to NA
    sep_gt[sep_gt == "."] <- NA

    # Count NAs before imputation
    na_count_before <- sum(is.na(sep_gt))

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
    imputed_matrix <- ifelse(is.na(genotype_matrix), mean_matrix, genotype_matrix)
    colnames(imputed_matrix) <- colnames(sep_gt)
  }
  # Display message
  message(paste0("mean was able to impute ", as.character(na_count_before)), " missing genotypes." )
  # Round to nearest integer if necessary
  imputed_matrix <- round(imputed_matrix)
  return(imputed_matrix)
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


### Change to work with process_vcf_in_batches and write new vcf ###
#' kNNImputation
#'
#' Function to execute the paralellized imputation of missing data using a k-nearest-neighbor algorithm. Imputation is done in chunks of SNP's though the genome. The size of the chunks needs to be chosen carefully, as larger chunks may give more accuracy to an extend (assuming that region very far apart in the genome are likely not neighbors any way, because they should be more different from closer regions), but will increase computation demand drastically. This implementation of the algorithm uses the annoy library (https://github.com/spotify/annoy) to detect the neighbors more efficiently. Neighboring SNP's with no data for the individual to be imputed, will be exculded. This can lead to some positions left missing, if there are individuals with large proportions of missing data. The number of neighbours k, also needs to be chosen wisely. Larger k's might give more accuracy but will also increase the computational needs, although not as drastically as for the chunk size. *I will carry out some more formal tests on this algortihm soon and will include more information about this here soon.* If you want to use this algorithm on your data, please use the \code{\link{imputeMissingData}} function which will do the operation on GPvcfR object.
#'
#' @param sep_gt A separated genotype matrix from a GPvcfR object.
#' @param k Number of nearest neighbours used for imputation, default: 3.
#' @param chunk_size Number of variants analyzed in on batch in the parallelization. Default: 1000. Increasing this might improve accuracy, but will substantially increase running time.
#' @param threads Number of threads used for the computation. Default is one less then available on the system.
#' @param write_log Logical, whether a log file of the process should be written to disk. This is adviced for imputing large data sets.
#' @param logfile Name of the log file, if write_log is true.
#'
#' @return A separated genotype matrix from a GPvcfR object, but with imputed missing values.
#' @export
#'
#' @examples
#' example_matrix <- matrix(c("0", "1", ".", "1", "1", ".", "0", "0", ".", "1", "1", ".", "0", "0", "0", "0", ".", "1", "0", "0", ".", "1", "1", ".", "0"), nrow = 5, byrow = TRUE)
#' kNNImputation(example_matrix, k = 3, chunk_size = 1000)
#' @export

kNNImputation <- function(sep_gt, k = 3, chunk_size = 1000, threads = NULL, write_log = FALSE, logfile = "logfile.txt") {
  # Convert missing values ("." entries) to NA
  sep_gt[sep_gt == "."] <- NA
  # Count NAs before imputation
  na_count_before <- sum(is.na(sep_gt))

  # Convert the character matrix to a numeric matrix, coercing NA where appropriate
  genotype_matrix <- matrix(as.numeric(sep_gt), nrow = nrow(sep_gt))

  # Determine the number of chunks and cores to use
  num_rows <- nrow(genotype_matrix)
  num_chunks <- ceiling(num_rows / chunk_size)
  # Determine the number of cores
  if (is.null(threads)){
    num_cores <- min(detectCores() - 1, num_chunks) # Reserve one core for the system
  } else {
    num_cores <- min(threads, num_chunks)
  }

  # Splitting matrix into chunks
  chunks <- list()
  for (i in seq_len(num_chunks)) {
    start_row <- (i - 1) * chunk_size + 1
    end_row <- min(i * chunk_size, num_rows)
    chunks[[i]] <- genotype_matrix[start_row:end_row, , drop = FALSE]
  }


  # Initiate parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

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
    # Adding identifiers to each matrix
    chunks <- lapply(1:length(chunks), function(i) {
      list(id = i, matrix = chunks[[i]])
    })
  }

  # Split matrix and perform imputation in parallel
  imputed_chunks <- foreach(chunk = chunks, .packages = "GenoPop") %dopar% {
    # Computation with and without log file
    if(write_log){
      mc <- chunk$matrix
      id <- chunk$id
      result <- knn_impute(mc, k)
      log_progress(paste0(Sys.time(), " Completed task ", id, " of ", num_chunks, "."), logfile)
    } else {
      result <- knn_impute(chunk, k)
    }
    return(result)
  }

  # Stop the cluster
  stopCluster(cl)

  # Combine the chunks back into a single matrix
  imputed_matrix <- do.call(rbind, imputed_chunks)

  imputed_matrix[is.na(imputed_matrix)] <- "."

  # Count NAs after imputation
  na_count_after <- sum(imputed_matrix == ".")

  # Display message
  message(paste0("kNN was able to impute ", as.character(na_count_before)), " missing genotypes." )

  colnames(imputed_matrix) <- colnames(sep_gt)
  return(imputed_matrix)
}

### Change to work with process_vcf_in_batches and write new vcf ###
#' rfImputation
#'
#' Missing data imputation using the random forest algorithm implemented in missForest R package. Computation is parallelized. Imputation is done in chunks of SNP's though the genome. The size of the chunks needs to be chosen carefully, as larger chunks may give more accuracy to an extend (assuming that region very far apart in the genome are likely not neighbors any way, because they should be more different from closer regions), but will increase computation demand drastically. *I will carry out some more formal tests on this algortihm soon and will include more information about this here soon.* If you want to use this algorithm on your data, please use the \code{\link{imputeMissingData}} function which will do the operation on GPvcfR object.
#'
#' @param sep_gt A separated genotype matrix from a myvcfR object.
#' @param maxiter The number of improvement iterations the random forest algorithm (missForest) runs.
#' @param ntree The number of decision trees in the random forest.
#' @param chunk_size Number of variants analyzed in on batch in the parallelization. Default: 1000. Increasing this might improve accuracy, but will substantially increase running time.
#' @param threads Number of threads used for the computation. Default is one less then available on the system.
#' @param write_log Logical, whether a log file of the process should be written to disk. This is adviced for imputing large data sets.
#' @param logfile Name of the log file, if write_log is true.
#'
#' @return A separated genotype matrix from a myvcfR object, but with imputed missing values.
#'
#' @importFrom missForest missForest
#' @importFrom Rsamtools bgzip indexTabix
#'
#' @examples
#' example_matrix <- matrix(c("0", "1", ".", "1", "1", ".", "0", "0", ".", "1", "1", ".", "0", "0", "0", "0", ".", "1", "0", "0", ".", "1", "1", ".", "0"), nrow = 5, byrow = TRUE)
#' rfImputation(example_matrix, maxiter = 10, ntree = 100, chunk_size = 1000)
#'
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

                           # Ensure imputed values are integers (0, 1, 2, etc.)
                           imputed <- round(imputed)
                           imputed[is.na(imputed)] <- 0  # Optionally handle any remaining NAs

                           # Convert the imputed matrix back to the VCF genotype format
                           # This step will depend on how your VCF encodes genotypes and may need adjustment
                           vcf_formatted_gt <- apply(imputed, c(1,2), function(x) paste0(x, "|", x))  # Simple example for diploid

                           # Combine the fix information with the imputed genotypes to get full VCF lines
                           full_vcf_lines <- cbind(fix, vcf_formatted_gt)

                           # Create a temporary file path in the specified directory with gzip compression
                           temp_file <- tempfile(pattern = ".imputed_batch_", tmpdir = temp_dir, fileext = ".vcf.gz")

                           # Write and compress the full VCF lines to the temporary file
                           write.table(full_vcf_lines, gzfile(temp_file), col.names = FALSE, row.names = FALSE, quote = FALSE)

                           # Return the path to the temporary file with an identifier (e.g., first position in the batch)
                           return(list(file = temp_file, index = index))
                         })

  # Assuming batch_results is a list of lists with 'file' and 'first_pos'
  # Sort the temporary file paths based on the first position in each batch to ensure correct order
  ordered_temp_files <- batch_results[order(sapply(batch_results, `[[`, "index"))]

  for (temp_info in ordered_temp_files) {
    temp_file <- temp_info$file
    # Read the compressed imputed data from each temporary file
    imputed_data <- read.table(gzfile(temp_file))

    # Append the imputed data to the final VCF file
    write.table(imputed_data, output_vcf, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

    # Delete the temporary file
    file.remove(temp_file)
  }
  zipped <- bgzip(output_vcf, overwrite = TRUE)
  file.remove(output_vcf)
  indexTabix(zipped, format = "vcf")
  message(paste("Imputation completed. Imputed VCF written to:", zipped))
  return(zipped)
}
