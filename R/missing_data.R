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

#' imputeMissingData
#'
#' Impute missing variants in genotype data stored in vcfR object.
#'
#' @param object A S4 object of class vcfR.
#' @param method Method used for missing data imputation. Available are "kNN", "rf", and "mean". (Default = "mean").
#' @param ... Additional parameters for different imputation methods. For more info look up the documentation of them:\code{\link{meanImputation}}, \code{\link{kNNImputation}}, \code{\link{rfImputation}}.
#'
#' @return A S4 object of the same class, but the slot imp_gt is now filled with the imputed genotype matrix.
#'
#' @examples
#' data("real", package = "GenoPop")
#' vcf <- imputeMissingData(real, method = "mean")
#'
#' @export

imputeMissingData <- function(object, method = "mean", ...) {
 # This is not removing anything, just utilizing that functions ability to prepare the data and get stats about the amount of missing data.
 object <- rmMissingData(object, threshold = 1)
 sep_gt <- object@sep_gt
 if (method == "mean") {
   # Extract additional arguments for 'mean'
   mean_args <- list(...)
   # Check if 'threshold' is provided, if not use a default value
   if (!"mode" %in% names(mean_args)) {
     message("No mode for mean imputation provided. Will use variant means per default.")
     mean_args$mode <- "variant"  # Default value
   }
  imputed_matrix <- meanImputation(sep_gt, mode = mean_args$mode)
 } else if (method == "kNN") {
   kNN_args <- list(...)
   if (!"k" %in% names(kNN_args)) {
     message("No k for kNN imputation provided. Will use k=3 per default.")
     kNN_args$k <- 3
   }
   if (!"write_log" %in% names(kNN_args)) {
     kNN_args$write_log <- FALSE
   }
   if (!"logfile" %in% names(kNN_args)) {
     kNN_args$logfile <- "log_kNN.txt"
   }
   if (!"chunk_size" %in% names(kNN_args)) {
     message("No chunk size provided for kNN imputation. Will use chunk_size=1000 per default.")
     kNN_args$chunk_size <- 1000
   }
   if (!"threads" %in% names(kNN_args)) {
     message("No number of threads provided for kNN imputation. Will use only 1 then.")
     kNN_args$threads <- 1
   }
  imputed_matrix <- kNNImputation(sep_gt, k = kNN_args$k, chunk_size = kNN_args$chunk_size, write_log = kNN_args$write_log, logfile = kNN_args$logfile, threads = kNN_args$threads)
 } else if (method == "rf") {
   rf_args <- list(...)
   if (!"maxiter" %in% names(rf_args)){
     message("No max number of refinment iterations provided for random forest imputation. Will use maxiter=10 per default.")
     rf_args$maxiter <- 10
   }
   if (!"ntree" %in% names(rf_args)){
     message("No number of trees in the forest provided for random forest imputation. Will use ntree=100 per default.")
     rf_args$ntree <- 100
   }
   if (!"chunk_size" %in% names(rf_args)){
     message("No chunk size provided for random forest imputation. Will use chunk_size=1000 per default.")
     rf_args$chunk_size <- 1000
   }
   if (!"write_log" %in% names(rf_args)) {
     rf_args$write_log <- FALSE
   }
   if (!"logfile" %in% names(rf_args)) {
     rf_args$logfile <- "log_rf.txt"
   }
   if (!"threads" %in% names(rf_args)) {
     message("No number of threads provided for kNN imputation. Will use only 1 then.")
     rf_args$threads <- 1
   }
  imputed_matrix <- rfImputation(sep_gt, maxiter = rf_args$maxiter, ntree = rf_args$ntree, chunk_size = rf_args$chunk_size, write_log = rf_args$write_log, logfile = rf_args$logfile, threads = rf_args$threads)
 } else {
   stop("Please provide a valid method! ('mean', 'kNN', 'rf')")
 }
 object@imp_gt <- imputed_matrix
 return(object)
}

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
#' @export
#'
#' @examples
#' example_matrix <- matrix(c("0", "1", ".", "1", "1", ".", "0", "0", ".", "1", "1", ".", "0", "0", "0", "0", ".", "1", "0", "0", ".", "1", "1", ".", "0"), nrow = 5, byrow = TRUE)
#' rfImputation(example_matrix, maxiter = 10, ntree = 100, chunk_size = 1000)
#'
#' @export

rfImputation <- function(sep_gt, maxiter = 10, ntree = 100, chunk_size = 1000, threads = NULL, write_log = FALSE, logfile = "logfile.txt") {

  # Convert missing values ("." entries) to NA
  sep_gt[sep_gt == "."] <- NA

  # Count NAs before imputation
  na_count_before <- sum(is.na(sep_gt))

  # Convert the character matrix to a numeric matrix, coercing NA where appropriate
  genotype_matrix <- matrix(as.numeric(sep_gt), nrow = nrow(sep_gt), ncol = ncol(sep_gt))

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

  # Initialize cluster
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  # Perform imputation in parallel
  imputed_chunks <- foreach(chunk = chunks, .packages = "missForest") %dopar% {
    # Do computation with or without writing progress to log file
    if (write_log) {
      mc <- chunk$matrix
      id <- chunk$id
      # Perform missForest imputation on each chunk
      imputed_chunk <- missForest(mc, maxiter = maxiter, ntree = ntree)$ximp
      # Ensure imputed values are integers
      imputed_chunk <- matrix(as.character(round(as.numeric(imputed_chunk))),
                              ncol = ncol(imputed_chunk),
                              nrow = nrow(imputed_chunk))
      log_progress(paste0(Sys.time(), " Completed task ", id, " of ", num_chunks, "."), logfile)
      return(imputed_chunk)
    } else {
      # Perform missForest imputation on each chunk
      imputed_chunk <- missForest(chunk, maxiter = maxiter, ntree = ntree)$ximp
      # Ensure imputed values are integers
      imputed_chunk <- matrix(as.character(round(as.numeric(imputed_chunk))),
                              ncol = ncol(imputed_chunk),
                              nrow = nrow(imputed_chunk))
      return(imputed_chunk)
    }
  }

  # Stop the cluster
  stopCluster(cl)

  # Combine the chunks back into a single matrix
  imputed_matrix <- do.call(rbind, imputed_chunks)

  # Convert NAs to "."
  imputed_matrix[is.na(imputed_matrix)] <- "."

  # Display message
  message(paste0("Random Forest was able to impute ", as.character(na_count_before)), " missing genotypes.")

  # Naming columns
  colnames(imputed_matrix) <- colnames(sep_gt)

  return(imputed_matrix)
}
