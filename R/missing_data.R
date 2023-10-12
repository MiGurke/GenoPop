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
#' vcf <- rmMissingData_int(real, 0.1)
#'
#' @export

rmMissingData_int <- function(object, threshold = 0.1) {
  # Check if needed slots are available, if not create them
  if (!inherits(object, "myVcfR") | !("ploidy" %in% slotNames(object)) | !("sep_gt" %in% slotNames(object))) {
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
  return(object)
}

setGeneric("rmMissingData", function(object, threshold = 0.1) {
  standardGeneric("rmMissingData")
})

setMethod("rmMissingData", "vcfR", function(object, threshold = 0.1) {
  rmMissingData_int(object, threshold)
})

#' imputeMissingData
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
#' vcf <- rmMissingData_int(real, 0.1)
#'
#' @export

imputeMissingData_int <- function(object, method = "mean") {
 NULL
}

setGeneric("imputeMissingData", function(object, method = "mean") {
  standardGeneric("imputeMissingData")
})

setMethod("imputeMissingData", "vcfR", function(object, method = "mean") {
  imputeMissingData_int(object, method)
})

#' meanImputation
#'
#' @param matrix Seperated genotype matrix from myvcfR object.
#'
#' @return Matrix with imputed missing data.
#'
#' @examples
#'
#' @export

meanImputation <- function(sep_gt) {
  # Convert missing values ("." entries) to NA
  sep_gt[sep_gt == "."] <- NA

  # Convert the character matrix to a numeric matrix, coercing NA where appropriate
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

  # Round to nearest integer if necessary
  imputed_matrix <- round(imputed_matrix)
}

#' knn_imputeR
#'
#' Description...
#'
#' @param data NumericMatrix, the data matrix.
#' @param k Integer, the number of neighbors.
#'
#' @return NumericMatrix, the imputed data matrix.
#' @export
#'
#' @examples
#' knn_imputeR(matrix, k = 3)
knn_imputeR <- function(data, k) {
  GenoPop::knn_impute(data, k)
}



#' kNNImputation
#'
#' @param matrix A seperated genotype matrix from a myvcfR object.
#'
#' @return A number, sum of x and y.
#' @export
#'
#' @examples
#' kNNImputation(matrix)
#'
#' @export

kNNImputation <- function(sep_gt, k = 3) {
  # Convert missing values ("." entries) to NA
  sep_gt[sep_gt == "."] <- NA

  # Convert the character matrix to a numeric matrix, coercing NA where appropriate
  genotype_matrix <- matrix(as.numeric(sep_gt), nrow = nrow(sep_gt))

  num_cores <- detectCores() - 1
  num_chunks <- num_cores
  # Splitting matrix into chunks
  chunk_size <- ceiling(nrow(genotype_matrix) / num_chunks)
  chunks <- list()

  for (i in seq_len(num_chunks)) {
    start_row <- (i - 1) * chunk_size + 1
    end_row <- min(i * chunk_size, nrow(genotype_matrix))
    chunks[[i]] <- genotype_matrix[start_row:end_row, , drop = FALSE]
  }

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  # Split matrix and perform imputation in parallel
  imputed_chunks <- foreach(chunk = chunks, .packages = "GenoPop") %dopar% {
    knn_impute(chunk, k)
  }

  # Stop the cluster
  stopCluster(cl)

  # Combine the chunks back into a single matrix
  imputed_matrix <- do.call(rbind, imputed_chunks)

  imputed_matrix[is.na(imputed_matrix)] <- "."
  colnames(imputed_matrix) <- colnames(sep_gt)
  return(imputed_matrix)
}


