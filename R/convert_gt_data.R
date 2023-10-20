#' filterBiallelicSNPs
#'
#' Filter for only biallelic SNPs in the data set.
#'
#' @param object A S4 object of class vcfR.
#'
#' @return A S4 object of the same class, but complex and multiallelic variantes are removed from the @fix and @gt slots.
#'
#' @examples
#' data("real", package = "GenoPop")
#' vcf <- filterBiallelicSNPs(real)
#'
#' @export

filterBiallelicSNPs <- function(object) {
# Extract fixed data
fixed_data <- as.data.frame(object@fix)

# Identify biallelic SNPs:
# - REF and ALT should be 1 base long (exclude indels)
# - ALT should not contain "," (exclude multi-allelic variants)
num_alleles <- sapply(strsplit(as.character(fixed_data$ALT), ","), length) + 1
biallelic_indices <- which(nchar(as.character(fixed_data$REF)) == 1 &
                             nchar(as.character(fixed_data$ALT)) == 1 &
                             num_alleles == 2)

new_gt <- object@gt[biallelic_indices,]
new_fix <- object@fix[biallelic_indices,]

message(paste0(as.character(nrow(object@gt) - nrow(new_gt)), " SNPs removed because they are not biallelic."))

# Subset VCF object to retain only biallelic SNPs
object@gt <- object@gt[biallelic_indices,]
object@fix <- object@fix[biallelic_indices,]

return(object)
}

#' calculatePloidyAndSepGT
#'
#' Calculate ploidy levels and separate the genotypes into a matrix according to the ploidy.
#'
#' @param object A S4 object of class vcfR.
#'
#' @return A S4 object of the same class but with following slots added:
#' * ploidy (integer)
#' * sep_gt (matrix)
#'
#' @examples
#' data("real", package = "GenoPop")
#' vcf <- calculatePloidyAndSepGT(real)
#' vcf@ploidy
#' head(vcf@sep_gt)
#'
#' @export

calculatePloidyAndSepGT <- function(object) {

  # Remove complex and multiallelic SNPs
  object <- filterBiallelicSNPs(object)
  gt_data <- object@gt
  # Extract genotype data and exclude the "FORMAT" column
  if ("FORMAT" %in% colnames(gt_data)) {
    gt_data <- gt_data[, -1, drop = FALSE]
  }

  # Convert the matrix to a character vector
  gt_vector <- as.character(as.vector(gt_data))

  # Check if there are multiple fields in the genotype data
  if (any(grepl(":", gt_vector))) {
    # Extract the GT part (assumes it's the first field)
    gt_vector <- sapply(strsplit(gt_vector, ":"), `[`, 1)
  }

  # Reshape the extracted GT vector back into the original matrix dimensions
  gt_data <- matrix(gt_vector, nrow = nrow(gt_data), ncol = ncol(gt_data),
                    dimnames = list(row.names(gt_data), colnames(gt_data)))


  # Detect separators for alleles (commonly '/' or '|')
  separators <- unique(gsub("[^/|]", "", gt_vector))

  # For simplicity, we'll assume one separator is used in the entire dataset
  separator <- separators[1]

  # Escape the pipe character if it's the separator
  if (separator == "|") {
    separator <- "\\|"
  }

  # Estimate ploidy level from the original genotype data
  example_gt <- gt_vector[1]  # Using the first variant as an example
  ploidy <- length(unlist(strsplit(example_gt, separator)))

  # Separate alleles into different columns
  separated_gt_data <- matrix(NA, nrow = nrow(gt_data), ncol = ploidy * ncol(gt_data))
  colnames(separated_gt_data) <- paste(rep(colnames(gt_data), each = ploidy),
                                       rep(1:ploidy, times = ncol(gt_data)),
                                       sep = "_")
  rownames(separated_gt_data) <- rownames(gt_data)

  for (i in seq_len(nrow(gt_data))) {
    # Get the genotype vector for the i-th variant
    variant_gt_vector <- gt_vector[seq(i, length(gt_vector), by = nrow(gt_data))]

    # Separate alleles and assign to the matrix
    separated_gt_data[i, ] <- unlist(strsplit(variant_gt_vector, separator))
  }

  # Create an object of the new class and populate the new slots
  new_object <- new("myVcfR", object)
  new_object@ploidy <- ploidy
  new_object@sep_gt <- separated_gt_data

  return(new_object)
}


#' seperateByPopulations
#'
#' Seperates a vcfR object into new objects per population. Needs to be done prior to calculating allele frequencies.
#'
#' @param object A S4 object of class vcfR.
#' @param pop_assignments A named vector. Elements are the population names and names are the individual name.
#'
#' @return A list containing one vcfR object per population.
#'
#' @examples
#' mys <- c("8449", "8128", "8779", "8816", "8823", "8157")
#' dav <- c("8213", "8241", "8232", "8224", "10165", "8221", "8813", "8825", "8182", "8187")
#'
#' individuals <- c(mys, dav)
#' pop_names <- c(rep("mys", length(mys)), rep("dav", length(dav)))
#' pop_assignments <- setNames(pop_names, individuals)
#'
#' data("real", package = "GenoPop")
#' vcfs <- seperateByPopulations(real, pop_assignments)
#'
#' @export

seperateByPopulations <- function(object, pop_assignments) {
  # Ensure names of pop_assignments match colnames of the vcf genotypic data
  if(!all(names(pop_assignments) %in% colnames(object@gt))) {
    stop("All individual names in pop_assignments must match those in the VCF.")
  }

  # Split individual names by population
  inds_by_pop <- split(names(pop_assignments), pop_assignments)

  # Initialize a list to store the separated vcf objects
  vcf_by_pop <- vector("list", length = length(inds_by_pop))
  names(vcf_by_pop) <- names(inds_by_pop)

  # Loop through each population, extract the individuals, and store in the list
  for (pop in names(inds_by_pop)) {
    inds <- inds_by_pop[[pop]]

    # Subset @gt slot
    gt_pop <- object@gt[, inds]
    gt_vector <- as.character(as.vector(gt_pop))
    # Check if there are multiple fields in the genotype data
    if (any(grepl(":", gt_vector))) {
      # Extract the GT part (assumes it's the first field)
      gt_vector <- sapply(strsplit(gt_vector, ":"), `[`, 1)
    }

    # Reshape the extracted GT vector back into the original matrix dimensions
    gt_values <- matrix(gt_vector, nrow = nrow(gt_pop), ncol = ncol(gt_pop),
                        dimnames = list(row.names(gt_pop), colnames(gt_pop)))

    # Identify rows with only reference alleles for all individuals
    non_ref_rows <- apply(gt_values, 1, function(x) !all(x %in% c("0/0", "0|0")))
    # Filter @gt for non-reference genotypes
    gt_pop <- gt_pop[non_ref_rows, ]
    # Common @meta slot
    meta_pop <- object@meta

    # Subset and filter @fix slot
    fix_pop <- object@fix[non_ref_rows, ]

    # Create a new vcfR object
    vcf_pop <- new("vcfR", gt = gt_pop, fix = fix_pop, meta = meta_pop)

    # Store the vcfR object in the list
    vcf_by_pop[[pop]] <- vcf_pop
  }


  return(vcf_by_pop)
}



#' calculateAlleleFreqs
#'
#' Calculate allele frequencies from a vcfR object. Will also add the slots ploidy and sep_gt, if not already present.
#'
#' @param object A S4 object of class vcfR.
#' @param missing_data Method to deal with missing data. Options are "remove", "impute", "none". Default is "none". In case of "remove", function needs the additional parameter threshold, which is the fraction of missing data in a variant that is still acceptable. In case of "impute" the function needs the additional parameters "method", with which the imputation method can be chosen. See \code{\link{imputeMissingData}}
#' @param ... Additional parameters for how to deal with missing data. For imputation see \code{\link{imputeMissingData}} and for removal see \code{\link{rmMissingData}}.
#'
#' @return A S4 object of the same class but with following slots added:
#' * allele_freqs (data frame)
#'
#'
#' @examples
#' data("real", package = "GenoPop")
#' vcf <- calculateAlleleFreqs(real, missing_data = "none")
#' vcf <- calculateAlleleFreqs(real, missing_data = "remove", threshold = 0.1)
#' vcf <- calculateAlleleFreqs(real, missing_data = "impute", method = "mean")
#' head(vcf@allele_freqs)
#'
#' @export

calculateAlleleFreqs <- function(object, missing_data = "none", ...) {

  # Check if needed slots are available, if not create them
  if (!inherits(object, "myVcfR") | !("ploidy" %in% slotNames(object)) | !("sep_gt" %in% slotNames(object))) {
    message("Calculating ploidy and separating genotype data...")
    object <- calculatePloidyAndSepGT(object)
  }

  if (missing_data == "remove") {
    # Extract additional arguments for 'remove'
    remove_args <- list(...)
    # Check if 'threshold' is provided, if not use a default value
    if (!"threshold" %in% names(remove_args)) {
      message("Threshold not provided for missing_data='remove'. Using default threshold of 0.1.")
      remove_args$threshold <- 0.1  # Default value
    }
    object <- rmMissingData(object, remove_args$threshold)
    separated_gt_data <- object@sep_gt

  } else if (missing_data == "impute") {
    # Extract additional arguments for 'impute'
    impute_args <- list(...)
    # Check if 'imputation_method' is provided, if not use a default value
    if (!"method" %in% names(impute_args)) {
      message("No valid imputation method provided for missing_data='impute'. Will use method='mean' per default.")
      impute_args$method = "mean"
    }
    # This is to calculate stats about the missing data without removing them.
    object <- rmMissingData(object, 1)
    object <- imputeMissingData(object, impute_args$method)
    separated_gt_data <- object@imp_gt
  } else if (missing_data == "none") {
    # This is to calculate stats about the missing data without removing them.
    object <- rmMissingData(object, 1)
    separated_gt_data <- object@sep_gt
  } else {
    stop("Please provide a valid way to deal with missing data! ('none', 'remove', or 'impute')")
  }

  # Calculate Allele Frequencies
  allele_frequencies_per_site <- vector("list", length = nrow(separated_gt_data))
  unique_alleles <- unique(as.vector(separated_gt_data))[unique(as.vector(separated_gt_data)) != "."]
  for (i in seq_len(nrow(separated_gt_data))) {
    alleles <- separated_gt_data[i, ]
    alleles <- alleles[alleles != "."]
    allele_counts <- table(factor(alleles, levels = unique_alleles))
    allele_frequencies <- allele_counts / sum(allele_counts)
    allele_frequencies_per_site[[i]] <- allele_frequencies
  }
  allele_frequencies_df <- do.call(rbind, allele_frequencies_per_site)

  object@allele_freqs <- allele_frequencies_df

  return(object)

}

#' calculateWindowedMetric
#'
#' Calculate one of the population genomics metrics of this package on a per window basis over a longer sequence or even whole chromsomes and genomes. Calculations are done in parallel.
#'
#' @param object An S4 object of type myVcfR.
#' @param window_size The size of the window for which Pi is calculated. (Default = 1000)
#' @param step_size The size of the step in between windows. (Default = 0)
#' @param min_var Minimum number of variants that must be present in a window to calculate the metric. Default is set to 2, because many metrics break if there is only one or none variant to work with.
#'
#' @return A data frame with four columns, the window chromosome, the window start and end postion, the number of variants in the window, and the value of the metric.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' calculateWindowedMetric(mys, Pi, window_size = 10000)
#'
#' @export

calculateWindowedMetric <- function(object, metricFunction, window_size = 1000, step_size = 0, min_var = 2) {
  # Convert the function name to a string for comparison
  metricFunctionName <- deparse(substitute(metricFunction))

  # Extract the positions and the chromosomes
  variant_positions <- as.numeric(object@fix[, 2])  # Assuming the second column contains positions
  chromosomes <- object@fix[, 1]  # Assuming the first column contains chromosome names

  # Prepare a list to hold results
  results <- list()

  # We will calculate the metric for each chromosome separately
  unique_chromosomes <- unique(chromosomes)
  for (chr in unique_chromosomes) {
    # Get indices of the variants on this chromosome
    chr_indices <- which(chromosomes == chr)
    chr_positions <- variant_positions[chr_indices]

    # Calculate the number of windows based on the step size and the positions of the variants
    num_windows <- ceiling((max(chr_positions) - min(chr_positions)) / (step_size + window_size))

    # Determine the number of cores
    num_cores <- min(detectCores() - 1, num_windows)  # Reserve one core for the system

    # Set up the parallel backend
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)

    # Prepare a list of arguments for each task
    tasks <- lapply(seq_len(num_windows), function(i) {
      start_pos <- (i - 1) * (window_size + step_size)
      end_pos <- start_pos + window_size

      # Get indices of variants within this window
      window_indices <- which(chr_positions >= start_pos & chr_positions <= end_pos)
      calculate <- TRUE
      imp <- FALSE
      # If there are less then the minimum number of variants in the window, don't calculate the metric.
      if (length(window_indices) < min_var) {
        calculate <- FALSE
      }
      if (!is.null(object@imp_gt) && nrow(object@imp_gt) != 0) {
        imp <- TRUE
      }
      # Data to be processed by the task
      list(
        calculate = calculate,
        imp = imp,
        object = object,
        chr_indices = chr_indices,
        window_indices = window_indices,
        metricFunction = metricFunction,
        metricFunctionName = metricFunctionName,
        window_size = window_size,
        start = start_pos,
        end = end_pos
      )
    })

    # Prepare the metric function to be exported.
    args_list <- list(metricFunction = metricFunction)
    # Perform the calculations in parallel
    window_results <- foreach(task = tasks, .packages = c("GenoPop")) %dopar% {
      if (!task$calculate) {
        return(list(chromosome = chr, window_start = start, window_end = end, num_variants = length(task$window_indices), metric = NA))  # No data for this window
      } else {
        # Create a 'window' object that contains only the data for this specific window
        window_object <- task$object  # This will be a deep copy of the object
        # Now, adjust the data in window_object to reflect only the current window
        window_object@fix <- task$object@fix[task$chr_indices[task$window_indices], ]
        window_object@allele_freqs <- task$object@allele_freqs[task$chr_indices[task$window_indices], ]
        window_object@sep_gt <- task$object@sep_gt[task$chr_indices[task$window_indices], ]
        if (task$imp){
          window_object@imp_gt <- task$object@imp_gt[task$chr_indices[task$window_indices], ]
        }
        start <- task$start
        end <- task$end
        current_metricFunction <- args_list$metricFunction
        # Calculate the metric
        if (task$metricFunctionName %in% c("Pi", "TajimasD")) {
          metric_value <- current_metricFunction(window_object, task$window_size)  # Also hand over the window size as sequence length
        } else {
          metric_value <- current_metricFunction(window_object)
        }

        return(list(chromosome = chr, window_start = start, window_end = end, num_variants = length(task$window_indices), metric = metric_value))

      }
    }

    # Stop the cluster
    stopCluster(cl)

    # Convert the list of lists into a data frame
    window_results_df <- do.call(rbind, lapply(window_results, as.data.frame))
    # Each chromosome's results are added to the main 'results' list
    results[[chr]] <- window_results_df
  }

  # Combine all results into a single data frame. This assumes all data frames have the same column structure.
  consolidated_results <- do.call(rbind, results)
  rownames(consolidated_results) <- 1:nrow(consolidated_results)
  colnames(consolidated_results) <- c("chromosome", "window_start", "window_end", "num_variants", metricFunctionName)
  return(consolidated_results)
}

#' calculateTwoPopWindowedMetric
#'
#' Calculates on the two populations metrics in this package for windows over a longer sequence or even whole chromosomes and genomes. Calculations are done in parallel.
#'
#' @param objects An S4 object of type myVcfR.
#' @param window_size The size of the window for which Pi is calculated. (Default = 1000)
#' @param step_size The size of the step in between windows. (Default = 0)
#' @param min_var Minimum number of variants that must be present in a window to calculate the metric. Default is set to 2, because many metrics break if there is only one or none variant to work with.
#'
#' @return A data frame with four columns, the window chromosome, the window start and end postion, the number of variants in the window, and the value of the metric.
#'
#' @examples
#' data("mys", package = "GenoPop")
#' data("dav", package = "GenoPop")
#' calculateTwoPopWindowedMetric(list(mys, dav), PrivateAlleles, window_size = 10000)
#'
#' @export

calculateTwoPopWindowedMetric <- function(objects, metricFunction, window_size = 1000, step_size = 0, min_var = 2) {
  # Convert the function name to a string for comparison
  metricFunctionName <- deparse(substitute(metricFunction))

  # Extract the positions and the chromosomes from the first object as a reference
  variant_positions <- as.numeric(objects[[1]]@fix[, 2])  # Assuming the second column contains positions
  chromosomes <- objects[[1]]@fix[, 1]  # Assuming the first column contains chromosome names

  # Prepare a list to hold results
  results <- list()

  # We will calculate the metric for each chromosome separately
  unique_chromosomes <- unique(chromosomes)
  for (chr in unique_chromosomes) {
    # Get indices of the variants on this chromosome from each object
    chr_indices_list <- lapply(objects, function(object) which(object@fix[, 1] == chr))
    chr_positions_list <- lapply(objects, function(object) as.numeric(object@fix[, 2]))

    # Calculate the number of windows based on the step size and the positions of the variants
    # This assumes that the chromosome structure is consistent across populations
    num_windows <- ceiling((max(unlist(chr_positions_list)) - min(unlist(chr_positions_list))) / (step_size + window_size))

    # Determine the number of cores
    num_cores <- min(detectCores() - 1, num_windows)  # Reserve one core for the system

    # Set up the parallel backend
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)

    # Prepare a list of arguments for each task
    tasks <- lapply(seq_len(num_windows), function(i) {
      start_pos <- (i - 1) * (window_size + step_size)
      end_pos <- start_pos + window_size
      calculate <- TRUE
      # Create a list to store whether we should calculate for each window, and the corresponding indices
      window_tasks <- list()
      for (index in seq_along(objects)) {
        # Get indices of variants within this window for each object
        window_indices <- which(chr_positions_list[[index]] >= start_pos & chr_positions_list[[index]] <= end_pos)
        # Check if there are enough variants to proceed with the calculations
        if (length(window_indices) < min_var) {
          calculate <- FALSE
        }


        window_tasks[[index]] <- list( window_indices = window_indices)
      }
      list(
        calculate = calculate,
        start = start_pos,
        end = end_pos,
        window_tasks = window_tasks,
        objects = objects,
        chr_indices_list = chr_indices_list,
        metricFunction = metricFunction
      )
    })

    # Perform the calculations in parallel
    window_results <- foreach(task = tasks, .packages = c("GenoPop")) %dopar% {
      if (!task$calculate) {
        if (identical(metricFunction, PrivateAlleles)){
          return(list(chromosome = chr, window_start = task$start, window_end = task$end, metric1 = NA, metric2 = NA))  # No data for this window
        } else {
          return(list(chromosome = chr, window_start = task$start, window_end = task$end, metric = NA))  # No data for this window
        }
      } else {
        # Create a 'window_objects' list that contains only the data for this specific window for each population
        window_objects <- lapply(seq_along(task$objects), function(index) {

          window_object <- task$objects[[index]]  # This will be a deep copy of the object

          # Adjust the data in window_object to reflect only the current window
          current_chr_indices <- task$chr_indices_list[[index]]
          current_window_indices <- task$window_tasks[[index]]$window_indices

          window_object@fix <- window_object@fix[current_chr_indices[current_window_indices], ]
          window_object@allele_freqs <- window_object@allele_freqs[current_chr_indices[current_window_indices], ]
          window_object@sep_gt <- window_object@sep_gt[current_chr_indices[current_window_indices], ]
          window_object@imp_gt <- window_object@imp_gt[current_chr_indices[current_window_indices], ]  # if available

          return(window_object)
        })

        # Calculate the metric
        metric_value <- metricFunction(window_objects)
        if (length(metric_value) == 1) {
          return(list(chromosome = chr, window_start = task$start, window_end = task$end, metric = metric_value))
        } else {
          return(list(chromosome = chr, window_start = task$start, window_end = task$end, metric1 = metric_value[1], metric2 = metric_value[2]))
        }
      }
    }
    # Stop the cluster
    stopCluster(cl)

    # Convert the list of lists into a data frame
    window_results_df <- do.call(rbind, lapply(window_results, as.data.frame))
    # Each chromosome's results are added to the main 'results' list
    results[[chr]] <- window_results_df
  }

  # Combine all results into a single data frame. This assumes all data frames have the same column structure.
  consolidated_results <- do.call(rbind, results)
  rownames(consolidated_results) <- 1:nrow(consolidated_results)
  if (identical(metricFunction, PrivateAlleles)){
    colnames(consolidated_results) <- c("chromosome", "window_start", "window_end", paste0(metricFunctionName, "_pop1"), paste0(metricFunctionName, "_pop2"))
  } else {
    colnames(consolidated_results) <- c("chromosome", "window_start", "window_end", metricFunctionName)
  }
  return(consolidated_results)
}


