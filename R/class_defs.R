#' @importClassesFrom vcfR vcfR
#' @importFrom methods new slotNames
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster




setClass("GPvcfR",
        contains = "vcfR",
        slots = c(
          ploidy = "numeric",
          sep_gt = "matrix",
          imp_gt = "matrix",
          allele_freqs = "matrix",
          missing_data = "list"
        )
)

setMethod("show", "GPvcfR", function(object) {
  cat("***** Object of Class GPvcfR *****\n")
  # Access and print information from inherited vcfR slots
  cat(paste(ncol(object@gt), "samples\n"))
  cat(paste(length(unique(object@fix[,1])), "CHROMs\n"))
  cat(paste(nrow(object@fix), "variants\n"))
  # Calculate and print the object size
  size_in_bytes <- object.size(object)
  size_in_mb <- size_in_bytes / (1024^2)  # Convert bytes to megabytes
  cat(sprintf("Object size: %.2f Mb\n", size_in_mb))

  # Access and print information from new slots in GPvcfR
  if(length(object@ploidy > 0)) {
    cat(paste("Ploidy:", object@ploidy, "\n"))
  }
  if (length(object@missing_data)) {
    cat(paste("Average fraction of missing data per individual:", round(mean(object@missing_data$fraction_per_individual), digits = 3)), "\n")
    cat(paste("Average fraction of missing data per variant:", round(mean(object@missing_data$fraction_per_variant), digits = 3)), "\n")
  }

  # ... add more as needed ...

  cat("***** ***** ***** *****\n")
})
