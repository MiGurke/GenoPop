#' @importClassesFrom vcfR vcfR
#' @importFrom methods new slotNames
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster




setClass("myVcfR",
        contains = "vcfR",
        slots = c(
          ploidy = "numeric",
          sep_gt = "matrix",
          imp_gt = "matrix",
          allele_freqs = "matrix",
          missing_data = "list"
        )
)
