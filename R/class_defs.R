#' @importClassesFrom vcfR vcfR
#' @importFrom methods new slotNames


setClass("myVcfR",
        contains = "vcfR",
        slots = c(
          ploidy = "numeric",
          sep_gt = "matrix",
          allele_freqs = "matrix",
          missing_data = "list"
        )
)
