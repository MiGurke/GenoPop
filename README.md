# GenoPop

A R package to perform several population genomics analyses directly on whole genome vcf files. 

Most important features are:
* Reading genotype data ready to use from vcf files using the vcfR package.
* Different methods to impute or remove missing data from the genotype matrix.
* Calculation of several commonly used population genomics metrics from the genotype data. 
* Window based analyses with those different metrics.
* Parallelization and optimization of heavy tasks to enable processing of whole genomes (imputation and window based analyses). 

## Installation Instructions for GenoPop

### Prerequisites

Before installing GenoPop, make sure you have R installed on your system. You can download and install R from [CRAN](https://cran.r-project.org/).

### Installing GenoPop from GitHub

To install the GenoPop package directly from GitHub, you will need the `devtools` package in R. If you don't have `devtools` installed, you can install it by running the following command in R:

```R
install.packages("devtools")
```
On clusters like curta or mfn it is easiest to create a new conda environment and then install devtools via conda with this command: 

```
conda install -c conda-forge r-devtools
```

Once devtools is installed, you can install GenoPop using the install_github function. Run the following commands in your R console:

```R
library(devtools)
install_github("https://github.com/MiGurke/GenoPop")
```
Then load the package..
```R
library(GenoPop)
```
and your ready to go!

### Dependencies

For proper compression of vcf's newly generated, you need to have tabix installed on your machine. All other dependencies are supposed to be handled by R itself. If that is causing any problems, these are the other R packages GenoPop depends on:
  * vcfR
  * methods
  * Rcpp
  * foreach
  * doParallel
  * parallel
  * missForest
They can also installed via conda, if needed.

## Getting started

To get your vcf formated data into R use the reading function of the vcfR package: 

```R
vcf <- read.vcfR( "example.vcf", verbose = FALSE )
```

It can read compressed and umcompressed vcf files. 

From there a good point to start your analysis is to use the calculateAlleleFreqs function to just populate the most important data and information slots from your vcf file automatically. 

```R
example <- calculateAlleleFreqs(vcf)
```

This will give you a nicely formatted genotype matrix (example@sep_gt), calculated allele frequencies (example@allele_freqs) and some interesting stats about the amount of missing data in you data set (example@missing_data).

From there it depends on you how to continue. You can fix some issue with missing data by imputing or removing it, directly calculate some stats like Fst, Pi, and the site frequency spectrum, or even do a window based analysis. Just have a look at the available functions in the man pages.

If you need further assistance or got suggestions for this package, feel free to open an issue on the GenoPop GitHub repo or contact me in any other way.



