# Welcome to the GenoPop Package Documentation

GenoPop is a R package designed to assist with population genomic analyses of data sets from non-model organisms or with low sequencing quality. It's created with the intention to simplify and streamline the analysis of large genomic datasets in VCF (Variant Call Format) files in a efficient manner, while handling problems of missing data. 

GenoPop can be divided into two parts. One part contains different genotype imputation methods to deal with missing data in the vcf file, and the second part contains several function to calculate commonly used population genomics metrics, like Fst, and Dxy.

This document will give an overview about GenoPops functionalities and usability. Detailed documentation of all the functions can be viewed in R and are also in the GenoPop_man.pdf document.

## Imputation Methods Overview

GenoPop includes a variety of imputation methods, which all work without a high quality reference panel to be applicable to data sets of non-model organisms. A key aspect of these methods is their approach to handling large genomic datasets. Recognizing that SNPs within a linkage block share the same evolutionary history, GenoPop employs the assumption that these SNPs exhibit more comparable patterns than those from different linkage blocks. This assumption justifies segmenting the dataset into smaller blocks for parallel processing. Essentially, GenoPop performs a batch-based imputation, where each batch contains SNPs likely to have similar characteristics due to their close proximity in the genome and therefore eventually shared linkage. This approach enhances the efficiency and accuracy of the imputation process, especially in large-scale genomic datasets.

### 1. Mean Imputation (`meanImputation`)

Mean Imputation is a straightforward approach where missing values are replaced with the mean (average) value of available data. This method can be applied in two ways:

- **By Variant:** The mean value for each variant (e.g., a specific SNP) is calculated using available data across all individuals. Missing values for that variant are then replaced with this mean.
- **By Individual:** The mean value for each individual is calculated using their available data across all variants. Missing values for that individual are replaced with their mean.

This method is fast and simple but assumes that the missing data are randomly distributed and that the mean is a reasonable estimate, which might not always be the case in genomic datasets.

### 2. k-Nearest Neighbors Imputation (`kNNImputation`)

k-Nearest Neighbors (kNN) is a more nuanced method that imputes missing data based on similarity to neighbors. In this context, 'neighbors' refer to other genomic data points that are most similar to the point with missing data. Key parameters include:

- **k (Number of Neighbors):** Determines how many neighboring data points are considered for imputing a missing value.
- **Max Iterations:** Controls how many times the algorithm iteratively refines its imputations.
- **Batch Size:** Impacts the number of variants processed at a time, affecting both accuracy and computational load.

kNN is more computationally intensive than Mean Imputation but often provides more accurate results, especially in datasets where patterns of similarity can be reliably identified.

 ### 3. Random Forest Imputation (`rfImputation`)

Random Forest Imputation uses an ensemble learning method, where multiple decision trees are used to predict missing values. Key aspects include:

- **Number of Trees (ntree):** The total count of decision trees used in the imputation process.
- **Max Iterations (maxiter):** Number of iterations for the algorithm to improve imputations.
- **Batch Size:** Similar to kNN, it affects the computational load and accuracy.

This method is especially powerful for complex datasets with intricate patterns, as it can capture non-linear relationships and interactions between variables. It's generally more computationally demanding than Mean Imputation but a little less then kNN.

---

Each of these methods has its strengths and is suitable for different types of genomic data and research needs. Mean Imputation is best for quick, general estimates, kNN for datasets where similarity patterns are strong and computational resources are not an issue, and Random Forest for capturing complex relationships in the data while using an acceptable amount of computational ressources.

## Population Genomics Metrics Overview

### General Functionality

Each function in this part of GenoPop is designed to calculate specific population genomics metrics directly from VCF (Variant Call Format) files. These functions are designed for efficiency and handle large genomic datasets by processing data in batches or windows. Most functions share a set of common parameters:

- **vcf_path**: Path to the VCF file.
- **threads**: Number of threads to use for parallel processing.
- **write_log**: Logical, indicating whether to write progress logs.
- **logfile**: Path to the log file where progress will be logged.
- **batch_size**: The number of variants to be processed in each batch.
- **window_size** (optional): Size of the window for windowed analysis in base pairs.
- **skip_size** (optional): Number of base pairs to skip between windows.
- **exclude_ind** (optional): Vector of individual IDs to exclude from the analysis.

In batch mode, the entire VCF file is processed at once to provide a general overview. In window mode, # the file is processed in sections to identify specific regions of interest. These functions typically return single metrics for batch mode or data frames detailing metrics per window.

### Metrics Overview

- **FixedSites**: Counts the number of sites fixed for the alternative allele. It helps identify regions with a complete fixation of an allele, potentially indicating selective sweeps or other evolutionary pressures.

- **SegregatingSites**: Counts the number of polymorphic sites, which are not fixed for the alternative allele. It's a measure of genetic variability within the population.

- **SingletonSites**: Counts the number of singleton sites, where a minor allele occurs only once in the sample. It can be an indicator of recent mutations.

- **PrivateAlleles**: Calculates the number of private alleles in two populations. Private alleles are present in one population but absent in another, providing insight into population differentiation.

- **ObservedHeterozygosity (Ho)**: Calculates the observed heterozygosity for each variant. It's a measure of genetic diversity within a population.

- **ExpectedHeterozygosity (He)**: Calculates the expected heterozygosity based on allele frequencies. It's a theoretical measure of how genetically diverse a population should be under random mating.

- **NucleotideDiversity (Pi)**: Measures the average number of nucleotide differences per site between two sequences. It's a key indicator of genetic diversity within a population.

- **Tajima's D**: A neutrality test comparing the number of segregating sites to the average number of nucleotide differences. It can suggest population expansion, selection, or bottlenecks.

- **Watterson's Theta**: A measure of genetic diversity within a population, based on the number of segregating sites.

- **Average Nucleotide Differences (Dxy)**: Measures the average number of nucleotide differences per site between two populations. It's a measure of genetic differentiation.

- **Fst**: The fixation index, measuring genetic differentiation between populations. It ranges from 0 (no differentiation) to 1 (complete differentiation).

- **OneDimSFS**: Calculates a one-dimensional site frequency spectrum, either folded or unfolded. It provides insights into allele frequency distributions within a population.

- **TwoDimSFS**: Calculates a two-dimensional site frequency spectrum for two populations. It's used to infer demographic history and population relationships.

Please note that this summary provides an overview of the functions and  their purposes. For complete understanding and appropriate usage, refer  to the detailed documentation of each function in GenoPop_man.pdf.

## Installation Instructions for GenoPop

### Prerequisites

Before installing GenoPop, make sure you have R installed on your system. You can download and install R from [CRAN](https://cran.r-project.org/).

### Installing GenoPop from GitHub

To install the GenoPop package directly from GitHub, you will need the `devtools` package in R. If you don't have `devtools` installed, you can install it by running the following command in R:

```R
install.packages("devtools")
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

### For the curta people (or others on HPC's)

R installations on curta are a bit tricky and annoying, so here is an easier way. (Which may also be easier in other HPC environments, if they have equivalent preinstalled modules.)

First download the package from github into your local directory using this command: 

```bash
git clone https://github.com/MiGurke/GenoPop.git
```

Then, in all of your slurm scripts using the package load a preinstalled R environment, that already includes all the dependencies for GenoPop. This is command you need to add for this on curta: 

```bash 
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
```

Then at the very beginning of each R script or command your starting on curta using GenoPop, you add these two lines (of course update the path to your download of the package):

```R
library(devtools)
devtools::load_all("/path/to/GenoPop/GenoPop.Rproj")
```

Then your also ready to go without the need to start a battle with conda or R to create your own installation ;)
