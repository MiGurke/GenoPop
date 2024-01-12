#include <Rcpp.h>
#include <cmath>  // for std::sqrt() and std::round()
using namespace Rcpp;
#include "../annoy/src/annoylib.h"
#include "../annoy/src/kissrandom.h"
using namespace Annoy;

// Function to calculate distance considering missing data
double calculate_distance(const NumericVector& instance1, const NumericVector& instance2) {
  double sum = 0.0;
  int valid_dims = 0;

  for(size_t i = 0; i < instance1.size(); ++i) {
    if(!R_IsNA(instance1[i]) && !R_IsNA(instance2[i])) {
      sum += (instance1[i] - instance2[i]) * (instance1[i] - instance2[i]);
      ++valid_dims;
    }
  }

  if(valid_dims == 0) {
    return std::nan("");
  }

  return std::sqrt(sum) / valid_dims;
}


std::vector<int> find_k_neighbors(Rcpp::NumericMatrix data, int target_row, int k) {
  int f = data.ncol();  // Number of features
  AnnoyIndex<int, double, Euclidean, Kiss32Random, AnnoyIndexSingleThreadedBuildPolicy> index(f);


  // Compute column means to replace NAs with as a temporary placeholder
  // This is to deal with data sets that have large amounts of missing data
  NumericVector means(f);
  for(int j = 0; j < f; ++j) {
    NumericVector col = data(_, j);
    double sum = 0.0;
    int count = 0;
    for(int i = 0; i < col.size(); ++i) {
      if(!R_IsNA(col[i])) {
        sum += col[i];
        ++count;
      }
    }
    means[j] = count > 0 ? sum / count : NA_REAL;
  }

  // Add all items to the index, replacing NAs with column means
  for(int i = 0; i < data.nrow(); ++i) {
    std::vector<double> item;
    for(int j = 0; j < f; ++j) {
      item.push_back(R_IsNA(data(i, j)) ? means[j] : data(i, j));
    }
    index.add_item(i, &item[0]);
  }



  // Build the index with 100 trees
  index.build(100);

  // Find the k + 1-nearest neighbors
  std::vector<int> neighbors;
  index.get_nns_by_item(target_row, k + 1, -1, &neighbors, nullptr);
  // Remove the first Neighbour, cause it'll be the target row itself
  neighbors.erase(neighbors.begin());

  // Print out the neighbors found (DEBUGGING)
  // Rcpp::Rcout << "Neighbors of item " << target_row << ": ";
  // for(int neighbor : neighbors) {
  //   Rcpp::Rcout << neighbor << " ";
  // }
  // Rcpp::Rcout << std::endl;

  return neighbors;
}

// [[Rcpp::export]]
NumericMatrix knn_impute(NumericMatrix data, int k, int maxiter) {
  int num_rows = data.nrow();
  int num_cols = data.ncol();

  NumericMatrix imputed_data(clone(data));
  for(int m = 0; m < maxiter; ++m) {
    bool any_na = false; // Flag to check if any NA values are left

    for(int i = 0; i < num_rows; ++i) {
      bool row_has_na = false;
      for(int j = 0; j < num_cols; ++j) {
        if(R_IsNA(imputed_data(i, j))) {
          row_has_na = true;
          any_na = true;
          break;
        }
      }

      if(row_has_na) {
        std::vector<int> neighbors = find_k_neighbors(imputed_data, i, k);

        for(int j = 0; j < num_cols; ++j) {
          if(R_IsNA(imputed_data(i, j))) {
            double weighted_sum = 0.0;
            double weight_total = 0.0;
            for(int neighbor_idx : neighbors) {
              if(!R_IsNA(imputed_data(neighbor_idx, j))) {
                double weight = 1.0 / (calculate_distance(imputed_data.row(i), imputed_data.row(neighbor_idx)) + 1e-8);
                weighted_sum += weight * imputed_data(neighbor_idx, j);
                weight_total += weight;
              }
            }
            // Only impute if the neighbors actually had any data for this individual/allele/column.
            // (Weight won't be increased from neighbours that have missing data at the position. See above.)
            // If not its kinda pointless to impute this. This is the reason, why  this algorithm may not impute all missing values.
            // Especially if there are individuals with high proportions of missing data.
            if(weight_total > 0) {
              imputed_data(i, j) = std::round(weighted_sum / weight_total);
            }
          }
        }
      }
    }

    // Break out of the loop if no NA values are left
    if(!any_na) {
      break;
    }
  }

  return imputed_data;
}

