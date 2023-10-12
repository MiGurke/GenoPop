#include <Rcpp.h>
#include <cmath>  // for std::sqrt() and std::round()
using namespace Rcpp;

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

std::vector<int> find_k_neighbors(NumericMatrix data, int target_row, int k, double max_missing = 0.10) {
  int num_rows = data.nrow();
  int num_cols = data.ncol();

  std::vector<std::pair<double, int>> distances;

  for(int i = 0; i < num_rows; ++i) {
    int missing_count = 0;

    for(int j = 0; j < num_cols; ++j) {
      missing_count += R_IsNA(data(i, j));
    }
    double missing_proportion = (double)missing_count / num_cols;

    if(i != target_row && missing_proportion <= max_missing) {
      double dist = calculate_distance(data.row(target_row), data.row(i));
      if(!std::isnan(dist)) {
        distances.push_back({dist, i});
      }
    }
  }

  std::sort(distances.begin(), distances.end());

  std::vector<int> neighbors;
  for(int i = 0; i < std::min(k, (int)distances.size()); ++i) {
    neighbors.push_back(distances[i].second);
  }
  return neighbors;
}

// [[Rcpp::export]]
NumericMatrix knn_impute(NumericMatrix data, int k) {
  int num_rows = data.nrow();
  int num_cols = data.ncol();

  NumericMatrix imputed_data(clone(data));

  for(int i = 0; i < num_rows; ++i) {
    for(int j = 0; j < num_cols; ++j) {
      if(R_IsNA(data(i, j))) {
        std::vector<int> neighbors = find_k_neighbors(data, i, k);

        double weighted_sum = 0.0;
        double weight_total = 0.0;
        for(int neighbor_idx : neighbors) {
          if(!R_IsNA(data(neighbor_idx, j))) {
            double weight = 1.0 / (calculate_distance(data.row(i), data.row(neighbor_idx)) + 1e-8);
            weighted_sum += weight * data(neighbor_idx, j);
            weight_total += weight;
          }
        }

        if(weight_total > 0) {
          imputed_data(i, j) = std::round(weighted_sum / weight_total);
        } else {
          Rcpp::Rcout << "Warning: No suitable neighbors found for imputation at (" << i << ", " << j << ").\n";
        }
      }
    }
  }
  return imputed_data;
}

