// LSM_PGD_Cpp.cpp
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <vector>

// Note: It's better practice to not use 'using namespace' in header-like contexts.
// By being explicit with Eigen::, we avoid issues in the auto-generated RcppExports.cpp
// using namespace Rcpp; // Rcpp functions are usually called with Rcpp:: prefix anyway.
// using namespace Eigen; // This was the source of the issue.

// Sigmoid function (elementwise)
inline double sigmoid(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

// [[Rcpp::export]]
Rcpp::List LSM_PGD_Cpp(const Eigen::Map<Eigen::MatrixXd>& A, Eigen::MatrixXd Z, Eigen::VectorXd alpha,
                       double step_size_z, double step_size_alpha,
                       int niter, bool trace = false) {
  int N = A.rows();
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(N);
  std::vector<double> obj;
  
  // Iterative update loop
  for (int iter = 0; iter < niter; iter++) {
    // Compute Theta = alpha_i + alpha_j + <Z_i,Z_j>
    Eigen::MatrixXd Theta = (alpha * ones.transpose()) + (ones * alpha.transpose()) + (Z * Z.transpose());
    
    // Compute Phat = sigmoid(Theta) elementwise
    Eigen::MatrixXd Phat = Theta.unaryExpr(&sigmoid);
    
    // Compute objective:
    double tmp_obj = 0.0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        // Check for valid probabilities to avoid log(0)
        if (Phat(i, j) > 1e-9 && (1.0 - Phat(i, j)) > 1e-9) {
           tmp_obj += A(i, j) * std::log(Phat(i, j)) + (1.0 - A(i, j)) * std::log(1.0 - Phat(i, j));
        }
      }
      // Exclude the diagonal term
      if ((1.0 - Phat(i, i)) > 1e-9) {
        tmp_obj -= std::log(1.0 - Phat(i, i));
      }
    }
    tmp_obj /= 2.0;
    
    if (trace) {
      Rcpp::Rcout << "Iteration " << iter + 1 << " objective: " << tmp_obj << "\n";
    }
    obj.push_back(tmp_obj);
    
    // Update Z: Z <- Z + 2 * step_size_z * (A - Phat) * Z
    Eigen::MatrixXd update_Z = (A - Phat) * Z;
    Z = Z + 2.0 * step_size_z * update_Z;
    
    // Update alpha: alpha <- alpha + 2 * step_size_alpha * (A - Phat) * ones
    Eigen::VectorXd update_alpha = (A - Phat) * ones;
    alpha = alpha + 2.0 * step_size_alpha * update_alpha;
    
    // Center Z by subtracting column means
    Eigen::RowVectorXd col_mean = Z.colwise().mean();
    for (int i = 0; i < N; i++) {
      Z.row(i) -= col_mean;
    }
  }
  
  // Final Theta and Phat computation after the loop
  Eigen::MatrixXd Theta_final = (alpha * ones.transpose()) + (ones * alpha.transpose()) + (Z * Z.transpose());
  Eigen::MatrixXd Phat_final = Theta_final.unaryExpr(&sigmoid);
  double final_obj = 0.0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (Phat_final(i, j) > 1e-9 && (1.0 - Phat_final(i, j)) > 1e-9) {
           final_obj += A(i, j) * std::log(Phat_final(i, j)) + (1.0 - A(i, j)) * std::log(1.0 - Phat_final(i, j));
        }
      }
      if ((1.0 - Phat_final(i, i)) > 1e-9) {
        final_obj -= std::log(1.0 - Phat_final(i, i));
      }
    }
  final_obj /= 2.0;
  obj.push_back(final_obj);
  
  return Rcpp::List::create(Rcpp::Named("Z") = Z,
                            Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("Phat") = Phat_final,
                            Rcpp::Named("obj") = obj);
}
