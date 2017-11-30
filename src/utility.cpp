#include "utility.hpp"


/** Check if two sets of eigenvalues are equal
 */
bool hf::utility::are_equal_eigenvalues(Eigen::VectorXd evals1, Eigen::VectorXd evals2, double tol) {
    return evals1.isApprox(evals2, tol);
}

/** Check if two eigenvectors are equal. This is the case if they are equal up to their sign.
 */
bool hf::utility::are_equal_eigenvectors(Eigen::VectorXd evec1, Eigen::VectorXd evec2, double tol) {
    return (evec1.isApprox(evec2, tol) || evec1.isApprox(-evec2, tol));
}

/** Check if two sets of eigenvectors are equal.
 */
bool hf::utility::are_equal_sets_eigenvectors(Eigen::MatrixXd evecs1, Eigen::MatrixXd evecs2, double tol) {
    auto dim = evecs1.cols();
    for (unsigned i = 0; i < dim; i++) {
        if (!are_equal_eigenvectors(evecs1.col(i), evecs2.col(i), tol)) {
            return false;
        }
    }
    return true;
}
