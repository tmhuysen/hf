#ifndef HF_UTILITY_HPP
#define HF_UTILITY_HPP

#include <Eigen/Dense>


namespace hf {
namespace utility {


/** Check if two sets of eigenvalues are equal
 */
bool are_equal_eigenvalues(Eigen::VectorXd evals1, Eigen::VectorXd evals2, double tol);


/** Check if two eigenvectors are equal. This is the case if they are equal up to their sign.
 */
bool are_equal_eigenvectors(Eigen::VectorXd evec1, Eigen::VectorXd evec2, double tol);


/** Check if two sets of eigenvectors are equal.
 */
bool are_equal_sets_eigenvectors(Eigen::MatrixXd evecs1, Eigen::MatrixXd evecs2, double tol);


}  // namespace utility
}  // namespace hf


#endif //HF_UTILITY_HPP
