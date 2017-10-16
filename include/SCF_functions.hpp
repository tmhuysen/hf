#ifndef HARTREE_FOCK_SCF_FUNCTIONS_HPP
#define HARTREE_FOCK_SCF_FUNCTIONS_HPP

#include <Eigen/Dense>

namespace HF {

/** Given S, calculate the matrix X=S^-1/2
 */
Eigen::MatrixXd calculate_X(Eigen::MatrixXd& S);

/** Given the coefficient matrix C, and the number of electrons N, calculate the RHF density matrix P
 */
Eigen::MatrixXd calculate_P(Eigen::MatrixXd& C, unsigned N);

}

#endif //HARTREE_FOCK_SCF_FUNCTIONS_HPP
