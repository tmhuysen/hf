#ifndef HF_SCF_FUNCTIONS_HPP
#define HF_SCF_FUNCTIONS_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


namespace HF {

/** Given the coefficient matrix C, and the number of electrons N, calculate the RHF density matrix P
 */
Eigen::MatrixXd calculate_P(Eigen::MatrixXd& C, unsigned N);


/** Given the density matrix P, and the two-electron integrals, calculate the G-matrix
 */
Eigen::MatrixXd calculate_G(Eigen::MatrixXd& P, Eigen::Tensor<double, 4>& tei);


/** Calculate the RHF energy based on the density matrix P, the core Hamiltonian H_core and the Fock matrix F
 *
 */
double calculate_electronic_energy(Eigen::MatrixXd& P, Eigen::MatrixXd& H_core, Eigen::MatrixXd& F);

} // namespace HF

#endif //HF_SCF_FUNCTIONS_HPP
