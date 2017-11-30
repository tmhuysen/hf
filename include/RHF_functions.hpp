#ifndef HF_RHF_FUNCTIONS_HPP
#define HF_RHF_FUNCTIONS_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

namespace hf {
namespace rhf {

/** Given a number of spatial orbitals K and a number of electrons N, calculated the index of the HOMO in the restricted case
 */
size_t HOMO_index(unsigned K, unsigned N);


/** Given a number of spatial orbitals K and a number of electrons N, calculated the index of the LUMO in the restricted case
 */
size_t LUMO_index(unsigned K, unsigned N);


/** Given the coefficient matrix C, and the number of electrons N, calculate the RHF density matrix P
 */
Eigen::MatrixXd calculate_P(Eigen::MatrixXd &C, unsigned N);


/** Given the density matrix P, and the two-electron integrals, calculate the G-matrix
 */
Eigen::MatrixXd calculate_G(Eigen::MatrixXd &P, Eigen::Tensor<double, 4> &tei);


/** Calculate the RHF energy based on the density matrix P, the core Hamiltonian H_core and the Fock matrix F
 *
 */
double calculate_electronic_energy(Eigen::MatrixXd &P, Eigen::MatrixXd &H_core, Eigen::MatrixXd &F);

} // namespace hf
} // namespace rhf

#endif // HF_RHF_FUNCTIONS_HPP
