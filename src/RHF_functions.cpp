#include "RHF_functions.hpp"

#include <iostream>


/** Given a number of spatial orbitals K and a number of electrons N, calculated the index of the HOMO in the restricted case
 */
size_t hf::rhf::HOMO_index(unsigned K, unsigned N) {
    if (N % 2 != 0) {
        throw std::invalid_argument("The unrestricted case is not supported.");
    }

    if (N > 2 * K) {
        throw std::invalid_argument("Cannot place that many electrons N in K spatial orbitals.");
    }

    return N / 2 - 1;  // Need to subtract 1 because computer indices start at 0
}


/** Given a number of spatial orbitals K and a number of electrons N, calculated the index of the LUMO in the restricted case
 */
size_t hf::rhf::LUMO_index(unsigned K, unsigned N) {
    if (N >= 2 * K) {
        throw std::invalid_argument("There is no LUMO for the given amount of electrons N and spatial orbitals K");
    }

    return HOMO_index(K, N) + 1;
}


/** Given the coefficient matrix C, and the number of electrons N, calculate the RHF density matrix P
 */
Eigen::MatrixXd hf::rhf::calculate_P(Eigen::MatrixXd& C, unsigned N) {
    auto nbf = C.cols();

    // Construct the occupancy matrix
    Eigen::MatrixXd O = Eigen::MatrixXd::Zero (nbf, nbf);
    O.topLeftCorner(N/2, N/2) = 2 * Eigen::MatrixXd::Identity (N/2, N/2);

    // P = C O C^dagger
    return C * O * C.adjoint();
}


/** Given the density matrix P, and the two-electron integrals, calculate the G-matrix
 */
Eigen::MatrixXd hf::rhf::calculate_G(Eigen::MatrixXd& P, Eigen::Tensor<double, 4>& tei_tensor) {
    // We will first have to convert the Eigen::MatrixXd P to an Eigen::Tensor<double, 2> P_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<double, 2>> P_tensor (P.data(), P.rows(), P.cols());

    // Specify the contraction pairs
    // To calculate G, we must perform two double contractions
    //      1. (mu nu|rho lambda) P(lambda rho)
    Eigen::array<Eigen::IndexPair<int>, 2> direct_contraction_pair = {Eigen::IndexPair<int>(3, 0), Eigen::IndexPair<int>(2, 1)};
    //      2. -0.5 (mu lambda|rho nu) P(lambda rho)
    Eigen::array<Eigen::IndexPair<int>, 2> exchange_contraction_pair = {Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1)};

    // Calculate both contractions (and incorporate prefactors)
    Eigen::Tensor<double, 2> direct_contraction = tei_tensor.contract(P_tensor, direct_contraction_pair);
    Eigen::Tensor<double, 2> exchange_contraction = -0.5 * tei_tensor.contract(P_tensor, exchange_contraction_pair);

    // The previous contractions are Eigen::Tensor<double 2> instances. In order to calculate the G matrix, we will convert them back into Eigen::MatrixXd instances.
    Eigen::Map<Eigen::MatrixXd> G1 (direct_contraction.data(), direct_contraction.dimension(0), direct_contraction.dimension(1));
    Eigen::Map<Eigen::MatrixXd> G2 (exchange_contraction.data(), exchange_contraction.dimension(0), exchange_contraction.dimension(1));

    // The result is the sum of both contractions, since prefactors were already incorporated
    return G1 + G2;
}


/** Calculate the RHF energy based on the density matrix P, the core Hamiltonian H_core and the Fock matrix F
 *
 */
double hf::rhf::calculate_electronic_energy(Eigen::MatrixXd& P, Eigen::MatrixXd& H_core, Eigen::MatrixXd& F) {
    // First, calculate the sum of H_core and F (this saves a contraction)
    Eigen::MatrixXd Z = H_core + F;

    // Convert the matrices Z and P to an Eigen::Tensor<double, 2> P_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<double, 2>> P_tensor (P.data(), P.rows(), P.cols());
    Eigen::TensorMap<Eigen::Tensor<double, 2>> Z_tensor (Z.data(), P.rows(), P.cols());

    // Specify the contraction pair
    // To calculate the electronic energy, we must perform a double contraction
    //      0.5 P(nu mu) Z(mu nu)
    Eigen::array<Eigen::IndexPair<int>, 2> contraction_pair = {Eigen::IndexPair<int>(0, 1), Eigen::IndexPair<int>(1, 0)};

    // Calculate the double contraction (with prefactor 0.5)
    Eigen::Tensor<double, 0> contraction = 0.5 * P_tensor.contract(Z_tensor, contraction_pair);

    // As the double contraction of two matrices is a scalar (a tensor of rank 0), we can access the value as (0).
    return contraction(0);
}
