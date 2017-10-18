#include "SCF_functions.hpp"

#include <iostream>


namespace HF {


/** Given the coefficient matrix C, and the number of electrons N, calculate the RHF density matrix P
 */
Eigen::MatrixXd calculate_P(Eigen::MatrixXd& C, unsigned N) {
    // In RHF, we should have an even number of electrons
    assert(N % 2 == 0);

    auto nbf = C.cols();

    // Construct the occupancy matrix
    Eigen::MatrixXd O = Eigen::MatrixXd::Zero (nbf, nbf);
    O.topLeftCorner(N/2, N/2) = Eigen::MatrixXd::Identity (N/2, N/2);

    // P = C O C^dagger
    return C * O * C.adjoint();
}


/** Given the density matrix P, and the two-electron integrals, calculate the G-matrix
 */
Eigen::MatrixXd calculate_G(Eigen::MatrixXd& P, Eigen::Tensor<double, 4>& tei_tensor) {
    // We will first have to convert the Eigen::MatrixXd P to an Eigen::Tensor<double, 2> P_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<double, 2>> P_tensor (P.data(), P.rows(), P.cols());

    std::cout << "P_tensor" << std::endl << P_tensor << std::endl << std::endl;

    // Specify the contraction pairs
    // To calculate G, we must perform two double contractions
    //      1. (mu nu|sigma lambda) P(lambda sigma)
    Eigen::array<Eigen::IndexPair<int>, 2> direct_contraction_pair = {Eigen::IndexPair<int>(3, 0), Eigen::IndexPair<int>(2, 1)};
    //      2. -0.5 (mu lambda|sigma nu) P(lambda sigma)
    Eigen::array<Eigen::IndexPair<int>, 2> exchange_contraction_pair = {Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1)};

    // Calculate both contractions (and incorporate prefactors)
    Eigen::Tensor<double, 2> direct_contraction = tei_tensor.contract(P_tensor, direct_contraction_pair);
    Eigen::Tensor<double, 2> exchange_contraction = -0.5 * tei_tensor.contract(P_tensor, exchange_contraction_pair);

    std::cout << "direct_contraction" << std::endl << direct_contraction << std::endl << std::endl;
    std::cout << "exchange_contraction" << std::endl << exchange_contraction << std::endl << std::endl;

    // The previous contractions are Eigen::Tensor<double 2> instances. In order to calculate the G matrix, we will convert them back into Eigen::MatrixXd instances.
    Eigen::Map<Eigen::MatrixXd> G1 (direct_contraction.data(), direct_contraction.dimension(0), direct_contraction.dimension(1));
    Eigen::Map<Eigen::MatrixXd> G2 (exchange_contraction.data(), exchange_contraction.dimension(0), exchange_contraction.dimension(1));

    // The result is the sum of both contractions, since prefactors were already incorporated
    return G1 + G2;
}


/** Calculate the RHF energy based on the density matrix P, the core Hamiltonian H_core and the Fock matrix F
 *
 */
double calculate_energy(Eigen::MatrixXd& P, Eigen::MatrixXd& H_core, Eigen::MatrixXd& F) {
    auto nbf = P.cols();

    double E = 0.0;

    // This is a very naive implementation to calculate E
    for (int mu = 0; mu < nbf; ++mu) {
        for (int nu = 0; nu < nbf; ++nu) {
            E += 0.5 * P(nu, mu) * (H_core(mu, nu) + F(mu, nu));
        }
    }

    return E;
}


} // namespace HF