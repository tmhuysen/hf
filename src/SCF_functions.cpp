#include "SCF_functions.hpp"

#include <iostream>


namespace HF {

/** Given S, calculate the matrix X=S^-1/2
 */
Eigen::MatrixXd calculate_X(Eigen::MatrixXd& S) {
    // Diagonalize S
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (S);
    Eigen::MatrixXd U = saes.eigenvectors();

    // s_invsqrt is the diagonal matrix consisting of the inverse square roots of the eigenvalues of S
    Eigen::MatrixXd s_invsqrt = saes.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal();

    // X = U * s_invsqrt * U^dagger
    return U * s_invsqrt * U.adjoint();
}


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
Eigen::MatrixXd calculate_G(Eigen::MatrixXd& P, Eigen::Tensor<double, 4>& tei) {
    auto nbf = P.cols();

    // Initialize G to a zero matrix
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero (nbf, nbf);

    // This is a very naive implementation to calculate G
    for (int mu = 0; mu < nbf; ++mu) {
        for (int nu = 0; nu < nbf; ++nu) {
            for (int lambda = 0; lambda < nbf; ++lambda) {
                for (int sigma = 0; sigma < nbf ; ++sigma) {
                    G(mu, nu) += P(lambda, sigma) * (tei(mu, lambda, nu, sigma) - 0.5 * tei(mu, sigma, lambda, nu));
                }
            }
        }
    }

    return G;
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