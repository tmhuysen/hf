#include "SCF_functions.hpp"

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
    auto nbf = C.cols();

    // Initialize P to a zero matrix
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero (nbf, nbf);

    // This is a very naive implementation to calculate P
    for (int mu = 0; mu < nbf; ++mu) {
        for (int nu = 0; nu < nbf; ++nu) {
            for (int a = 0; a < N / 2; ++a) {
                P(mu, nu) += 2 * C(mu, a) * C.conjugate()(nu, a);
            }
        }

    }

}




}