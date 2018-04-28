// This file is part of GQCG-hf.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-hf is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-hf is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-hf.  If not, see <http://www.gnu.org/licenses/>.
#include "DIISSCFSolver.hpp"


namespace hf {
namespace rhf {
namespace solver {


/*
 * CONSTRUCTORS
 */

/**
 *  Constructor to initialize the const S, H_core, g, calculateP, calculateG, threshold and maximum_number_of_iterations.
 */
DIISSCFSolver::DIISSCFSolver(const Eigen::MatrixXd S, const Eigen::MatrixXd H_core, const Eigen::Tensor<double, 4> g,
                             const hf::DensityFunction calculateP, const hf::TwoElectronMatrixFunction calculateG,
                             double threshold, size_t maximum_number_of_iterations, size_t maximum_subspace_dimension) :
    maximum_subspace_dimension(maximum_subspace_dimension),
    BaseSCFSolver(S, H_core, g, calculateP, calculateG, threshold, maximum_number_of_iterations)
{}



/*
 * PUBLIC METHODS
 */

/**
 *  Execute the SCF procedure using DIIS (direct inversion of the iterative subspace).
 *
 *  If successful, it sets
 *      - @member is_converged to true
 *      - @member C_canonical
 *      - @member orbital_energies
 */
void DIISSCFSolver::solve() {

    // Initialize error_matrix and fock_matrix deques for the DIIS procedure
    this->fock_matrix_deque = {};
    this->error_matrix_deque = {};


    // Solve the generalized eigenvalue problem for H_core to obtain a guess for the density matrix P
    //  H_core should be self-adjoint
    //  S should be positive definite
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes0 (this->H_core,this->S);
    Eigen::MatrixXd C = gsaes0.eigenvectors();
    Eigen::MatrixXd P = this->calculateP(C);


    size_t iteration_counter = 1;
    while (!this->is_converged) {
        // Calculate the G-matrix
        Eigen::MatrixXd G = this->calculateG(P, this->g);
        // Calculate the Fock matrix
        Eigen::MatrixXd f_AO = this->H_core + G;


        // Update deques for the DIIS procedure
        this->fock_matrix_deque.emplace_back(f_AO);

        Eigen::MatrixXd error_matrix = f_AO * P * this->S - this->S * P * f_AO;
        this->error_matrix_deque.emplace_back(error_matrix);


        // Collapse the subspace, if needed
        size_t n = error_matrix_deque.size();  // n is the current subspace dimension
        if (n == this->maximum_subspace_dimension) {

            // Initialize the augmented B matrix
            Eigen::MatrixXd B = -1 * Eigen::MatrixXd::Ones(n+1,n+1);  // +1 for the multiplier
            B(n,n) = 0;

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    // B(i,j) = Tr(e_i^T e_j)
                    B(i,j) = (this->error_matrix_deque[i].transpose() * this->error_matrix_deque[j]).trace();
                }
            }

            // Initialize the RHS of the system of equations
            Eigen::VectorXd b = Eigen::VectorXd::Zero(n+1);  // +1 for the multiplier
            b(n) = -1;  // the last entry of b is accessed through n: dimension of b is n+1 - 1 because of computers


            // Solve the DIIS non-linear equations
            Eigen::VectorXd y = B.inverse() * b;


            // Use the coefficients that are in y to construct 'the best' Fock matrix as a linear combination of previous Fock matrices
            f_AO = Eigen::MatrixXd::Zero(this->S.cols(),this->S.cols());
            for (size_t i = 0; i < n; i++) {  // n is the dimension of the subspace (not equal to the size of the augmented B matrix)
                f_AO += y(i) * this->fock_matrix_deque[i];
            }

            // Remove the oldest entries, which means that we collapse every iteration once the dimension is large enough
            this->fock_matrix_deque.pop_front();
            this->error_matrix_deque.pop_front();
        }  // DIIS


        // Solve the Roothaan equation (generalized eigenvalue problem)
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes (f_AO, this->S);
        C = gsaes.eigenvectors();

        // Calculate an improved density matrix P from the improved coefficient matrix C
        Eigen::MatrixXd P_previous = P; // We will store the previous density matrix
        P = this->calculateP(C);

        // Check for convergence on the density matrix P
        if ((P - P_previous).norm() <= this->threshold) {
            this->is_converged = true;

            // After the SCF procedure, we end up with canonical spatial orbitals, i.e. the Fock matrix should be diagonal in this basis
            // Let's check if this is the case, within double float precision
            Eigen::MatrixXd f_SO = libwint::transformations::transform_AO_to_SO(f_AO, C);
            assert(f_SO.isDiagonal());

            // Furthermore, add the orbital energies and the coefficient matrix to (this)
            this->orbital_energies = gsaes.eigenvalues();
            this->C_canonical = C;
        }

        else {
            iteration_counter ++;
            // If we reach more than this->maximum_number_of_iterations, the system is considered not to be converging
            if (iteration_counter >= this->maximum_number_of_iterations) {
                throw std::runtime_error("The SCF procedure did not converge.");
            }
        }
    }  // SCF cycle loop
}


}  // namespace solver
}  // namespace rhf
}  // namespace hf
