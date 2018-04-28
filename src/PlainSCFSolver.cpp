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
#include "PlainSCFSolver.hpp"



namespace hf {
namespace rhf {
namespace solver {


/*
 * CONSTRUCTORS
 */

/**
 *  Constructor based on a given (AO) overlap matrix @param S, one-electron integrals @param H_core, two-electron
 *  integrals @param g, the function @param calculateP, the function calculateG, an SCF convergence threshold
 *  @param threshold and a @param maximum_number_of_iterations.
 */
PlainSCFSolver::PlainSCFSolver(const Eigen::MatrixXd S, const Eigen::MatrixXd H_core, const Eigen::Tensor<double ,4> g,
                               const hf::DensityFunction calculateP, const hf::TwoElectronMatrixFunction calculateG,
                               double threshold, size_t maximum_number_of_iterations
                               ) :
    BaseSCFSolver(S, H_core, g, calculateP, calculateG, threshold, maximum_number_of_iterations)
{}



/*
 * PUBLIC METHODS
 */

/**
 *  Execute the SCF procedure.
 *
 *  If successful, it sets
 *      - @member is_converged to true
 *      - @member C_canonical
 *      - @member orbital_energies
 */
void PlainSCFSolver::solve() {
    // Solve the generalized eigenvalue problem for H_core to obtain an initial guess for the density matrix P
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

            std::cout<<std::endl<<"SCF ITERATIONS : "<<iteration_counter<<std::endl;
        } else {  // not converged yet
            iteration_counter ++;

            // If we reach more than this->maximum_number_of_iterations, the system is considered not to be converging
            if (iteration_counter >= this->maximum_number_of_iterations) {
                throw std::runtime_error("The SCF procedure did not converge.");
            }
        }
    }  // SCF cycle loop
}



} // solver
} // rhf
} // hf