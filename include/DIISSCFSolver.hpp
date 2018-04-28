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
#ifndef HF_DIISSCFSOLVER_HPP
#define HF_DIISSCFSOLVER_HPP


#include <deque>
#include "BaseSCFSolver.hpp"



namespace hf {
namespace rhf {
namespace solver {


class DIISSCFSolver: public BaseSCFSolver {
private:
    const size_t maximum_subspace_dimension = 6;  // maximum subspace size before collapsing

    std::deque<Eigen::MatrixXd> fock_matrix_deque;  // deque of fock matrices used in the DIIS algorithm
    std::deque<Eigen::MatrixXd> error_matrix_deque;  // deque of error matrices used in the DIIS algorithm



public:
    // CONSTRUCTOR
    /**
     *  Constructor based on a given (AO) overlap matrix @param S, one-electron integrals @param H_core, two-electron
     *  integrals @param g, the function @param calculateP, the function calculateG, an SCF convergence threshold
     *  @param threshold and a @param maximum_number_of_iterations.
     */
    explicit DIISSCFSolver(const Eigen::MatrixXd S, const Eigen::MatrixXd H_core, const Eigen::Tensor<double ,4> g,
                            const hf::DensityFunction calculateP,
                            const hf::TwoElectronMatrixFunction calculateG,
                            double threshold = 1.0e-6, size_t maximum_number_of_iterations = 128);


    // DESTRUCTOR
    ~DIISSCFSolver() override = default;


    // PUBLIC METHODS
    /**
     *  Execute the SCF procedure using DIIS (direct inversion of the iterative subspace).
     *
     *  If successful, it sets
     *      - @member is_converged to true
     *      - @member C_canonical
     *      - @member orbital_energies
     */
    void solve() override;
};


}  // namespace solver
}  // namespace rhf
}  // namespace hf


#endif  // HF_DIISSCFSOLVER_HPP
