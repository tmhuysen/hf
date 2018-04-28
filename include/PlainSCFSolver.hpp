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
#ifndef HF_PLAINSCFSOLVER_HPP
#define HF_PLAINSCFSOLVER_HPP

#include "BaseSCFSolver.hpp"



namespace hf {
namespace rhf {
namespace solver {


class PlainSCFSolver: public BaseSCFSolver {
public:
    // CONSTRUCTOR
    /**
     *  Constructor based on a given (AO) overlap matrix @param S, one-electron integrals @param H_core, two-electron
     *  integrals @param g, the function @param calculateP, the function calculateG, an SCF convergence threshold
     *  @param threshold and a @param maximum_number_of_iterations.
     */
    PlainSCFSolver(const Eigen::MatrixXd S, const Eigen::MatrixXd H_core, const Eigen::Tensor<double ,4> g,
                   const hf::DensityFunction calculateP, const hf::TwoElectronMatrixFunction calculateG,
                   double threshold = 1.0e-6, size_t maximum_number_of_iterations = 128);

    // DESTRUCTOR
    ~PlainSCFSolver() override = default;


    // PUBLIC PURE VIRTUAL METHODS
    /**
     *  Execute the SCF procedure.
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



#endif  // HF_PLAINSCFSOLVER_HPP
