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
#include "BaseSCFSolver.hpp"




namespace hf {
namespace rhf {
namespace solver {


/*
 *  CONSTRUCTORS
 */

/**
 *  Protected constructor to initialize all parameters.
 */
BaseSCFSolver::BaseSCFSolver(const Eigen::MatrixXd S, const Eigen::MatrixXd H_core, const Eigen::Tensor<double ,4> g,
                             const hf::DensityFunction calculateP, const hf::TwoElectronMatrixFunction calculateG,
                             double threshold, size_t maximum_number_of_iterations) :
    S (S),
    H_core (H_core),
    g (g),
    calculateP (calculateP),
    calculateG (calculateG),
    threshold (threshold),
    maximum_number_of_iterations(maximum_number_of_iterations)
{}



/*
 *  GETTERS
 */

Eigen::VectorXd BaseSCFSolver::get_orbital_energies() const {

    if (this->is_converged) {
        return this->orbital_energies;
    } else {
        throw std::logic_error("The SCF procedure hasn't converged yet and you are trying to get the orbital_energies.");
    }
}

Eigen::MatrixXd BaseSCFSolver::get_C_canonical() const {

    if (this->is_converged) {
        return this->C_canonical;
    } else {
        throw std::logic_error("The SCF procedure hasn't converged yet and you are trying to get the canonical matrix.");
    }
}


}  // namespace solver
}  // namespace rhf
}  // namespace hf
