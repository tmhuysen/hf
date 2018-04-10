#include "BaseSCFSolver.hpp"


namespace hf {
namespace rhf {
namespace solver {

/*
 *  CONSTRUCTORS
 */

/**
 *  Protected constructor to initialize the const S, H_core, g, calculateP, calculateG, threshold and maximum_number_of_iterations.
 */
BaseSCFSolver::BaseSCFSolver(const Eigen::MatrixXd S, const Eigen::MatrixXd H_core, const Eigen::Tensor<double ,4> g, const hf::DensityFunction calculateP,
                             const hf::TwoElectronMatrixFunction calculateG, double threshold,
                             size_t maximum_number_of_iterations) : S(S), H_core(H_core), g(g), calculateP(calculateP),
                                                                    calculateG(calculateG),threshold(threshold),
                                                                    maximum_number_of_iterations(maximum_number_of_iterations) {

}

/*
 *  GETTERS
 */

Eigen::VectorXd BaseSCFSolver::get_orbital_energies() const {

    if(this->is_converged){
        return this->orbital_energies;
    }else{
        throw std::runtime_error("The SCF procedure hasn't converged yet and you are trying to get the orbital_energies.");
    }
}

Eigen::MatrixXd BaseSCFSolver::get_C_canonical() const {

    if(this->is_converged){
        return this->C_canonical;
    }else{
        throw std::runtime_error("The SCF procedure hasn't converged yet and you are trying to get the canonical matrix.");
    }
}



} // solver
} // rhf
} // hf