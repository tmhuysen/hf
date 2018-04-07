#include "BaseSCFSolver.hpp"



namespace hf {
namespace rhf {
namespace solver {

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