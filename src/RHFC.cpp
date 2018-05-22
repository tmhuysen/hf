#include "RHFC.hpp"

namespace hf {
namespace rhf {


/*
 *  PRIVATE METHODS
 */

/**
 *  Given the RHF AO density matrix @param: P calculate the GA (mulliken related matrix)
 */
Eigen::MatrixXd RHFC::calculateGA(std::vector<size_t> AO_set) const {
    Eigen::MatrixXd GA = Eigen::MatrixXd::Zero(this->K,this->K);
    for(int i=0;i<this->K;i++){
        for(size_t p : AO_set){
            GA(i,p) += 0.5 * this->ao_basis.get_S()(i,p);
            GA(p,i) += 0.5 * this->ao_basis.get_S()(p,i);
        }
    }
    return GA;


}


/**
 *  Given the RHF AO density matrix @param: P calculate the mulliken population for a set of AO's
 */
double RHFC::calculatePopulation(const Eigen::MatrixXd &P, std::vector<size_t> AO_set) const {
    Eigen::MatrixXd PS = P*this->ao_basis.get_S();
    double population = 0;
    for(size_t p : AO_set){
        population += PS(p,p);
    }
    return population;
}



/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given libwint::AOBasis @param: ao_basis, a number of electrons @param: N and an SCF-cycle @param: scf_threshold, MAX_CYCLES
 */
RHFC::RHFC(const libwint::Molecule &molecule, const libwint::AOBasis &ao_basis, double scf_threshold, size_t MAX_CYCLES)
        : RHF(molecule, ao_basis, scf_threshold, MAX_CYCLES) {

}



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the restricted Hartree-Fock equations (i.e. the Roothaan-Hall equations) (with a mulliken modifier)
 */
void RHFC::solve(hf::rhf::solver::SCFSolverType solver_type, std::vector<size_t> AO_set, double multiplier) {
    // Calculate H_core
    Eigen::MatrixXd H_core = this->ao_basis.get_T() + this->ao_basis.get_V();
    Eigen::MatrixXd mod_core = H_core - multiplier*this->calculateGA(AO_set);
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes0 (mod_core,this->ao_basis.get_S());
    this->C_guess = gsaes0.eigenvectors();
    hf::DensityFunction calculateP = [this] (const Eigen::MatrixXd& x) { return this->calculateP(x);};
    hf::TwoElectronMatrixFunction calculateG = [this] (const Eigen::MatrixXd & x, const Eigen::Tensor<double, 4> & y) { return this->calculateG(x,y);};

    switch (solver_type) {

        case hf::rhf::solver::SCFSolverType::PLAIN: {
            auto plain_solver = new hf::rhf::solver::PlainSCFSolver(this->ao_basis.get_S(), mod_core, this->ao_basis.get_g(),
                                                                    calculateP, calculateG, this->scf_threshold,
                                                                    this->MAX_NUMBER_OF_SCF_CYCLES);
            plain_solver->solve(C_guess);
            this->SCF_solver_ptr = plain_solver;  // prevent data from going out of scope
            // we are only assigning this->eigensolver_ptr now, because
            // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }

        case hf::rhf::solver::SCFSolverType::DIIS: {
            auto DIIS_solver = new hf::rhf::solver::DIISSCFSolver(this->ao_basis.get_S(), mod_core, this->ao_basis.get_g(),
                                                                  calculateP, calculateG, this->scf_threshold,
                                                                  this->MAX_NUMBER_OF_SCF_CYCLES);
            DIIS_solver->solve(C_guess);
            this->SCF_solver_ptr = DIIS_solver;  // prevent data from going out of scope
            // we are only assigning this->eigensolver_ptr now, because
            // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }

    }
    this->is_converged = true;
    Eigen::MatrixXd P = calculateP(get_C_canonical());
    this->electronic_energy=calculateElectronicEnergy(P,H_core,H_core+calculateG(P,this->ao_basis.get_g()));
    this->population_set = calculatePopulation(P, AO_set);
}

void RHFC::solve(hf::rhf::solver::SCFSolverType solver_type, std::vector<size_t> AO_set, double multiplier, Eigen::MatrixXd C_guess) {
    // Calculate H_core
    Eigen::MatrixXd H_core = this->ao_basis.get_T() + this->ao_basis.get_V();
    Eigen::MatrixXd mod_core = H_core - multiplier*this->calculateGA(AO_set);
    this->C_guess = C_guess;
    hf::DensityFunction calculateP = [this] (const Eigen::MatrixXd& x) { return this->calculateP(x);};
    hf::TwoElectronMatrixFunction calculateG = [this] (const Eigen::MatrixXd & x, const Eigen::Tensor<double, 4> & y) { return this->calculateG(x,y);};
    switch (solver_type) {

        case hf::rhf::solver::SCFSolverType::PLAIN: {
            auto plain_solver = new hf::rhf::solver::PlainSCFSolver(this->ao_basis.get_S(), mod_core, this->ao_basis.get_g(),
                                                                    calculateP, calculateG, this->scf_threshold,
                                                                    this->MAX_NUMBER_OF_SCF_CYCLES);
            plain_solver->solve(C_guess);
            this->SCF_solver_ptr = plain_solver;  // prevent data from going out of scope
            // we are only assigning this->eigensolver_ptr now, because
            // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }

        case hf::rhf::solver::SCFSolverType::DIIS: {
            auto DIIS_solver = new hf::rhf::solver::DIISSCFSolver(this->ao_basis.get_S(), mod_core, this->ao_basis.get_g(),
                                                                  calculateP, calculateG, this->scf_threshold,
                                                                  this->MAX_NUMBER_OF_SCF_CYCLES);
            DIIS_solver->solve(C_guess);
            this->SCF_solver_ptr = DIIS_solver;  // prevent data from going out of scope
            // we are only assigning this->eigensolver_ptr now, because
            // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }

    }
    this->is_converged = true;
    Eigen::MatrixXd P = calculateP(get_C_canonical());
    this->electronic_energy=calculateElectronicEnergy(P,H_core,H_core+calculateG(P,this->ao_basis.get_g()));
    this->population_set = calculatePopulation(P, AO_set);
}



}  // namespace rhf
}  // namespace hf