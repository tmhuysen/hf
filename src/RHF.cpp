#include "RHF.hpp"


namespace hf {
namespace rhf {


/*
 *  PRIVATE METHODS
 */

/**
 *  Given a coefficient matrix @param: C, and the number of electrons N, calculate the RHF AO density matrix P
 */
Eigen::MatrixXd RHF::calculateP(const Eigen::MatrixXd& C) const {

    // Construct the occupancy matrix O
    Eigen::MatrixXd O = Eigen::MatrixXd::Zero (this->K, this->K);
    O.topLeftCorner(this->N/2, this->N/2) = 2 * Eigen::MatrixXd::Identity(this->N/2, this->N/2);


    // P = C * O * C^dagger
    return C * O * C.adjoint();
}


/**
 *  Given the RHF AO density matrix @param: P, and the two-electron integrals @param g in chemist's notation, calculate the RHF G-matrix
 */
Eigen::MatrixXd RHF::calculateG(const Eigen::MatrixXd& P, const Eigen::Tensor<double, 4>& g) const {

    // We will first have to convert the Eigen::MatrixXd P to an Eigen::Tensor<const double, 2> P_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> P_tensor (P.data(), P.rows(), P.cols());

    // Specify the contraction pairs
    // To calculate G, we must perform two double contractions
    //      1. (mu nu|rho lambda) P(lambda rho)
    Eigen::array<Eigen::IndexPair<int>, 2> direct_contraction_pair = {Eigen::IndexPair<int>(3, 0), Eigen::IndexPair<int>(2, 1)};
    //      2. -0.5 (mu lambda|rho nu) P(lambda rho)
    Eigen::array<Eigen::IndexPair<int>, 2> exchange_contraction_pair = {Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1)};

    // Calculate both contractions (and incorporate prefactors)
    Eigen::Tensor<double, 2> direct_contraction = g.contract(P_tensor, direct_contraction_pair);
    Eigen::Tensor<double, 2> exchange_contraction = -0.5 * g.contract(P_tensor, exchange_contraction_pair);

    // The previous contractions are Eigen::Tensor<double 2> instances. In order to calculate the G matrix, we will convert them back into Eigen::MatrixXd instances.
    Eigen::Map<Eigen::MatrixXd> G1 (direct_contraction.data(), direct_contraction.dimension(0), direct_contraction.dimension(1));
    Eigen::Map<Eigen::MatrixXd> G2 (exchange_contraction.data(), exchange_contraction.dimension(0), exchange_contraction.dimension(1));

    // The result is the sum of both contractions, since prefactors were already incorporated
    return G1 + G2;
}


/**
 *  Calculate the RHF electronic energy based on the RHF AO density matrix @param: P, the core Hamiltonian @param: H_core and the Fock matrix @param: F
 */
double RHF::calculateElectronicEnergy(const Eigen::MatrixXd& P, const Eigen::MatrixXd& H_core, const Eigen::MatrixXd& F) const {

    // First, calculate the sum of H_core and F (this saves a contraction)
    Eigen::MatrixXd Z = H_core + F;

    // Convert the matrices Z and P to an Eigen::Tensor<double, 2> P_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> P_tensor (P.data(), P.rows(), P.cols());
    Eigen::TensorMap<Eigen::Tensor<double, 2>> Z_tensor (Z.data(), P.rows(), P.cols());

    // Specify the contraction pair
    // To calculate the electronic energy, we must perform a double contraction
    //      0.5 P(nu mu) Z(mu nu)
    Eigen::array<Eigen::IndexPair<int>, 2> contraction_pair = {Eigen::IndexPair<int>(0, 1), Eigen::IndexPair<int>(1, 0)};

    // Calculate the double contraction (with prefactor 0.5)
    Eigen::Tensor<double, 0> contraction = 0.5 * P_tensor.contract(Z_tensor, contraction_pair);

    // As the double contraction of two matrices is a scalar (a tensor of rank 0), we can access the value as (0).
    return contraction(0);
}



/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given libwint::AOBasis @param: ao_basis, a number of electrons @param: N and an SCF-cycle @param: scf_threshold, MAX_CYCLES
 */
RHF::RHF(const libwint::Molecule& molecule, const libwint::AOBasis& ao_basis, double scf_threshold, size_t MAX_CYCLES) :
    scf_threshold (scf_threshold),
    ao_basis (ao_basis),
    molecule (molecule),
    K (this->ao_basis.calculateNumberOfBasisFunctions()),
    N (this->molecule.get_N()),
    MAX_NUMBER_OF_SCF_CYCLES (MAX_CYCLES)
{

    if ((this->N % 2) != 0) {
        throw std::invalid_argument("The given molecule has an odd number of electrons, which is not compatible with an RFH calculation.");
    }

    if (this->N > 2 * this->K) {
        throw std::invalid_argument("There are too many electrons in the molecule for the given number of spatial orbitals in the AOBasis.");
    }
}



/*
 *  GETTERS
 */

Eigen::VectorXd RHF::get_orbital_energies() const {

    if (!(this->is_converged)) {
        throw std::runtime_error("The RHF procedure isn't converged yet and you are trying to get orbital energies.");
    }

    return this->SCF_solver_ptr->get_orbital_energies();
}

double RHF::get_orbital_energies(size_t index) const {

    if (!(this->is_converged)) {
        throw std::runtime_error("The RHF procedure isn't converged yet and you are trying to get orbital energies.");
    }

    return this->SCF_solver_ptr->get_orbital_energies()(index);
}

Eigen::MatrixXd RHF::get_C_canonical() const {

    if (!(this->is_converged)) {
        throw std::runtime_error("The RHF procedure isn't converged yet and you are trying to get the coefficient matrix.");
    }

    return this->SCF_solver_ptr->get_C_canonical();
}

double RHF::get_electronic_energy() const {

    if (!(this->is_converged)) {
        throw std::runtime_error("The RHF procedure isn't converged yet and you are trying to get the electronic energy.");
    }

    return this->electronic_energy;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the restricted Hartree-Fock equations (i.e. the Roothaan-Hall equations)
 */
void RHF::solve(hf::rhf::solver::SCFSolverType solver_type) {

    // Calculate H_core
    Eigen::MatrixXd H_core = this->ao_basis.get_T() + this->ao_basis.get_V();
    hf::DensityFunction calculateP = [this] (const Eigen::MatrixXd& x) { return this->calculateP(x);};
    hf::TwoElectronMatrixFunction calculateG = [this] (const Eigen::MatrixXd & x, const Eigen::Tensor<double, 4> & y) { return this->calculateG(x,y);};

    switch (solver_type) {

        case hf::rhf::solver::SCFSolverType::PLAIN: {
            auto plain_solver = new hf::rhf::solver::PlainSCFSolver(this->ao_basis.get_S(), H_core, this->ao_basis.get_g(),
                                                                    calculateP, calculateG, this->scf_threshold,
                                                                    this->MAX_NUMBER_OF_SCF_CYCLES);
            plain_solver->solve();
            this->SCF_solver_ptr = plain_solver;  // prevent data from going out of scope
            // we are only assigning this->eigensolver_ptr now, because
            // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }

        case hf::rhf::solver::SCFSolverType::DIIS: {
            auto DIIS_solver = new hf::rhf::solver::PlainSCFSolver(ao_basis.get_S(), H_core, ao_basis.get_g(),
                                                                   calculateP, calculateG, this->scf_threshold,
                                                                   this->MAX_NUMBER_OF_SCF_CYCLES);
            DIIS_solver->solve();
            this->SCF_solver_ptr = DIIS_solver;  // prevent data from going out of scope
            // we are only assigning this->eigensolver_ptr now, because
            // this->solveMatrixEigenvalueProblem only accepts BaseMatrixSolver*
            break;
        }

    }
    this->is_converged = true;
    Eigen::MatrixXd P = calculateP(get_C_canonical());
    this->electronic_energy=calculateElectronicEnergy(P,H_core,H_core+calculateG(P,this->ao_basis.get_g()));





}


/**
 *  Given a number of spatial orbitals @param: K and a number of electrons @param: N, calculated the index of the HOMO in the restricted case
 */
size_t RHF::HOMOIndex() const {

    return this->N / 2 - 1;  // need to subtract 1 because computer indices start at 0
}


/**
 *  Given a number of spatial orbitals @param: K and a number of electrons @param: N, calculated the index of the LUMO in the restricted case
 */
size_t RHF::LUMOIndex() const {

    if (this->N >= 2 * this->K) {
        throw std::invalid_argument("There is no LUMO for the given amount of electrons N and spatial orbitals K");
    }

    return this->HOMOIndex() + 1;
}



}  // namespace rhf
}  // namespace hf
