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
 *  Constructor based on a given libwint::AOBasis @param: ao_basis, a number of electrons @param: N and an SCF-cycle @param: scf_threshold
 */
RHF::RHF(const libwint::Molecule& molecule, const libwint::AOBasis& ao_basis, double scf_threshold) :
    scf_threshold (scf_threshold),
    ao_basis (ao_basis),
    molecule (molecule),
    K (this->ao_basis.calculateNumberOfBasisFunctions()),
    N (this->molecule.get_N())
{

    if ((this->N % 2) != 0) {
        throw std::invalid_argument("The given molecule has an odd number of electrons, which is not compatible with an RFH calculation.");
    }

    if (this->N > 2 * this->K) {
        throw std::invalid_argument("There are too many electrons in the molecule for the given number of spatial orbitals in the AOBasis.");
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the restricted Hartree-Fock equations (i.e. the Roothaan-Hall equations)
 */
void RHF::solve() {

    // Calculate H_core
    Eigen::MatrixXd H_core = this->ao_basis.get_T() + this->ao_basis.get_V();

    // Solve the generalized eigenvalue problem for H_core to obtain a guess for the density matrix P
    //  H_core should be self-adjoint
    //  S should be positive definite
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes0 (H_core, this->ao_basis.get_S());
    Eigen::MatrixXd C = gsaes0.eigenvectors();
    Eigen::MatrixXd P = this->calculateP(C);

    // Initialize the loop parameters
    bool converged = false;
    size_t iteration_counter = 1;

    while (!converged) {
        // Calculate the G-matrix
        Eigen::MatrixXd G = this->calculateG(P, this->ao_basis.get_g());

        // Calculate the Fock matrix
        Eigen::MatrixXd f_AO = H_core + G;

        // Solve the Roothaan equation (generalized eigenvalue problem)
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes (f_AO, this->ao_basis.get_S());
        C = gsaes.eigenvectors();

        // Calculate an improved density matrix P from the improved coefficient matrix C
        Eigen::MatrixXd P_previous = P; // We will store the previous density matrix
        P = this->calculateP(C);

        // Check for convergence on the density matrix P
        if ((P - P_previous).norm() <= this->scf_threshold) {
            converged = true;

            // After the SCF procedure, we end up with canonical spatial orbitals, i.e. the Fock matrix should be diagonal in this basis
            // Let's check if this is the case, within double float precision
            Eigen::MatrixXd f_SO = libwint::transformations::transform_AO_to_SO(f_AO, C);
            assert(f_SO.isDiagonal());

            // After the calculation has converged, calculate the electronic energy
            this->electronic_energy = this->calculateElectronicEnergy(P, H_core, f_AO);

            // Furthermore, add the orbital energies and the coefficient matrix to (this)
            this->orbital_energies = gsaes.eigenvalues();
            this->C_canonical = C;
        }
        else {  // not converged yet
            iteration_counter ++;

            // If we reach more than this->MAX_NUMBER_OF_SCF_CYCLES, the system is considered not to be converging
            if (iteration_counter >= this->MAX_NUMBER_OF_SCF_CYCLES) {
                throw std::runtime_error("The SCF procedure did not converge.");
            }
        }
    }  // SCF cycle loop
}


/**
 *  Given a number of spatial orbitals @param: K and a number of electrons @param: N, calculated the index of the HOMO in the restricted case
 */
size_t RHF::HOMOIndex() {

    return this->N / 2 - 1;  // need to subtract 1 because computer indices start at 0
}


/**
 *  Given a number of spatial orbitals @param: K and a number of electrons @param: N, calculated the index of the LUMO in the restricted case
 */
size_t RHF::LUMOIndex() {

    if (this->N >= 2 * this->K) {
        throw std::invalid_argument("There is no LUMO for the given amount of electrons N and spatial orbitals K");
    }

    return this->HOMOIndex() + 1;
}



}  // namespace rhf
}  // namespace hf
