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
 *  Constructor based on a @param molecule, @param ao_basis, @param scf_threshold and a @param maximum_number_of_iterations.
 */
RHF::RHF(const libwint::Molecule& molecule, const libwint::AOBasis& ao_basis, double scf_threshold, size_t maximum_number_of_iterations) :
    scf_threshold (scf_threshold),
    ao_basis (ao_basis),
    molecule (molecule),
    K (this->ao_basis.calculateNumberOfBasisFunctions()),
    N (this->molecule.get_N()),
    maximum_number_of_iterations (maximum_number_of_iterations)
{

    if ((this->N % 2) != 0) {
        throw std::invalid_argument("The given molecule has an odd number of electrons, which is not compatible with an RFH calculation.");
    }

    if (this->N > 2 * this->K) {
        throw std::invalid_argument("There are too many electrons in the molecule for the given number of spatial orbitals in the AOBasis.");
    }
}



/*
 *  DESTRUCTORS
 */
RHF::~RHF() {
    delete this->SCF_solver_ptr;
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
 *  Solve the restricted Hartree-Fock equations (i.e. the Roothaan-Hall equations). On default, a plain SCF solver is used.
 */
void RHF::solve(hf::rhf::solver::SCFSolverType solver_type) {

    // Before everything, we need to calculate H_core
    Eigen::MatrixXd H_core = this->ao_basis.get_T() + this->ao_basis.get_V();

    hf::DensityFunction calculateP = [this] (const Eigen::MatrixXd& x) { return this->calculateP(x); };
    hf::TwoElectronMatrixFunction calculateG = [this] (const Eigen::MatrixXd& x, const Eigen::Tensor<double, 4>& y) { return this->calculateG(x,y); };

    switch (solver_type) {

        case hf::rhf::solver::SCFSolverType::PLAIN: {
            this->SCF_solver_ptr = new hf::rhf::solver::PlainSCFSolver(this->ao_basis.get_S(), H_core, this->ao_basis.get_g(),
                                                                       calculateP, calculateG, this->scf_threshold,
                                                                       this->maximum_number_of_iterations);
            this->SCF_solver_ptr->solve();
            break;
        }

        case hf::rhf::solver::SCFSolverType::DIIS: {
            this->SCF_solver_ptr = new hf::rhf::solver::DIISSCFSolver(this->ao_basis.get_S(), H_core, this->ao_basis.get_g(),
                                                                      calculateP, calculateG, this->scf_threshold,
                                                                      this->maximum_number_of_iterations);
            this->SCF_solver_ptr->solve();
            break;
        }

    }

    this->is_converged = true;

    // Calculate the converged electronic energy
    Eigen::MatrixXd P = this->calculateP(this->SCF_solver_ptr->get_C_canonical());  // AO density matrix
    Eigen::MatrixXd F = H_core + this->calculateG(P, this->ao_basis.get_g());  // AO Fock matrix
    this->electronic_energy = calculateElectronicEnergy(P, H_core, F);
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
