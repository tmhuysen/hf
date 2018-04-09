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
#ifndef HF_RHF_HPP
#define HF_RHF_HPP


#include <string>
#include <libwint.hpp>
#include "SCFSolverType.hpp"
#include "PlainSCFSolver.hpp"
#include "DIISSCFSolver.hpp"
#include "common.hpp"

namespace hf {
namespace rhf {


class RHF {
private:
    const size_t MAX_NUMBER_OF_SCF_CYCLES;
    bool is_converged = false;
    const double scf_threshold;  // convergence threshold for the SCF procedure

    hf::rhf::solver::BaseSCFSolver* SCF_solver_ptr = nullptr;
    const libwint::AOBasis& ao_basis;
    const libwint::Molecule& molecule;

    const size_t K;  // shortcut to this->ao_basis.calculateNumberOfBasisFunctions()
    const size_t N;  // shortcut to this->molecule.get_N()


    //Eigen::VectorXd orbital_energies;  // energies of the spatial orbitals (i.e. eigenvalues of the Fock operator)
    //Eigen::MatrixXd C_canonical;  // coefficient matrix linking the spatial orbitals to the underlying basis set
                                  // every column represents a spatial orbital in terms of its AOs
                                  // the coefficient matrix is canonical, which means that the Fock matrix in this basis is diagonal

    double electronic_energy;  // the converged energy



    // Methods
    /**
     *  Given a coefficient matrix @param: C, and the number of electrons N, calculate the RHF AO density matrix P
     */
    Eigen::MatrixXd calculateP(const Eigen::MatrixXd& C) const;

    /**
     *  Given the RHF AO density matrix @param: P, and the two-electron integrals @param g in chemist's notation, calculate the RHF G-matrix
     */
    Eigen::MatrixXd calculateG(const Eigen::MatrixXd& P, const Eigen::Tensor<double, 4>& g) const;

    /**
     *  Calculate the RHF electronic energy based on the RHF AO density matrix @param: P, the core Hamiltonian @param: H_core and the Fock matrix @param: F
     */
    double calculateElectronicEnergy(const Eigen::MatrixXd& P, const Eigen::MatrixXd& H_core, const Eigen::MatrixXd& F) const;



public:
    // Constructors
    /**
     *  Constructor based on a given libwint::AOBasis @param: ao_basis, a molecule @param: molecule and an SCF-cycle @param: scf_threshold, MAX_CYCLES
     */
    RHF(const libwint::Molecule& molecule, const libwint::AOBasis& ao_basis, double scf_threshold = 1e-6, size_t MAX_CYCLES = 128);


    // DESTRUCTOR
    ~RHF(){delete(this->SCF_solver_ptr);}


    // Getters
    Eigen::VectorXd get_orbital_energies() const;
    double get_orbital_energies(size_t index) const;
    Eigen::MatrixXd get_C_canonical() const;
    double get_electronic_energy() const;


    // Methods
    /**
     *  Solve the restricted Hartree-Fock equations (i.e. the Roothaan-Hall equations)
     */
    void solve(hf::solver::SCFSolverType solver_type);

    /**
     *  Given a number of spatial orbitals @param: K and a number of electrons @param: N, calculated the index of the HOMO in the restricted case
     */
    size_t HOMOIndex() const;

    /**
     *  Given a number of spatial orbitals @param: K and a number of electrons @param: N, calculated the index of the LUMO in the restricted case
     */
    size_t LUMOIndex() const;
};



} // namespace hf
} // namespace rhf



#endif  // HF_RHF_HPP
