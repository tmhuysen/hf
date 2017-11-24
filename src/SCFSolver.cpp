#include "SCFSolver.hpp"
#include "SCF_functions.hpp"


namespace HF {

    /** Constructor based on a given libwrp::Molecule molecule and a threshold
     *
     * This automatically starts the SCF procedure
     */
    SCFSolver::SCFSolver(libwrp::Molecule& molecule, double threshold, std::string& basis_name) :
        molecule(molecule), threshold(threshold), basis_name(basis_name)
    {
        assert(this->molecule.nelec % 2 == 0);  // We have only implemented a restricted Hartree-Fock algorithm

        // Automatically start the SCF procedure
        std::cout << "\nStarting a Hartree-Fock calculation." << std::endl;
        libint2::initialize();

        // Calculate all one- and two-electron integrals, based on my libint wrapper
        libwrp::Basis basis (this->molecule, this->basis_name);

        basis.compute_overlap_integrals();
        basis.compute_kinetic_integrals();
        basis.compute_nuclear_integrals();
        basis.compute_two_electron_integrals();

        // Calculate H_core
        Eigen::MatrixXd H_core = basis.T + basis.V;

        // Solve the generalized eigenvalue problem for H_core to obtain a guess for the density matrix P
        //  H_core should be self-adjoint
        //  S should be positive definite
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes0 (H_core, basis.S);
        Eigen::MatrixXd C = gsaes0.eigenvectors();
        Eigen::MatrixXd P = HF::calculate_P(C, this->molecule.nelec);

        // Initialize the loop parameters
        bool converged = false;
        unsigned iteration_counter = 1;

        while ((! converged) || iteration_counter > this->MAX_NO_ITERATIONS) {
            // Calculate the G-matrix
            Eigen::MatrixXd G = HF::calculate_G(P, basis.tei);

            // Calculate the Fock matrix
            Eigen::MatrixXd F = H_core + G;

            // Solve the Roothaan equation (generalized eigenvalue problem)
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes (F, basis.S);
            C = gsaes.eigenvectors();

            // Calculate an improved density matrix P from the improved coefficient matrix C
            Eigen::MatrixXd P_previous = P; // We will store the previous density matrix
            P = HF::calculate_P(C, this->molecule.nelec);

            // Check for convergence on the density matrix P
            if ((P - P_previous).norm() <= this->threshold) {

                // After the SCF procedure, we end up with canonical spatial orbitals, i.e. the Fock matrix should be diagonal in this basis
                // Let's check if this is the case, within double float precision
                Eigen::MatrixXd f_MO = C.adjoint() * F * C;  // FIXME: we can use the libwrp function for this
                assert(f_MO.isDiagonal());

                converged = true;
                std::cout << "The SCF algorithm has converged after " << iteration_counter << " iterations.\n" << std::endl;

                // After the calculation has converged, calculate the energy as the sum of the electronic energy and the internuclear repulsion energy
                this->energy = HF::calculate_electronic_energy(P, H_core, F) + this->molecule.internuclear_repulsion();

                // Furthermore, add the orbital energies and the coefficient matrix to (this)
                this->orbital_energies = gsaes.eigenvalues();
                this->C_canonical = C;
            }

            // Update the iteration number
            iteration_counter ++;

        } // SCF cycle loop

        libint2::finalize();
    } // constructor

} // namespace HF
