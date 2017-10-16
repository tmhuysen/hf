#include "SCFSolver.hpp"
#include "SCF_functions.hpp"


namespace HF {

    /** Constructor based on a given Wrapper::Molecule molecule and a threshold
     *
     * This automatically starts the SCF procedure
     */
    SCFSolver::SCFSolver(Wrapper::Molecule& molecule, double threshold, std::string& basis_name) :
        molecule(molecule), threshold(threshold), basis_name(basis_name)
    {
        // Automatically start the SCF procedure
        libint2::initialize();

        // Calculate all one- and two-electron integrals, based on my libint wrapper
        Wrapper::Basis basis (this->molecule, this->basis_name);
        auto nbf = basis.nbf();

        Eigen::MatrixXd S = basis.compute_overlap_integrals();
        Eigen::MatrixXd T = basis.compute_kinetic_integrals();
        Eigen::MatrixXd V = basis.compute_nuclear_integrals();
        Eigen::Tensor<double, 4> tei = basis.compute_two_electron_integrals();

        // Calculate Hcore
        Eigen::MatrixXd H_core = T + V;

        // Diagonalize S to obtain X=S^-1/2
        Eigen::MatrixXd X = HF::calculate_X(S);

        // Diagonalize Hcore to obtain a guess for the density matrix P
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes0 (H_core);
        Eigen::MatrixXd C = saes0.eigenvectors();
        Eigen::MatrixXd P = HF::calculate_P(C, this->molecule.nelec);

        // Initialize the loop parameters
        bool converged = false;
        unsigned iteration_number = 1;

        while ((! converged) || iteration_number > this->MAX_NO_ITERATIONS ) {
            // Calculate the G-matrix
            Eigen::MatrixXd G = HF::calculate_G(P, tei);

            // Calculate the Fock matrix
            Eigen::MatrixXd F = H_core + G;

            // Transform the Fock matrix: F'=X^dagger*F*X
            Eigen::MatrixXd F_ = X.adjoint() * F * X;

            // Diagonalize F' to obtain C' and e
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (F_);
            Eigen::MatrixXd C_ = saes.eigenvectors();

            // Calculate the improved coefficient matrix C
            C = X * C_;

            // Calculate an improved density matrix P from the improved coefficient matrix C
            Eigen::MatrixXd P_previous = P; // We will store the previous density matrix
            P = HF::calculate_P(C, this->molecule.nelec);

            // Check for convergence on the density matrix P
            if ((P - P_previous).norm() <= this->threshold) {
                converged = true;

                // After the calculation has converged, calculate the energy
                this->energy = HF::calculate_energy(P, H_core, F);
            }
        } // SCF cycle loop

        libint2::finalize();
    } // constructor



} // namespace HF

