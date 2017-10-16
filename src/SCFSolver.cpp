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
        Eigen::MatrixXd C_0 = saes0.eigenvectors();
        Eigen::MatrixXd P_0 = calculate_P(C_0, this->molecule.nelec);
        
        // Initialize the loop parameters
        bool converged = false;
        unsigned iteration_number = 1;

        while ((! converged) || iteration_number > this->MAX_NO_ITERATIONS ) {
            //


            converged = true;

        }


        libint2::finalize();

    }



}

