#ifndef HF_SCFSOLVER_HPP
#define HF_SCFSOLVER_HPP

#include <string>

#include <libwrp.hpp>


namespace HF {

class SCFSolver {
public:
    Wrapper::Molecule molecule;                 // A molecule instance
                                                //      This has properties atoms, and nelec
    double threshold;                           // Convergence threshold for the SCF procedure
    const unsigned MAX_NO_ITERATIONS = 128;
    std::string basis_name;                     // e.g. "STO-3G" or "6-31G"

    Eigen::VectorXd orbital_energies;           // Energies of the spatial orbitals (i.e. eigenvalues of the diagonal Fock operator)
    Eigen::MatrixXd C_canonical;                // Coefficient matrix linking the spatial orbitals to the underlying (Gaussian) basis set
                                                //      Every column represents a spatial orbital in terms of its AOs
                                                //      The coefficient matrix is canonical, which means that the Fock matrix in this basis is diagonal

    double energy;                              // The converged energy


    // Constructors
    /** Constructor based on a given Wrapper::Molecule molecule and a threshold
     *
     * This automatically starts the restricted SCF procedure
     */
    SCFSolver(Wrapper::Molecule& molecule, double threshold, std::string& basis_name);
};

} // namespace HF


#endif //HF_SCFSOLVER_HPP
