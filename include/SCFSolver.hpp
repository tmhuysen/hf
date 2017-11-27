#ifndef HF_SCFSOLVER_HPP
#define HF_SCFSOLVER_HPP

#include <libwrp.hpp>
#include <string>


namespace hf {

class SCFSolver {
public:
    libwrp::Basis& basis;                       // A reference to a Basis object. Post_SCF methods transform the integrals from AO basis to SO basis, so it's better not to re-calculate them.
                                                // Contains: .S, .T, .V, .tei

    double threshold;                           // Convergence threshold for the SCF procedure
    const unsigned MAX_NO_ITERATIONS = 128;

    Eigen::VectorXd orbital_energies;           // Energies of the spatial orbitals (i.e. eigenvalues of the diagonal Fock operator)
    Eigen::MatrixXd C_canonical;                // Coefficient matrix linking the spatial orbitals to the underlying (Gaussian) basis set
                                                //      Every column represents a spatial orbital in terms of its AOs
                                                //      The coefficient matrix is canonical, which means that the Fock matrix in this basis is diagonal

    double energy;                              // The converged energy


    /** Constructor based on a given libwrp::Basis and an SCF-cycle threshold
     *
     *      This automatically starts the restricted SCF procedure
     */
    SCFSolver(libwrp::Basis& basis, double threshold);
};

} // namespace hf


#endif //HF_SCFSOLVER_HPP
