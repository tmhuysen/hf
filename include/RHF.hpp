#ifndef HF_RHF_HPP
#define HF_RHF_HPP

#include <libwrp.hpp>
#include <string>

namespace hf {
namespace rhf {

class RHF {
public:
    libwrp::Basis &basis;                       // A reference to a Basis object. Post-SCF methods transform the integrals from AO basis to SO basis, so it's better not to re-calculate them.
    // Contains: .S, .T, .V, .tei

    double threshold;                           // Convergence threshold for the SCF procedure
    const size_t MAX_NO_ITERATIONS = 128;

    Eigen::VectorXd orbital_energies;           // Energies of the spatial orbitals (i.e. eigenvalues of the diagonal Fock operator)
    Eigen::MatrixXd C_canonical;                // Coefficient matrix linking the spatial orbitals to the underlying (Gaussian) basis set
    //      Every column represents a spatial orbital in terms of its AOs
    //      The coefficient matrix is canonical, which means that the Fock matrix in this basis is diagonal

    double energy;                              // The converged energy


    /** Constructor based on a given libwrp::Basis and an SCF-cycle threshold
     *
     *      This automatically starts the restricted SCF procedure
     */
    RHF(libwrp::Basis &basis, double threshold);
};

} // namespace hf
} // namespace rhf


#endif // HF_RHF_HPP
