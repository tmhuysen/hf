#ifndef HARTREE_FOCK_SCFSOLVER_HPP
#define HARTREE_FOCK_SCFSOLVER_HPP

#include <string>

#include <libint-wrapper.hpp>

namespace HF {

class SCFSolver {
public:
    Wrapper::Molecule molecule;                 // A molecule instance
                                                //      This has properties atoms, and nelec
    double threshold;                           // Convergence threshold for the SCF procedure
    constexpr unsigned MAX_NO_ITERATIONS = 128;
    std::string basis_name;                     // e.g. "STO-3G" or "6-31G"

    double E;                                   // The converged energy


    // Constructors
    /** Constructor based on a given Wrapper::Molecule molecule and a threshold
     *
     * This automatically starts the SCF procedure
     */
    SCFSolver(Wrapper::Molecule& molecule, double threshold, std::string& basis_name);
};

} // namespace HF


#endif //HARTREE_FOCK_SCFSOLVER_HPP
