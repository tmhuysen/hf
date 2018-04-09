#ifndef HF_PLAINSCFSOLVER_HPP
#define HF_PLAINSCFSOLVER_HPP

#include "BaseSCFSolver.hpp"

namespace hf {
namespace rhf {
namespace solver {

class PlainSCFSolver: public BaseSCFSolver {
public:
    // CONSTRUCTOR
    /**
     *   Constructor based on the dimension @param dim of the eigenvalue problem.
     */
    explicit PlainSCFSolver(Eigen::MatrixXd S, double threshold = 1e-6, size_t maximum_number_of_iterations = 128);


    // DESTRUCTOR
    ~PlainSCFSolver() override = default;


    // PUBLIC PURE VIRTUAL METHODS
    /**
     *  Execute the SCF procedure.
     *
     *  If successful, it sets
     *      - @member is_converged to true
     *      - @member C_canonical
     *      - @member orbital_energies
     */
     void solve() override;



};

} // solver
} // rhf
} // hf

#endif //HF_PLAINSCFSOLVER_HPP
