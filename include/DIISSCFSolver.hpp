#ifndef HF_DIISSCFSOLVER_HPP
#define HF_DIISSCFSOLVER_HPP


#include <deque>
#include "BaseSCFSolver.hpp"


namespace hf {
namespace rhf {
namespace solver {

class DIISSCFSolver: public BaseSCFSolver {
private:
    std::deque<Eigen::MatrixXd> fock_vector;
    std::deque<Eigen::MatrixXd> error_vector;

    const size_t max_error_size = 5;
public:
    // CONSTRUCTOR
    /**
     *   Constructor based on the dimension @param dim of the eigenvalue problem.
     */
    explicit DIISSCFSolver(const Eigen::MatrixXd S, const Eigen::MatrixXd H_core, const Eigen::Tensor<double ,4> g,
                            const hf::DensityFunction calculateP,
                            const hf::TwoElectronMatrixFunction calculateG,
                            double threshold = 1e-6, size_t maximum_number_of_iterations = 128);

    // DESTRUCTOR
    ~DIISSCFSolver() override = default;


    // PUBLIC METHODS
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


#endif //HF_DIISSCFSOLVER_HPP
