#ifndef HF_BASESCFSOLVER_HPP
#define HF_BASESCFSOLVER_HPP



#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>




namespace hf {
namespace rhf {
namespace solver {


class BaseSCFSolver {
protected:
    const size_t maximum_number_of_iterations;
    const double threshold;

    const Eigen::MatrixXd S;
    std::function<Eigen::MatrixXd(const Eigen::VectorXd &)> calculateP;
    std::function<Eigen::MatrixXd(const Eigen::VectorXd &, const Eigen::Tensor<double, 4> &)> calculateG;

    bool is_converged = false;
    Eigen::VectorXd orbital_energies;
    Eigen::MatrixXd C_canonical;


public:
    // DESTRUCTOR
    virtual ~BaseSCFSolver() = default;


    // PUBLIC PURE VIRTUAL METHODS
    /**
     *  Execute the SCF procedure.
     *
     *  If successful, it sets
     *      - @member is_converged to true
     *      - @member C_canonical
     *      - @member orbital_energies
     */
    virtual void solve() = 0;


    // GETTERS
    Eigen::VectorXd get_orbital_energies() const;
    Eigen::MatrixXd get_C_canonical() const;
};


} // solver
} // rhf
} // hf
#endif //HF_BASESCFSOLVER_HPP
