#ifndef HF_BASESCFSOLVER_HPP
#define HF_BASESCFSOLVER_HPP



#include "common.hpp"
#include <libwint.hpp>



namespace hf {
namespace rhf {
namespace solver {


class BaseSCFSolver {
protected:
    const size_t maximum_number_of_iterations;
    const double threshold;

    const Eigen::MatrixXd S;  // Overlap integral matrix of the basis functions
    const Eigen::MatrixXd H_core;  // One electron integral matrix of the basis functions
    const Eigen::Tensor<double ,4> g;  // Two electron integral matrix of the basis functions
    const hf::DensityFunction calculateP;  // Function that calculates the density
    const hf::TwoElectronMatrixFunction calculateG;  // Function that contracts the Two electron integrals with the density

    bool is_converged = false;
    Eigen::VectorXd orbital_energies;  // energies of the spatial orbitals (i.e. eigenvalues of the Fock operator)
    Eigen::MatrixXd C_canonical;  // coefficient matrix linking the spatial orbitals to the underlying basis set
    // every column represents a spatial orbital in terms of its AOs
    // the coefficient matrix is canonical, which means that the Fock matrix in this basis is diagonal


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor to initialize the const @member dim by @param dim.
     */
    explicit BaseSCFSolver(const Eigen::MatrixXd S, const Eigen::MatrixXd H_core, const Eigen::Tensor<double ,4> g,
                           const hf::DensityFunction calculateP,
                           const hf::TwoElectronMatrixFunction calculateG,
                           double threshold = 1e-6, size_t maximum_number_of_iterations = 128);


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
