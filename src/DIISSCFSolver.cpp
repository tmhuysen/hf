#include "DIISSCFSolver.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

namespace hf {
namespace rhf {
namespace solver {


/*
 * CONSTRUCTORS
 */

/**
 *  Constructor to initialize the const S, H_core, g, calculateP, calculateG, threshold and maximum_number_of_iterations.
 */
DIISSCFSolver::DIISSCFSolver(const Eigen::MatrixXd S, const Eigen::MatrixXd H_core, const Eigen::Tensor<double, 4> g,
                             const hf::DensityFunction calculateP, const hf::TwoElectronMatrixFunction calculateG,
                             double threshold, size_t maximum_number_of_iterations) :
        BaseSCFSolver(S, H_core, g, calculateP, calculateG, threshold, maximum_number_of_iterations) {
}



/*
 * PUBLIC METHODS
 */

/**
 *  Execute the SCF procedure. direct inversion of the iterative subspace.
 *
 *  If successful, it sets
 *      - @member is_converged to true
 *      - @member C_canonical
 *      - @member orbital_energies
 */
void DIISSCFSolver::solve() {
    // Initialize error and fock vector for the DIIS procedure
    this->fock_vector = {};
    this->error_vector = {};
    // Solve the generalized eigenvalue problem for H_core to obtain a guess for the density matrix P
    //  H_core should be self-adjoint
    //  S should be positive definite
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes0 (this->H_core,this->S);
    Eigen::MatrixXd C = gsaes0.eigenvectors();
    Eigen::MatrixXd P = this->calculateP(C);
    //this->fock_vector.emplace_back(H_core);  // add fock matrix
    //this->error_vector.emplace_back((H_core*P*this->S - this->S*P*H_core));
    double E = this->calculateElectronicEnergy(P,H_core,H_core);
    size_t iteration_counter = 1;
    while (!this->is_converged) {
        std::cout<<std::setprecision(16);
        std::cout<<iteration_counter<<" E :"<<E<<std::endl;
        // Calculate the G-matrix
        Eigen::MatrixXd G = this->calculateG(P, this->g);
        // Calculate the Fock matrix
        Eigen::MatrixXd f_AO = H_core + G;

        double E_prev = E;

        // Fill deques for DIIS procedure
        this->fock_vector.emplace_back(P);  // add fock matrix
        //this->error_vector.emplace_back((f_AO*P*this->S - this->S*P*f_AO));  // add error matrix
        this->error_vector.emplace_back((f_AO*P*this->S - this->S*P*f_AO));  // add error matrix
        if(true){  // Collapse subspace
            //  Initialize B matrix, representation off all errors
            Eigen::MatrixXd B = -1*Eigen::MatrixXd::Ones(error_vector.size()+1,error_vector.size()+1);  // +1 for the multiplier
            B(error_vector.size(),error_vector.size()) = 0;  // last address of the matrix is 0

            for(size_t i = 0; i<error_vector.size();i++){
                for(size_t j = 0; j < error_vector.size();j++){
                    B(i,j) = (this->error_vector[i]*this->error_vector[j]).trace();
                }
            }
            Eigen::VectorXd b = Eigen::VectorXd::Zero(error_vector.size()+1);  // +1 for the multiplier
            b(error_vector.size()) = -1;  // last address is -1
            Eigen::VectorXd coefficients = B.inverse()*b; // calculate the coefficients
            // Recombine previous fock matrix into improved fock matrix
            P = Eigen::MatrixXd::Zero(this->S.cols(),this->S.cols());
            for(size_t i = 0; i<error_vector.size();i++){
                P += coefficients[i]*fock_vector[i];
            }
            if(error_vector.size()==this->max_error_size){
                this->fock_vector.pop_front();
                this->error_vector.pop_front();
            }
            Eigen::MatrixXd G = this->calculateG(P, this->g);
            // Calculate the Fock matrix
            f_AO = H_core + G;

        }
        // Solve the Roothaan equation (generalized eigenvalue problem)
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes (f_AO, this->S);
        C = gsaes.eigenvectors();

        // Calculate an improved density matrix P from the improved coefficient matrix C
        Eigen::MatrixXd P_previous = P; // We will store the previous density matrix
        P = this->calculateP(C);
        E = this->calculateElectronicEnergy(P,H_core,f_AO);

        // Check for convergence on the density matrix P
        if (std::abs(E_prev- E) <= this->threshold) {
            this->is_converged = true;

            // After the SCF procedure, we end up with canonical spatial orbitals, i.e. the Fock matrix should be diagonal in this basis
            // Let's check if this is the case, within double float precision
            Eigen::MatrixXd f_SO = libwint::transformations::transform_AO_to_SO(f_AO, C);
            assert(f_SO.isDiagonal());

            // Furthermore, add the orbital energies and the coefficient matrix to (this)
            this->orbital_energies = gsaes.eigenvalues();
            this->C_canonical = C;

            std::cout<<std::endl<<"SCF ITERATIONS : "<<iteration_counter<<std::endl;
        }
        else {
            iteration_counter ++;
            // If we reach more than this->maximum_number_of_iterations, the system is considered not to be converging
            if (iteration_counter >= this->maximum_number_of_iterations) {
                throw std::runtime_error("The SCF procedure did not converge.");
            }
        }
    }  // SCF cycle loop
}

double DIISSCFSolver::calculateElectronicEnergy(const Eigen::MatrixXd& P, const Eigen::MatrixXd& H_core, const Eigen::MatrixXd& F) const{

    // First, calculate the sum of H_core and F (this saves a contraction)
    Eigen::MatrixXd Z = H_core + F;

    // Convert the matrices Z and P to an Eigen::Tensor<double, 2> P_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> P_tensor (P.data(), P.rows(), P.cols());
    Eigen::TensorMap<Eigen::Tensor<double, 2>> Z_tensor (Z.data(), P.rows(), P.cols());

    // Specify the contraction pair
    // To calculate the electronic energy, we must perform a double contraction
    //      0.5 P(nu mu) Z(mu nu)
    Eigen::array<Eigen::IndexPair<int>, 2> contraction_pair = {Eigen::IndexPair<int>(0, 1), Eigen::IndexPair<int>(1, 0)};

    // Calculate the double contraction (with prefactor 0.5)
    Eigen::Tensor<double, 0> contraction = 0.5 * P_tensor.contract(Z_tensor, contraction_pair);

    // As the double contraction of two matrices is a scalar (a tensor of rank 0), we can access the value as (0).
    return contraction(0);

}


} // solver
} // rhf
} // hf