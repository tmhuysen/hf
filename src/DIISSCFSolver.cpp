#include "DIISSCFSolver.hpp"
namespace hf {
namespace rhf {
namespace solver {


DIISSCFSolver::DIISSCFSolver(const Eigen::MatrixXd S, const Eigen::MatrixXd H_core, const Eigen::Tensor<double, 4> g,
                             const hf::DensityFunction calculateP, const hf::TwoElectronMatrixFunction calculateG,
                             double threshold, size_t maximum_number_of_iterations) : BaseSCFSolver(S, H_core, g,
                                                                                                    calculateP,
                                                                                                    calculateG,
                                                                                                    threshold,
                                                                                                    maximum_number_of_iterations) {
    this->fock_vector = {};
    this->error_vector = {};

}

void DIISSCFSolver::solve() {
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes0 (this->H_core,this->S);
    Eigen::MatrixXd C = gsaes0.eigenvectors();
    Eigen::MatrixXd P = this->calculateP(C);


    size_t iteration_counter = 1;
    while (!this->is_converged) {
        // Calculate the G-matrix
        Eigen::MatrixXd G = this->calculateG(P, this->g);
        // Calculate the Fock matrix
        Eigen::MatrixXd f_AO = H_core + G;

        // Fill deques for DIIS procedure
        this->fock_vector.emplace_back(f_AO);  // add fock matrix
        this->error_vector.emplace_back((f_AO*P*this->S - this->S*P*f_AO));  // add error matrix
        if(error_vector.size()==this->max_error_size){  // Collapse subspace
            //  Initialize B matrix, representation off all errors
            Eigen::MatrixXd B = -1*Eigen::MatrixXd::Ones(this->max_error_size+1,this->max_error_size+1);  // +1 for the multiplier
            B(this->max_error_size,this->max_error_size) = 0;  // last address of the matrix is 0

            for(size_t i = 0; i<this->max_error_size;i++){
                for(size_t j = 0; j < this->max_error_size;j++){
                    B(i,j) = (this->error_vector[i]*this->error_vector[j]).trace();
                }
            }
            Eigen::VectorXd b = Eigen::VectorXd::Zero(this->max_error_size+1);  // +1 for the multiplier
            b(this->max_error_size) = -1;  // last address is -1
            Eigen::VectorXd coefficients = B.inverse()*b; // calculate the coefficients
            // Recombine previous fock matrix into improved fock matrix
            f_AO = Eigen::MatrixXd::Zero(this->S.cols(),this->S.cols());
            for(size_t i = 0; i<max_error_size;i++){
                f_AO += coefficients[i]*fock_vector[i];
            }

            // Remove the oldest entries. So we collapse every iteration
            this->fock_vector.pop_front();
            this->error_vector.pop_front();
        }
        // Solve the Roothaan equation (generalized eigenvalue problem)
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gsaes (f_AO, this->S);
        C = gsaes.eigenvectors();

        // Calculate an improved density matrix P from the improved coefficient matrix C
        Eigen::MatrixXd P_previous = P; // We will store the previous density matrix
        P = this->calculateP(C);

        // Check for convergence on the density matrix P
        if ((P - P_previous).norm() <= this->threshold) {
            this->is_converged = true;

            // After the SCF procedure, we end up with canonical spatial orbitals, i.e. the Fock matrix should be diagonal in this basis
            // Let's check if this is the case, within double float precision
            Eigen::MatrixXd f_SO = libwint::transformations::transform_AO_to_SO(f_AO, C);
            assert(f_SO.isDiagonal());

            // Furthermore, add the orbital energies and the coefficient matrix to (this)
            this->orbital_energies = gsaes.eigenvalues();
            this->C_canonical = C;

            std::cout<<std::endl<<"ITERATINOS : "<<iteration_counter<<std::endl;
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
} // solver
} // rhf
} // hf