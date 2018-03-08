#ifndef HF_RHF_HPP
#define HF_RHF_HPP


#include <string>

#include <libwint.hpp>



namespace hf {
namespace rhf {


class RHF {
private:
    static const size_t MAX_NUMBER_OF_SCF_CYCLES = 128;
    const double scf_threshold;  // convergence threshold for the SCF procedure
    const libwint::AOBasis& ao_basis;
    const size_t N;  // number of electrons

    Eigen::VectorXd orbital_energies;  // energies of the spatial orbitals (i.e. eigenvalues of the Fock operator)
    Eigen::MatrixXd C_canonical;  // coefficient matrix linking the spatial orbitals to the underlying basis set
                                  // every column represents a spatial orbital in terms of its AOs
                                  // the coefficient matrix is canonical, which means that the Fock matrix in this basis is diagonal

    double electronic_energy;  // the converged energy


    // Methods
    /**
     *  Given a coefficient matrix @param: C, and the number of electrons N, calculate the RHF AO density matrix P
     */
    Eigen::MatrixXd calculateP(const Eigen::MatrixXd& C) const;

    /**
     *  Given the RHF AO density matrix @param: P, and the two-electron integrals @param g in chemist's notation, calculate the RHF G-matrix
     */
    Eigen::MatrixXd calculateG(const Eigen::MatrixXd& P, const Eigen::Tensor<double, 4>& g) const;

    /**
     *  Calculate the RHF electronic energy based on the RHF AO density matrix @param: P, the core Hamiltonian @param: H_core and the Fock matrix @param: F
     */
    double calculateElectronicEnergy(const Eigen::MatrixXd& P, const Eigen::MatrixXd& H_core, const Eigen::MatrixXd& F) const;



public:
    // Constructors
    /**
     *  Constructor based on a given libwint::AOBasis @param: ao_basis, a number of electrons @param: N and an SCF-cycle @param: scf_threshold
     */
    RHF(const libwint::AOBasis& ao_basis, size_t N, double scf_threshold);


    // Getters
    Eigen::VectorXd get_orbital_energies() const { return this->orbital_energies; }
    Eigen::MatrixXd get_C_canonical() const { return this->C_canonical; }
    double get_electronic_energy() const { return this->electronic_energy; }


    // Methods
    /**
     *  Solve the restricted Hartree-Fock equations (i.e. the Roothaan-Hall equations)
     */
    void solve();

    /**
     *  Given a number of spatial orbitals @param: K and a number of electrons @param: N, calculated the index of the HOMO in the restricted case
     */
    static size_t HOMO_index(size_t K, size_t N);

    /**
     *  Given a number of spatial orbitals @param: K and a number of electrons @param: N, calculated the index of the LUMO in the restricted case
     */
    static size_t LUMO_index(size_t K, size_t N);
};



} // namespace hf
} // namespace rhf



#endif // HF_RHF_HPP
