#ifndef HF_RHFC_HPP
#define HF_RHFC_HPP

#include "RHF.hpp"


namespace hf {
namespace rhf {


class RHFC : public RHF {
private:
    double population_set;

    /**
     *  Calculates GA for AO_set
     */
    Eigen::MatrixXd calculateGA(std::vector<size_t> AO_set) const;

    /**
     *  Given the RHF AO density matrix @param: P calculate the mulliken population
     */
    double calculatePopulation(const Eigen::MatrixXd& P,std::vector<size_t> AO_set) const;




public:
    // Constructors
    /**
     *  Constructor based on a given libwint::AOBasis @param: ao_basis, a molecule @param: molecule and an SCF-cycle @param: scf_threshold, MAX_CYCLES
     */
    RHFC(const libwint::Molecule& molecule, const libwint::AOBasis& ao_basis, double scf_threshold = 1e-6, size_t MAX_CYCLES = 128);


    // DESTRUCTOR
    ~RHFC(){};


    // Methods
    /**
     *  Solve the restricted Hartree-Fock equations (i.e. the Roothaan-Hall equations)
     */
    void solve(hf::rhf::solver::SCFSolverType solver_type, std::vector<size_t> AO_set, double constraint );
    
    
    // GETTERS
    double getPopulation_set() const { return population_set; }


};



} // namespace hf
} // namespace rhf


#endif //HF_RHFC_HPP
