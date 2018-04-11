#define BOOST_TEST_MODULE "SCFSolver"

#include "RHFC.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain




BOOST_AUTO_TEST_CASE ( NN_test_2 ) {
    // Create a Molecule and an AOBasis
    libwint::Molecule N2 ("../tests/ref_data/NN_2.xyz");
    libwint::AOBasis ao_basis (N2, "STO-3G");
    ao_basis.calculateIntegrals();

    // Do the SCF cycle
    hf::rhf::RHFC rhfc (N2, ao_basis, 1.0e-15, 1000);
    std::vector<size_t> AO_set = {0,1,2,3,4};
    rhfc.solve( hf::rhf::solver::SCFSolverType::DIIS,AO_set,0);
    double total_energy = rhfc.get_electronic_energy() + N2.calculateInternuclearRepulsionEnergy();
    double pop = rhfc.getPopulation_set();
    std::cout<<std::setprecision(16);
    std::cout << total_energy << std::endl;
    std::cout << pop -7 << std::endl;
}

