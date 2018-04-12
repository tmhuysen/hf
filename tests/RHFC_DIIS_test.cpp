#define BOOST_TEST_MODULE "SCFSolver"

#include "RHFC.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain




BOOST_AUTO_TEST_CASE ( constrained_NN_test) {
    // Create a Molecule and an AOBasis
    libwint::Molecule N2 ("../tests/ref_data/NN.xyz");
    libwint::AOBasis ao_basis (N2, "STO-3G");
    ao_basis.calculateIntegrals();

    // Do the SCF cycle
    hf::rhf::RHFC rhfc (N2, ao_basis, 1.0e-12, 10000);
    std::vector<size_t> AO_set = {0,1,2,3,4};
    for(double i = -1;i<1.0;i += 0.1){
        rhfc.solve( hf::rhf::solver::SCFSolverType::PLAIN,AO_set,i);
        double total_energy = rhfc.get_electronic_energy() + N2.calculateInternuclearRepulsionEnergy();
        double pop = rhfc.getPopulation_set();
        std::cout<<std::setprecision(9);
        std::cout <<  i << " & " << -pop + 7 << " & " << total_energy << " \\\\" << std::endl;


    }

}

BOOST_AUTO_TEST_CASE ( constrained_CO_test) {
    // Create a Molecule and an AOBasis
    libwint::Molecule CO ("../tests/ref_data/CO.xyz");
    libwint::AOBasis ao_basis (CO, "STO-3G");
    ao_basis.calculateIntegrals();

    // Do the SCF cycle
    hf::rhf::RHFC rhfc (CO, ao_basis, 1.0e-12, 10000);
    std::vector<size_t> AO_set = {0,1,2,3,4};
    for(double i = -1;i<1.0;i += 0.1){
        rhfc.solve( hf::rhf::solver::SCFSolverType::DIIS,AO_set,i);
        double total_energy = rhfc.get_electronic_energy() + CO.calculateInternuclearRepulsionEnergy();
        double pop = rhfc.getPopulation_set();
        std::cout<<std::setprecision(9);
        std::cout <<  i << " & " << -pop + 6 << " & " << total_energy << " \\\\" << std::endl;


    }

}

BOOST_AUTO_TEST_CASE ( constrained_FF_test) {
    // Create a Molecule and an AOBasis
    libwint::Molecule FF ("../tests/ref_data/FF.xyz");
    libwint::AOBasis ao_basis (FF, "STO-3G");
    ao_basis.calculateIntegrals();

    // Do the SCF cycle
    hf::rhf::RHFC rhfc (FF, ao_basis, 1.0e-12, 100000);
    std::vector<size_t> AO_set = {0,1,2,3,4};
    for(double i = -1;i<1.0;i += 0.1){
        rhfc.solve( hf::rhf::solver::SCFSolverType::PLAIN,AO_set,i);
        double total_energy = rhfc.get_electronic_energy() + FF.calculateInternuclearRepulsionEnergy();
        double pop = rhfc.getPopulation_set();
        std::cout<<std::setprecision(9);
        std::cout <<  i << " & " << -pop + 9 << " & " << total_energy << " \\\\" << std::endl;


    }

}