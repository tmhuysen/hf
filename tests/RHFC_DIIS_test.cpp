#define BOOST_TEST_MODULE "SCFSolver"

#include "RHFC.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

BOOST_AUTO_TEST_CASE ( constrained_CO_test) {
    // Create a Molecule and an AOBasis
    libwint::Molecule CO ("../tests/ref_data/CO.xyz");
    libwint::AOBasis ao_basis (CO, "STO-3G");
    ao_basis.calculateIntegrals();
    // Initialize the ref data form Self-consistent methods constrained to a fixed number of particles in a given fragment and its relation to the electronegativity equalization method
    // Andrés Cedillo • Dimitri Van Neck • Patrick Bultinck
    Eigen::Matrix<double,21,3> CO_data;
    CO_data <<  -1.0 , 1.73  , -110.530475 ,
                -0.9 , 1.62  , -110.634766 ,
                -0.8 , 1.50  , -110.737606 ,
                -0.7 , 1.37  , -110.836596 ,
                -0.6 , 1.23  , -110.929219 ,
                -0.5 , 1.08  , -111.012983 ,
                -0.4 , 0.91  , -111.085538 ,
                -0.3 , 0.75  , -111.144760 ,
                -0.2 , 0.57  , -111.188802 ,
                -0.1 , 0.39  , -111.216114 ,
                0    , 0.20  , -111.225446 ,
                0.1  , 0.01  , -111.215842 ,
                0.2  , -0.19 , -111.186619 ,
                0.3  , -0.38 , -111.137356 ,
                0.4  , -0.58 , -111.067872 ,
                0.5  , -0.78 , -110.978210  ,
                0.6  , -0.98 , -110.868621 ,
                0.7  , -1.18 , -110.739544 ,
                0.8  , -1.38 , -110.591599 ,
                0.9  , -1.57 , -110.425574 ,
                1.0  , -1.77 , -110.242423 ;
    // Initialize the constrained RHF.
    hf::rhf::RHFC rhfc (CO, ao_basis, 1.0e-12, 10000);
    // Pick a set of AO's
    std::vector<size_t> AO_set = {0,1,2,3,4};
    std::cout<<std::setprecision(9);
    // Iterate over multipliers
    size_t iteration = 0;
    for(double i = -1;i<1.0;i += 0.1){
        rhfc.solve( hf::rhf::solver::SCFSolverType::DIIS,AO_set,i);
        double total_energy = rhfc.get_electronic_energy() + CO.calculateInternuclearRepulsionEnergy();
        double population = -rhfc.get_population_set() + 6;
        std::cout<<std::setprecision(9);
        BOOST_CHECK(std::abs(total_energy - CO_data(iteration,2)) < 1.0e-2);
        BOOST_CHECK(std::abs(population - CO_data(iteration,1)) < 1.0e-2);
        iteration++;


    }

}

/* was for creating data sets
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
        double pop = rhfc.get_population_set();
        std::cout<<std::setprecision(9);
        std::cout <<  i << " & " << -pop + 9 << " & " << total_energy << " \\\\" << std::endl;


    }

}


BOOST_AUTO_TEST_CASE ( constrained_NN_test) {
    // Create a Molecule and an AOBasis
    libwint::Molecule N2 ("../tests/ref_data/NN.xyz");
    libwint::AOBasis ao_basis (N2, "STO-3G");
    ao_basis.calculateIntegrals();

    // Do the SCF cycle
    hf::rhf::RHFC rhfc (N2, ao_basis, 1.0e-12, 10000);
    std::vector<size_t> AO_set = {0,1,2,3,4};
    for(double i = -1;i<1.0;i += 0.1){
        rhfc.solve( hf::rhf::solver::SCFSolverType::PLAIN,AO_set,i);  // DIIS converges to incorrect basis.
        double total_energy = rhfc.get_electronic_energy() + N2.calculateInternuclearRepulsionEnergy();
        double pop = rhfc.get_population_set();
        std::cout<<std::setprecision(9);
        std::cout <<  i << " & " << -pop + 7 << " & " << total_energy << " \\\\" << std::endl;


    }

}
*/