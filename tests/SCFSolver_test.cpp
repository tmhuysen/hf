#define BOOST_TEST_MODULE "SCFSolver"

#include "SCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( h2o_sto3g ) {

    // Specify some data
    constexpr auto xyzfilename = "../../docs/h2o.xyz"; // Anticipate an out-of source build, so we need one level higher in directories
    Wrapper::Molecule water (xyzfilename);
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";

    // Do the SCF cycle
    HF::SCFSolver scf (water, threshold, basis_name);

    // Check the energy
    std::cout << scf.energy << std::endl;
}