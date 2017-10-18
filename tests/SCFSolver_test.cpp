#define BOOST_TEST_MODULE "SCFSolver"

#include "SCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( h2_sto3g_szabo ) {
    // In this test case, we will follow section 3.5.2 in Szabo.

    // Specify the data
    const auto xyzfilename = "../../docs/h2.xyz";
    std::string basis_name = "STO-3G";
    double threshold = 1.0e-06;

    // Create a Molecule and a Basis
    Wrapper::Molecule h2 (xyzfilename);

    // Do the SCF cycle
    HF::SCFSolver scf (h2, threshold, basis_name);


    // The converged coefficient matrix is listed as
    Eigen::MatrixXd P_converged_ref (2, 2);
    double p11 = 1 / (1 + 0.6593);
    P_converged_ref << p11, p11,
                       p11, p11;

    std::cout << "P_converged_ref:" << std::endl << P_converged_ref << std::endl << std::endl;

    // Check the energy
    BOOST_CHECK(std::abs(scf.energy - (-1.1167)) < 1.0e-04); // Reference data from Szabo
}


BOOST_AUTO_TEST_CASE ( h2o_sto3g ) {

    // Specify some data
    constexpr auto xyzfilename = "../../docs/h2o.xyz"; // Anticipate an out-of source build, so we need one level higher in directories
    Wrapper::Molecule water (xyzfilename);
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";

    // Do the SCF cycle
    HF::SCFSolver scf (water, threshold, basis_name);

    // Check the energy
    BOOST_CHECK(std::abs(scf.energy - (-74.942080)) < 1.0e-06); // Reference data from horton
}
