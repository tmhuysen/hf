#define BOOST_TEST_MODULE "SCFSolver"

#include "SCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( h2_sto3g_szabo ) {
    // In this test case, we will follow section 3.5.2 in Szabo.

    // Specify the data
    const std::string xyzfilename = "../tests/reference/h2.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
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

    // Check the converged density matrix
    // FIXME: Check the converged density matrix

    // Check the energy
    BOOST_CHECK(std::abs(scf.energy - (-1.1167)) < 1.0e-04); // Reference data from Szabo
}


BOOST_AUTO_TEST_CASE ( h2o_sto3g ) {

    // Specify some data
    const std::string xyzfilename = "../tests/reference/h2o.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    Wrapper::Molecule water (xyzfilename);
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";

    // Do the SCF cycle
    HF::SCFSolver scf (water, threshold, basis_name);

    // Check the energy
    BOOST_CHECK(std::abs(scf.energy - (-74.942080)) < 1.0e-06); // Reference data from horton
}


BOOST_AUTO_TEST_CASE ( crawdad_h2o_sto3g ) {

    // This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.

    // Specify the input file, energy threshold and basis set
    const std::string xyzfilename = "../tests/reference/h2o_crawdad.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    Wrapper::Molecule water (xyzfilename);
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";

    // Check if the internuclear distance between O and H is really 1.1 A (= 2.07869 bohr), as specified in the text
    BOOST_CHECK(std::abs(Wrapper::distance(water.atoms[0], water.atoms[1]) - 2.07869) < 1.0e-3);

    // Do the SCF cycle
    HF::SCFSolver scf (water, threshold, basis_name);

    // Check the energy
    BOOST_CHECK(std::abs(scf.energy - (-74.9420799281920)) < threshold);  // Reference data from crawdad
}


BOOST_AUTO_TEST_CASE ( crawdad_ch4_sto3g ) {

    // This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.

    // Specify the input file, energy threshold and basis set
    const std::string xyzfilename = "../tests/reference/ch4_crawdad.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    Wrapper::Molecule methane (xyzfilename);
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";

    // Check if the internuclear distance between C and H is really around 2.05 bohr, which is the bond distance Wikipedia (108.7 pm) specifies
    BOOST_CHECK(std::abs(Wrapper::distance(methane.atoms[0], methane.atoms[1]) - 2.05) < 1.0e-1);

    // Do the SCF cycle
    HF::SCFSolver scf (methane, threshold, basis_name);

    // Check the energy
    BOOST_CHECK(std::abs(scf.energy - (-39.726850324347)) < threshold);
}
