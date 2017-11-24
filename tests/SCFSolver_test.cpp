#define BOOST_TEST_MODULE "SCFSolver"

#include "SCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


/** Check if two sets of eigenvalues are equal
 */
bool are_equal_eigenvalues(Eigen::VectorXd evals1, Eigen::VectorXd evals2, double tol) {
    return evals1.isApprox(evals2, tol);
}

/** Check if two eigenvectors are equal. This is the case if they are equal up to their sign.
 */
bool are_equal_eigenvectors(Eigen::VectorXd evec1, Eigen::VectorXd evec2, double tol) {
    return (evec1.isApprox(evec2, tol) || evec1.isApprox(-evec2, tol));
}

/** Check if two sets of eigenvectors are equal.
 */
bool are_equal_sets_eigenvectors(Eigen::MatrixXd evecs1, Eigen::MatrixXd evecs2, double tol) {
    auto dim = evecs1.cols();
    for (unsigned i = 0; i < dim; i++) {
        if (! are_equal_eigenvectors(evecs1.col(i), evecs2.col(i), tol)) {
            return false;
        }
    }
    return true;
}



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


    // Supply the reference data from HORTON
    double ref_energy = -74.942080055631;
    Eigen::VectorXd ref_orbital_energies (7);  // 7 BF in STO-3G for water
    ref_orbital_energies << -20.26289322, -1.20969863, -0.54796582, -0.43652631, -0.38758791, 0.47762043, 0.5881361;

    Eigen::MatrixXd C_ref (7, 7);
    C_ref << -9.94434594e-01, -2.39158997e-01,  3.61117086e-17, -9.36837259e-02,  3.73303682e-31, -1.11639152e-01, -9.04958229e-17,
             -2.40970260e-02,  8.85736467e-01, -1.62817254e-16,  4.79589270e-01, -1.93821120e-30,  6.69575233e-01,  5.16088339e-16,
              1.59542752e-18,  5.29309704e-17, -6.07288675e-01, -1.49717339e-16,  8.94470461e-17, -8.85143477e-16,  9.19231270e-01,
             -3.16155527e-03,  8.58957413e-02,  2.89059171e-16, -7.47426286e-01,  2.81871324e-30,  7.38494291e-01,  6.90314422e-16,
              6.65079968e-35,  1.16150362e-32, -2.22044605e-16, -4.06685146e-30, -1.00000000e+00, -1.78495825e-31,  2.22044605e-16,
              4.59373756e-03,  1.44038811e-01, -4.52995183e-01, -3.29475784e-01,  2.16823939e-16, -7.09847234e-01, -7.32462496e-01,
              4.59373756e-03,  1.44038811e-01,  4.52995183e-01, -3.29475784e-01, -2.16823939e-16, -7.09847234e-01,  7.32462496e-01;


    // Do the SCF cycle
    HF::SCFSolver scf (water, threshold, basis_name);


    // Check the calculated results with the reference
    BOOST_CHECK(std::abs(scf.energy - ref_energy) < 1.0e-06);
    BOOST_CHECK(are_equal_eigenvalues(ref_orbital_energies, scf.orbital_energies, 1.0e-06));
    BOOST_CHECK(are_equal_sets_eigenvectors(C_ref, scf.C_canonical, 1.0e-05));
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
