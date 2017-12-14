#define BOOST_TEST_MODULE "SCFSolver"

#include "RHF.hpp"
#include "utility.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( reference ) {

    // Specify some data, create a Molecule and a Basis
    const std::string xyzfilename = "../tests/reference/h2_szabo.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    std::string basis_name = "STO-3G";
    double threshold = 1.0e-06;

    libwrp::Molecule h2 (xyzfilename);
    libwrp::Basis basis (h2, basis_name);

    // Create an SCFSolver instance - this automatically performs the SCF cycle (not that this is needed in this test)
    hf::rhf::RHF rhf (basis, threshold);

    // Test if a reference to the Basis object is made
    BOOST_CHECK_EQUAL(&basis, &rhf.basis);
}


BOOST_AUTO_TEST_CASE ( h2_sto3g_szabo ) {
    // In this test case, we will follow section 3.5.2 in Szabo.

    // Specify the data
    const std::string xyzfilename = "../tests/reference/h2_szabo.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    std::string basis_name = "STO-3G";
    double threshold = 1.0e-06;

    // Create a Molecule and a Basis
    libwrp::Molecule h2 (xyzfilename);
    libwrp::Basis basis (h2, basis_name);

    // Do the SCF cycle
    hf::rhf::RHF rhf (basis, threshold);


    // The converged coefficient matrix is listed as
    Eigen::MatrixXd P_converged_ref (2, 2);
    double p11 = 1 / (1 + 0.6593);
    P_converged_ref << p11, p11,
                       p11, p11;

    // Check the energy
    BOOST_CHECK(std::abs(rhf.energy - (-1.1167)) < 1.0e-04); // Reference data from Szabo
}


BOOST_AUTO_TEST_CASE ( h2o_sto3g ) {

    // Specify some data
    const std::string xyzfilename = "../tests/reference/h2o.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";

    libwrp::Molecule water (xyzfilename);
    libwrp::Basis basis (water, basis_name);


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
    hf::rhf::RHF rhf (basis, threshold);


    // Check the calculated results with the reference
    BOOST_CHECK(std::abs(rhf.energy - ref_energy) < 1.0e-06);
    BOOST_CHECK(hf::utility::are_equal_eigenvalues(ref_orbital_energies, rhf.orbital_energies, 1.0e-06));
    BOOST_CHECK(hf::utility::are_equal_sets_eigenvectors(C_ref, rhf.C_canonical, 1.0e-05));
}


BOOST_AUTO_TEST_CASE ( crawdad_h2o_sto3g ) {

    // This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.

    // Specify the input file, energy threshold and basis set
    const std::string xyzfilename = "../tests/reference/h2o_crawdad.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";

    libwrp::Molecule water (xyzfilename);
    libwrp::Basis basis (water, basis_name);


    // Check if the internuclear distance between O and H is really 1.1 A (= 2.07869 bohr), as specified in the text
    BOOST_CHECK(std::abs(libwrp::distance(water.atoms[0], water.atoms[1]) - 2.07869) < 1.0e-3);

    // Do the SCF cycle
    hf::rhf::RHF rhf (basis, threshold);

    // Check the energy
    BOOST_CHECK(std::abs(rhf.energy - (-74.9420799281920)) < threshold);  // Reference data from crawdad
}


BOOST_AUTO_TEST_CASE ( crawdad_ch4_sto3g ) {

    // This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.

    // Specify the input file, energy threshold and basis set
    const std::string xyzfilename = "../tests/reference/ch4_crawdad.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";

    libwrp::Molecule methane (xyzfilename);
    libwrp::Basis basis (methane, basis_name);


    // Check if the internuclear distance between C and H is really around 2.05 bohr, which is the bond distance Wikipedia (108.7 pm) specifies
    BOOST_CHECK(std::abs(libwrp::distance(methane.atoms[0], methane.atoms[1]) - 2.05) < 1.0e-1);

    // Do the SCF cycle
    hf::rhf::RHF rhf (basis, threshold);

    // Check the energy
    BOOST_CHECK(std::abs(rhf.energy - (-39.726850324347)) < threshold);
}


BOOST_AUTO_TEST_CASE ( h2_sto6g ) {

    // We have some reference data from olsens: H2 with HF/STO-6G orbitals
    double E_el_rhf_ref = -1.838434256;

    // Test the H2 results
    const std::string xyzfilename = "../tests/reference/h2_olsens.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    std::string basis_name = "STO-6G";
    double threshold = 1.0e-06;

    libwrp::Molecule h2 (xyzfilename);
    libwrp::Basis basis (h2, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    double E_el_rhf = rhf.energy - rhf.basis.molecule.internuclear_repulsion();  // we only have reference data for the electronic repulsion
    BOOST_CHECK(std::abs(E_el_rhf_ref - E_el_rhf) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( h2_631gdp ) {

    // We have some reference data from olsens: H2 with HF/6-31G** orbitals
    double E_el_rhf_ref = -1.84444667247;

    // Test the H2 results
    const std::string xyzfilename = "../tests/reference/h2_olsens.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    std::string basis_name = "6-31g**";
    double threshold = 1.0e-06;

    libwrp::Molecule h2 (xyzfilename);
    libwrp::Basis basis (h2, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    double E_el_rhf = rhf.energy - rhf.basis.molecule.internuclear_repulsion();  // we only have reference data for the electronic repulsion
    BOOST_CHECK(std::abs(E_el_rhf_ref - E_el_rhf) < 1.0e-06);
    std::cout << E_el_rhf << std::endl;
}


BOOST_AUTO_TEST_CASE ( lih_sto6g ) {

    // We have some reference data from olsens: LiH with HF/STO-6G orbitals
    double E_el_rhf_ref = -8.9472891719;

    // Test the LiH results
    const std::string xyzfilename = "../tests/reference/lih_olsens.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    std::string basis_name = "sto-6g";
    double threshold = 1.0e-06;

    libwrp::Molecule lih (xyzfilename);
    libwrp::Basis basis (lih, basis_name);
    hf::rhf::RHF rhf (basis, threshold);

    double E_el_rhf = rhf.energy - rhf.basis.molecule.internuclear_repulsion();  // we only have reference data for the electronic repulsion
    BOOST_CHECK(std::abs(E_el_rhf_ref - E_el_rhf) < 1.0e-06);
}