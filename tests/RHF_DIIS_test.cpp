#define BOOST_TEST_MODULE "SCFSolver"

#include "hf.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( constructor ) {

    // Check if we only accept even numbers of electrons
    libwint::Molecule h2_cation ("../tests/ref_data/h2_szabo.xyz", +1);  // H2+
    libwint::AOBasis ao_basis1 (h2_cation, "STO-3G");
    ao_basis1.calculateOverlapIntegrals();  // need to calculate one of the integrals to be able to access the number of basis functions

    BOOST_CHECK_THROW(hf::rhf::RHF (h2_cation, ao_basis1, 1.0e-06), std::invalid_argument);


    // Check if we don't accept molecules with too many electrons
    libwint::Molecule h2_many_electrons ("../tests/ref_data/h2_szabo.xyz", -10);  // H2 10-
    libwint::AOBasis ao_basis2 (h2_many_electrons, "STO-3G");
    ao_basis2.calculateOverlapIntegrals();  // need to calculate one of the integrals to be able to access the number of basis functions

    BOOST_CHECK_THROW(hf::rhf::RHF (h2_many_electrons, ao_basis2, 1.0e-06), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( h2_sto3g_szabo ) {

    // In this test case, we will follow section 3.5.2 in Szabo.
    double ref_total_energy = -1.1167;


    // Create a Molecule and an AOBasis
    libwint::Molecule h2 ("../tests/ref_data/h2_szabo.xyz");
    libwint::AOBasis ao_basis (h2, "STO-3G");
    ao_basis.calculateIntegrals();

    // Do the SCF cycle
    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve( hf::rhf::solver::SCFSolverType::DIIS);
    double total_energy = rhf.get_electronic_energy() + h2.calculateInternuclearRepulsionEnergy();

    std::cout << total_energy << std::endl;
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-04);
}


BOOST_AUTO_TEST_CASE ( h2o_sto3g_horton ) {

    // We have some reference data from horton
    double ref_total_energy = -74.942080055631;

    Eigen::VectorXd ref_orbital_energies (7);  // the STO-3G basisset has 7 basis functions for water
    ref_orbital_energies << -20.26289322, -1.20969863, -0.54796582, -0.43652631, -0.38758791, 0.47762043, 0.5881361;

    Eigen::MatrixXd ref_C (7, 7);
    ref_C << -9.94434594e-01, -2.39158997e-01,  3.61117086e-17, -9.36837259e-02,  3.73303682e-31, -1.11639152e-01, -9.04958229e-17,
             -2.40970260e-02,  8.85736467e-01, -1.62817254e-16,  4.79589270e-01, -1.93821120e-30,  6.69575233e-01,  5.16088339e-16,
              1.59542752e-18,  5.29309704e-17, -6.07288675e-01, -1.49717339e-16,  8.94470461e-17, -8.85143477e-16,  9.19231270e-01,
             -3.16155527e-03,  8.58957413e-02,  2.89059171e-16, -7.47426286e-01,  2.81871324e-30,  7.38494291e-01,  6.90314422e-16,
              6.65079968e-35,  1.16150362e-32, -2.22044605e-16, -4.06685146e-30, -1.00000000e+00, -1.78495825e-31,  2.22044605e-16,
              4.59373756e-03,  1.44038811e-01, -4.52995183e-01, -3.29475784e-01,  2.16823939e-16, -7.09847234e-01, -7.32462496e-01,
              4.59373756e-03,  1.44038811e-01,  4.52995183e-01, -3.29475784e-01, -2.16823939e-16, -7.09847234e-01,  7.32462496e-01;


    // Do our own RHF calculation
    libwint::Molecule water ("../tests/ref_data/h2o.xyz");
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (water, ao_basis, 1.0e-06);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);

    double total_energy = rhf.get_electronic_energy() + water.calculateInternuclearRepulsionEnergy();


    // Check the calculated results with the reference
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
    BOOST_CHECK(cpputil::linalg::areEqualEigenvalues(ref_orbital_energies, rhf.get_orbital_energies(), 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqualSetsOfEigenvectors(ref_C, rhf.get_C_canonical(), 1.0e-05));
}


BOOST_AUTO_TEST_CASE ( crawdad_h2o_sto3g ) {

    // This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.
    double ref_total_energy = -74.9420799281920;


    // Do our own RHF calculation
    libwint::Molecule water ("../tests/ref_data/h2o_crawdad.xyz");
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateIntegrals();


    // Check if the internuclear distance between O and H is really 1.1 A (= 2.07869 bohr), as specified in the text
    BOOST_REQUIRE(std::abs(water.calculateInternuclearDistance(0, 1) - 2.07869) < 1.0e-4);


    // Do the SCF cycle
    hf::rhf::RHF rhf (water, ao_basis, 1.0e-06);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    double total_energy = rhf.get_electronic_energy() + water.calculateInternuclearRepulsionEnergy();


    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( crawdad_ch4_sto3g ) {

    // This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.
    double ref_total_energy = -39.726850324347;

    // Do our own RHF calculation
    libwint::Molecule methane ("../tests/ref_data/ch4_crawdad.xyz");
    libwint::AOBasis ao_basis (methane, "STO-3G");
    ao_basis.calculateIntegrals();

    // Check if the internuclear distance between C and H is really around 2.05 bohr, which is the bond distance Wikipedia (108.7 pm) specifies
    BOOST_CHECK(std::abs(methane.calculateInternuclearDistance(0, 1) - 2.05) < 1.0e-1);


    // Do the SCF cycle
    hf::rhf::RHF rhf (methane, ao_basis, 1.0e-06);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);
    double total_energy = rhf.get_electronic_energy() + methane.calculateInternuclearRepulsionEnergy();


    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( h2_sto6g ) {

    // We have some reference data from olsens: H2@RHF//STO-6G orbitals
    double ref_electronic_energy = -1.838434256;


    // Do our own RHF calculation
    libwint::Molecule h2 ("../tests/ref_data/h2_olsens.xyz");
    libwint::AOBasis ao_basis (h2, "STO-6G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);


    BOOST_CHECK(std::abs(rhf.get_electronic_energy() - ref_electronic_energy) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( h2_631gdp ) {

    // We have some reference data from olsens: H2@RHF//6-31G** orbitals
    double ref_electronic_energy = -1.84444667247;


    // Do our own RHF calculation
    libwint::Molecule h2 ("../tests/ref_data/h2_olsens.xyz");
    libwint::AOBasis ao_basis (h2, "6-31g**");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (h2, ao_basis, 1.0e-06);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);


    BOOST_CHECK(std::abs(rhf.get_electronic_energy() - ref_electronic_energy) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( lih_sto6g ) {

    // We have some reference data from olsens: LiH@RHF//STO-6G orbitals
    double ref_electronic_energy = -8.9472891719;


    // Do our own RHF calculation
    libwint::Molecule lih ("../tests/ref_data/lih_olsens.xyz");
    libwint::AOBasis ao_basis (lih, "STO-6G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (lih, ao_basis, 1.0e-06);
    rhf.solve(hf::rhf::solver::SCFSolverType::DIIS);


    BOOST_CHECK(std::abs(rhf.get_electronic_energy() - ref_electronic_energy) < 1.0e-06);
}



BOOST_AUTO_TEST_CASE ( homo ) {

    // Create an RHF object to test the HOMOIndex function on
    libwint::Molecule water ("../tests/ref_data/h2o_crawdad.xyz");
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateOverlapIntegrals();  // need to calculate one of the integrals to be able to access the number of basis functions

    hf::rhf::RHF rhf (water, ao_basis, 1.0e-06);


    // In this case, K=7 and N=10, so the index of the HOMO should be 4
    BOOST_CHECK_EQUAL(rhf.HOMOIndex(), 4);
}


BOOST_AUTO_TEST_CASE ( lumo ) {

    // Create an RHF object to test the LUMOIndex function on
    libwint::Molecule water ("../tests/ref_data/h2o_crawdad.xyz");
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateOverlapIntegrals();  // need to calculate one of the integrals to be able to access the number of basis functions

    hf::rhf::RHF rhf (water, ao_basis, 1.0e-06);


    // In this case, K=7 and N=10, so the index of the LUMO should be 5
    BOOST_CHECK_EQUAL(rhf.LUMOIndex(), 5);
}

BOOST_AUTO_TEST_CASE ( covergence_test ) {

    // Test to see far apart NO converges

    // Do our own RHF calculation
    libwint::Molecule lih ("../tests/ref_data/NO.xyz",1);
    libwint::AOBasis ao_basis (lih, "STO-6G");
    ao_basis.calculateIntegrals();

    hf::rhf::RHF rhf (lih, ao_basis, 1.0e-06);
    // DIIS should converge
    BOOST_CHECK_NO_THROW(rhf.solve(hf::rhf::solver::SCFSolverType::DIIS));

}