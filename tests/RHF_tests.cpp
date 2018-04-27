#define BOOST_TEST_MODULE "RHFbase"

#include "hf.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


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