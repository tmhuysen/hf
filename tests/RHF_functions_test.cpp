#define BOOST_TEST_MODULE "SCF_functions"

#include "RHF_functions.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( homo ) {

        size_t K1 = 4;
        size_t K2 = 3;

        size_t N1 = 2;
        size_t N2 = 4;
        size_t N3 = 6;
        size_t N4 = 8;
        size_t N5 = 10;
        size_t N_odd = 3;

        BOOST_CHECK_EQUAL(hf::rhf::HOMO_index(K1, N1), 0);
        BOOST_CHECK_EQUAL(hf::rhf::HOMO_index(K1, N2), 1);
        BOOST_CHECK_EQUAL(hf::rhf::HOMO_index(K1, N3), 2);
        BOOST_CHECK_EQUAL(hf::rhf::HOMO_index(K1, N4), 3);
        BOOST_REQUIRE_THROW(hf::rhf::HOMO_index(K1, N5), std::invalid_argument);  // Cannot place more than 8 electrons in 4 orbitals
        BOOST_REQUIRE_THROW(hf::rhf::HOMO_index(K1, N_odd), std::invalid_argument);  // The unrestricted case is not supported

        BOOST_CHECK_EQUAL(hf::rhf::HOMO_index(K2, N1), 0);
        BOOST_CHECK_EQUAL(hf::rhf::HOMO_index(K2, N2), 1);
        BOOST_CHECK_EQUAL(hf::rhf::HOMO_index(K2, N3), 2);
        BOOST_REQUIRE_THROW(hf::rhf::HOMO_index(K2, N4), std::invalid_argument);  // Cannot place more than 6 electrons in 3 orbitals
        BOOST_REQUIRE_THROW(hf::rhf::HOMO_index(K2, N_odd), std::invalid_argument);  // The unrestricted case is not supported
}


BOOST_AUTO_TEST_CASE ( lumo ) {

        size_t K1 = 4;
        size_t K2 = 3;

        size_t N1 = 2;
        size_t N2 = 4;
        size_t N3 = 6;
        size_t N4 = 8;
        size_t N_odd = 3;

        BOOST_CHECK_EQUAL(hf::rhf::LUMO_index(K1, N1), 1);
        BOOST_CHECK_EQUAL(hf::rhf::LUMO_index(K1, N2), 2);
        BOOST_CHECK_EQUAL(hf::rhf::LUMO_index(K1, N3), 3);
        BOOST_REQUIRE_THROW(hf::rhf::LUMO_index(K1, N4), std::invalid_argument);  // There is no lumo for 8 electrons in 4 spatial orbitals
        BOOST_REQUIRE_THROW(hf::rhf::LUMO_index(K1, N_odd), std::invalid_argument);  // The unrestricted case is not supported

        BOOST_CHECK_EQUAL(hf::rhf::LUMO_index(K2, N1), 1);
        BOOST_CHECK_EQUAL(hf::rhf::LUMO_index(K2, N2), 2);
        BOOST_REQUIRE_THROW(hf::rhf::LUMO_index(K2, N3), std::invalid_argument);  // There is no LUMO for 6 electrons in 3 spatial orbitals
        BOOST_REQUIRE_THROW(hf::rhf::LUMO_index(K2, N_odd), std::invalid_argument);  // The unrestricted case is not supported
}
