#define BOOST_TEST_MODULE "SCF_functions"

#include "SCF_functions.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( homo ) {

        unsigned K1 = 4;
        unsigned K2 = 3;

        unsigned N1 = 2;
        unsigned N2 = 4;
        unsigned N3 = 6;
        unsigned N4 = 8;
        unsigned N5 = 10;
        unsigned N_odd = 3;

        BOOST_CHECK_EQUAL(HF::HOMO_index(K1, N1), 0);
        BOOST_CHECK_EQUAL(HF::HOMO_index(K1, N2), 1);
        BOOST_CHECK_EQUAL(HF::HOMO_index(K1, N3), 2);
        BOOST_CHECK_EQUAL(HF::HOMO_index(K1, N4), 3);
        BOOST_REQUIRE_THROW(HF::HOMO_index(K1, N5), std::invalid_argument);  // Cannot place more than 8 electrons in 4 orbitals
        BOOST_REQUIRE_THROW(HF::HOMO_index(K1, N_odd), std::invalid_argument);  // The unrestricted case is not supported

        BOOST_CHECK_EQUAL(HF::HOMO_index(K2, N1), 0);
        BOOST_CHECK_EQUAL(HF::HOMO_index(K2, N2), 1);
        BOOST_CHECK_EQUAL(HF::HOMO_index(K2, N3), 2);
        BOOST_REQUIRE_THROW(HF::HOMO_index(K2, N4), std::invalid_argument);  // Cannot place more than 6 electrons in 3 orbitals
        BOOST_REQUIRE_THROW(HF::HOMO_index(K2, N_odd), std::invalid_argument);  // The unrestricted case is not supported
}


BOOST_AUTO_TEST_CASE ( lumo ) {

        unsigned K1 = 4;
        unsigned K2 = 3;

        unsigned N1 = 2;
        unsigned N2 = 4;
        unsigned N3 = 6;
        unsigned N4 = 8;
        unsigned N5 = 10;
        unsigned N_odd = 3;

        BOOST_CHECK_EQUAL(HF::LUMO_index(K1, N1), 1);
        BOOST_CHECK_EQUAL(HF::LUMO_index(K1, N2), 2);
        BOOST_CHECK_EQUAL(HF::LUMO_index(K1, N3), 3);
        BOOST_REQUIRE_THROW(HF::LUMO_index(K1, N4), std::invalid_argument);  // There is no lumo for 8 electrons in 4 spatial orbitals
        BOOST_REQUIRE_THROW(HF::LUMO_index(K1, N_odd), std::invalid_argument);  // The unrestricted case is not supported

        BOOST_CHECK_EQUAL(HF::LUMO_index(K2, N1), 1);
        BOOST_CHECK_EQUAL(HF::LUMO_index(K2, N2), 2);
        BOOST_REQUIRE_THROW(HF::LUMO_index(K2, N3), std::invalid_argument);  // There is no LUMO for 6 electrons in 3 spatial orbitals
        BOOST_REQUIRE_THROW(HF::LUMO_index(K2, N_odd), std::invalid_argument);  // The unrestricted case is not supported
}
