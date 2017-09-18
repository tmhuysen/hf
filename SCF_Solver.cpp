#include <libint2.hpp>
#include <libint-eigen.hpp>

int main () {
    libint2::initialize();

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/n.xyz";
    std::string basis_name = "STO-3G";


    Molecule particle (xyzfilename);
    // 2. CALCULATIONS
    Basis basis (particle, basis_name);

    auto S = basis.compute_overlap_integrals();
    auto T = basis.compute_kinetic_integrals();
    auto V = basis.compute_nuclear_integrals();

    auto tei = basis.compute_two_electron_integrals();


    // Finalize libint2
    libint2::finalize();
    return 0;
}