from horton import *

# Specify some data
water = IOData.from_file('../../../docs/h2o.xyz')
threshold = 10 ** (-6)
basis_name = "STO-3G"


# Create a Gaussian basis set
gobasis = get_gobasis(water.coordinates, water.numbers, 'STO-3G')


# Compute Gaussian integrals
olp = gobasis.compute_overlap()
kin = gobasis.compute_kinetic()
na = gobasis.compute_nuclear_attraction(water.coordinates, water.pseudo_numbers)
er = gobasis.compute_electron_repulsion()

# Create alpha orbitals
orb_alpha = Orbitals(gobasis.nbasis)

# Initial guess
one = kin + na
guess_core_hamiltonian(olp, one, orb_alpha)


print orb_alpha.coeffs

# Construct the restricted HF effective Hamiltonian
external = {'nn': compute_nucnuc(water.coordinates, water.pseudo_numbers)}

print external
terms = [
    RTwoIndexTerm(kin, 'kin'),
    RDirectTerm(er, 'hartree'),
    RExchangeTerm(er, 'x_hf'),
    RTwoIndexTerm(na, 'ne'),
]
ham = REffHam(terms, external)

# Decide how to occupy the orbitals (5 alpha electrons and 5 beta electrons)
occ_model = AufbauOccModel(5)

# Converge WFN with plain SCF
scf_solver = PlainSCFSolver(threshold)
scf_solver(ham, olp, occ_model, orb_alpha)
