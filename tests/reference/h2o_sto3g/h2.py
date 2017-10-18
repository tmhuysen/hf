from horton import *
import numpy as np
np.set_printoptions(linewidth=150)

# Specify some data
water = IOData.from_file('../h2.xyz')
threshold = 10 ** (-6)
basis_name = "STO-3G"


# Create a Gaussian basis set
gobasis = get_gobasis(water.coordinates, water.numbers, basis_name)

# Compute the one- and two-electron integrals
S = gobasis.compute_overlap()
T = gobasis.compute_kinetic()
V = gobasis.compute_nuclear_attraction(water.coordinates, water.pseudo_numbers)
tei = gobasis.compute_electron_repulsion()

# Calculate H_core
H_core = T + V
print "H_core: \n {}".format(H_core)

# Create alpha orbitals
orb_alpha = Orbitals(gobasis.nbasis)
guess_core_hamiltonian(S, H_core, orb_alpha)
print "Initial guess for C: \n {}".format(orb_alpha.coeffs)

# Construct the restricted HF effective Hamiltonian
external = {'Nuclear repulsion': compute_nucnuc(water.coordinates, water.pseudo_numbers)}

terms = [
    RTwoIndexTerm(T, 'Kinetic energy'),
    RTwoIndexTerm(V, 'Nuclear attraction'),
    RDirectTerm(tei, 'HF coulomb'),
    RExchangeTerm(tei, 'HF exchange'),
]
hamiltonian = REffHam(terms, external)

# Decide how to occupy the orbitals (5 alpha electrons and 5 beta electrons)
occ_model = AufbauOccModel(1)

# Converge WFN with plain SCF
scf_solver = PlainSCFSolver(threshold)
scf_solver(hamiltonian, S, occ_model, orb_alpha)
