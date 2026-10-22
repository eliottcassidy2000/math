# Same spectrum and blue count, different metagraph black degree

Status: FINITE-EXACT explicit counterexample, proving non-determination.

Use lexicographic unordered pairs(i,j), i<j, with bit1 meaning i→j and0
meaning j→i. The almost-regular order6 tournaments

    A: 000110001000000
    B: 000110011000000

have the same skew characteristic polynomial

    lambda^6+15lambda^4+15lambda^2+1

and hence the same rho=66637/27755. Both are self-converse. Yet their complete
Hamiltonian/automorphism counts are

    A: H=41, |Aut|=1, h=H/|Aut|=41, blue b=9, deficit h-b=32;
    B: H=45, |Aut|=3, h=H/|Aut|=15, blue b=9, deficit h-b=6.

Thus their fundamental black degrees(h-b)/2 are16 and3. B is the known
order6 deficit-six center C3[K2] (its three two-vertex modules can also be
recovered directly; this identification is not needed for the counterexample).

All quantities are checked by enumerating the720 vertex permutations separately
for each explicit graph. Blue labelled paths are those for which reversal of
path positions induces an anti-automorphism; divide by the free Aut action.
The labelled blue counts are9 and27 respectively. Every automorphism fixing
an ordered Hamilton path fixes all vertices, explaining why the division is
exact. Characteristic polynomials and rho are obtained by exact symbolic
matrix arithmetic, not rounded spectral comparisons.

Consequence: even the pair (skew spectrum, number of blue tilings) does not
determine total Hamiltonian count, number of Hamilton-word orbits, or fundamental
black degree, already within strong self-converse almost-regular tournaments.
No scalar spectral normalization can completely classify metagraph wiring.
This is an obstruction to an exact classifier, not a refutation of a suitably
qualified asymptotic approximation with additional endpoint/score invariants.

Reproduce: python3 check_isospectral_metagraph.py.
