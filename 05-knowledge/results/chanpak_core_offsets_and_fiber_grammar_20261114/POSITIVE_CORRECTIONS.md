# Positive finite-inventory corrections

This is a window corollary of the inherited SCC multiplicity law,
not a new scalar Burnside formula. Use formal H labels [a][b]=[ab].
Separate the strong atom series into B of nontrivial classes of size
at most M, and U of the remaining strong classes. Include singletons
via z, and put G=(1-z-B)^(-1). The exact ordinary class series is

    V = sum_{j>=0} U^j G^(j+1).

A term has exactly j omitted-size SCCs and retained words in the
j+1 gaps; all contributions are nonnegative. If Q is the exact
minimum strong Hamiltonian count at order M+1 (M>=2), Busch's
minimum does not decrease with order. Thus all j>r terms vanish
whenever h<Q^(r+1) OR n<(M+1)(r+1). Equality at the corner really
can contribute, by concatenating r+1 copies of an extremal strong
(M+1)-tournament. See
[Busch](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).
This is sharp for the generic two-cut window, not an assertion
that individual pairs outside it always have a missing term.

For h<Q² the only correction to the retained catalog is linear:

    V_h = G_h + sum_{qde=h} U_q G_d G_e,

where the subscript extracts the H label (preserving z).
Adding a newly catalogued batch C to G updates it exactly by
sum_{j>=0} C^j G^(j+1). If the direct atom catalog is complete only
for certain sizes/counts, a witness per count is not a substitute
for its true class multiplicities.

For M=8, Q=75; all odd fibers through 73 are exact for every
order from that strong catalog alone. For M=12, Q=625. A fixed
finite strong catalog cannot account for A000568 at unbounded n:
its word series has at most exponential growth, whereas the number
of ordinary classes is at least 2^(n(n-1)/2)/n!. The full recurrence
is exact, not a claim to beat the unrestricted strong census cost.
Related inherited work: [fixed-fiber SCC laws](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md)
and [divisor quotient and Burnside](../session239_sequence_engines_20261110/proof_bridges/PROOF.md).
