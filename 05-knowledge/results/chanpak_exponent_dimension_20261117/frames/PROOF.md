# A closed Hamiltonian frame has an exact positive mean surplus

**Proved for all orders n >= 3; literal finite replay through 7.** This is a
first-moment law for lines, not an individual-line inequality. The fixed
Hamilton path setting and closed/open complement pairing are inherited
from THM-643 and THM-4729. The full-cube forest expansion is inherited from
THM-076/848. THM-589 already gives the mean of the all-tilings ensemble
through the second moment of the ordinary random tournament.

## Two adjacent strata in the same cube

Orient 1 -> ... -> n; all other signs are independent fair bits.
A tiling is closed if n -> 1 and open if 1 -> n. The nonpath complement
matches each closed tiling to its open partner. Write m=n-2r,
mu=n!/2^(n-1), and

    A(m,r) = sum_{c=1}^r 2^c binom(r-1,c-1) binom(m,c)
           = [z^r] ((1+z)/(1-z))^m       (r>=1).

Mean Hamilton-path counts for all frames (M), for closed frames (C), and
for open frames (O) are

    M/mu = 1+sum_{1<=r< n/2} A(n-2r,r)/(n)_{2r},
    C/mu = 1+sum_{1<=r< n/2} (n/(n-2r)) A(n-2r,r)/(n)_{2r},
    O = 2M-C.

Consequently the exact, strictly positive signed count drift averaged
over *literal* closed-to-open lines is O-C, where

    (C-O)/mu = sum_{1<=r<n/2} 4r A(n-2r,r)/((n-2r)(n)_{2r})
             = 8/[n(n-1)] + O(n^-3).

The expectation uses line multiplicity on the merged metagraph, not a
uniform draw among its class vertices or among distinct neighbors. A blue
line ends once in each color; that fact does not entitle one to apply the
same mean conditioned on blue. There is no pointwise monotonicity: already
one of the 32 size-five lines has smaller H at its closed end.

## Proof

For a permutation, its path indicator is a product of n-1 factors
(1+x_ab)/2 with x_ab the signed arc, antisymmetric in a,b. In the sum over
permutations, the coefficient of an arc set vanishes unless it is a
linear forest with each nontrivial component even in edge length. In
fact, contract each component to a block. Odd length forces cancellation
between the two traversals. Even length gives twice the coherent traversal
sign. With k selected edges and c nontrivial components the coefficient
in units of mu is the traversal sign times 2^c/(n)_k. This recovers the
forest law in THM-076/848.

Conditional expectation kills each monomial using a still-fair edge. On a
coherent fixed path or cycle every surviving monomial has positive sign.
With k=2r fixed edges, choose c positive even run lengths in
binom(r-1,c-1) ways. On the noncyclic path the gaps place these runs in
binom(n-k,c) ways. On the circle the c marked runs are separated by
positive gaps of total length n-k, giving
(n/c) binom(n-k-1,c-1) = n/(n-k) binom(n-k,c)
placements. We never count the full cycle as a forest. This proves
the formulas. Reversing only the closing edge instead picks the other
half of the noncyclic family, hence O=2M-C.

For the error bound, when c<=r, (m)_c/m <= (m+c-1)!/m! for m>=1
(and the term is zero if c>m). Thus every term of A(m,r)/(m (n)_{2r})
is at most 2^c binom(r-1,c-1)/[c! (n)_{r+1}]. The r>=2 tail is
bounded by sum 8r 3^(r-1)/(n)_{r+1}, O(n^-3) since
(n)_{r+1} >= (n/2)^(r+1). The r=1 term is exactly 8/(n(n-1)).

As another formulation in the full random tournament, let Y denote the
number of directed Hamilton cycles counted up to rotation. Index-symmetry
and planting give E(H^2)/E(H)=M and E(HY)/E(Y)=C; hence
Cov(H,Y)/(E(H)E(Y)) = C/mu-1, and its excess over
Var(H)/E(H)^2 is (C-O)/(2mu). In particular the extra diagonal which
closes the frame has a response of order n^-2, while conditioning on a
whole Hamiltonian frame has order n^-1. These are relative responses.

## Exponent, dimension and physical choice

The two appearances of powers of two are different: 2^c counts the two
ways to traverse each surviving even component *before* reading its
orientation product, whereas the number of tilings is 2^m for m
independent free arc positions (here m=binom(n-1,2), not the m=n-2r
in the coefficient formula). Conditional averaging projects to the
forests supported on the fixed positions; the circle multiplies each
placement rank by an explicit boundary factor. This is a concrete
instance of the archive's Boolean support tomography (THM-2383); it
does not identify the full law from a scalar or certify any prescribed
odd outcome. The 20261005 algorithmic Collatz measure makes the analogous
distinction between uniform address mass and mass that remembers an
actual source integer.
