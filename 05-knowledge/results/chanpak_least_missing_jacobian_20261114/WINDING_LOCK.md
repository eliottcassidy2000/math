# Relative winding and a hidden congruence lock on a physical apex chart

There is a sharp distinction between symmetry group size and the
interaction of two symmetric components with a reader. Here is a
physical instance using a pair of strong six-vertex cores with
the same number of Hamiltonian paths, the same full path-cover
polynomial, the same score sequence, and isomorphic automorphism
groups. One core traps every count in 3Z under arbitrary ordinal
one-apex continuations. The other realizes all three residue
classes already at the smallest chart, including prime 109.

## An exact gcd test for all apex cuts on a given core

For a tournament B define L_ij to be the sum of H(B−P) over
all even-order directed paths P starting at i and ending at j.
An apex that beats exactly U has augmented count
$$
 K(U)=H(B)+2D(U),\qquad D(U)=\sum_{i\in U,j\notin U}L_{ij}.
$$
Indeed, each odd cycle through the apex becomes one even-order
path of this type and has packing weight two. The remaining
packing has weight H(B−P). Packings avoiding the apex give
the baseline H(B).

Let r_i=Σ_j L_ij and s_ij=L_ij+L_ji for i<j. The greatest
common divisor of all apex half-increments D(U) is exactly
$$
             \gcd(\{r_i\},\{s_{ij}\}).
$$
One direction is the quadratic identity
D(U)=Σ_{i∈U}r_i−Σ_{i<j both in U}s_ij. The other follows by
taking singleton masks, then recovering s_ij as
D({i})+D({j})−D({i,j}). Thus all masks are locked modulo a
specified number exactly when their one- and two-contact
tests are locked. A gcd of one is only a test against a
global common divisor, not a claim about filling residues
or absolute values for every modulus in a particular family.

## Same symmetry group, different winding

Index two triangles by X_i,Y_i for i∈Z/3Z with indices
0,1,2,3,4,5. In both cores, X_i→X_(i+1),
Y_i→X_i, and X_i→Y_j for i≠j. In B_+ orient
Y_i→Y_(i+1), and in B_− orient Y_i→Y_(i−1).
The outgoing masks are

    B_+ = (50,44,25,17,34,12)
    B_− = (50,44,25,33,10,20).

They are strong: from X to Y use an unmatched edge and
from Y to X use a matched edge. Their score classes are
exactly X and Y (scores 3 and 2). Each automorphism must
preserve both, act on each by a rotation, and preserve
the matching, hence the automorphism group in both cases
is exactly simultaneous C3. Matching also makes the
relative direction of the two rings invariant under
isomorphism, so the two cores are not isomorphic.

Both have H=45 and P_B(2)=432 by the replay. Since the
ordered cover polynomial on six vertices is monic and has
only even degrees 6,4,2, these evaluations give exactly
(P_B(t)=t^6+16t^4+28t^2) in each case.

For B_− the even end/starting-path kernel, computed by
subset DP, is

    0 6 3 8 8 8
    3 0 6 8 8 8
    6 3 0 8 8 8
    7 4 4 0 3 6
    4 7 4 6 0 3
    4 4 7 3 6 0

The row sums are 33,33,33,24,24,24 and column sums are
24,24,24,33,33,33. Every pair sum is divisible by 3.
Thus K(U) is divisible by 3 for every U. For B_+, the
singletons X_0 and Y_0 give 109 and 95 respectively,
and the empty mask gives 45. These three residues are
distinct modulo 3. This is stronger than simply testing
that its gcd is one.

## The locked core stays locked after any ordinal growth

For B_−, the odd starting-path row totals are
39,39,39,33,33,33, and the odd ending-path column totals
are 33,33,33,39,39,39. In particular all of the following
are divisible by three for every U:

1. ordinary path count and two-cover count;
2. augmented apex count;
3. the unsigned and signed starting sums on U;
4. the unsigned and signed ending sums on the complement.

The endpoint numbers follow from the same path/complement DP.
To see the assertion is stable, use the
[ordinal transfer](../chanpak_selective_target_20261114/CUT_STATE.md):
appending an ordinal block gives the new path and two-cover
counts as products; each new boundary contribution is a
sum of products of one coordinate from **each** block;
the new augmented count is
(h_RK_L+h_L(K_R-h_R)+P_LE_R+S_LT_R), with E_R the
weighted end total and T_R its signed difference. The
ending-state recursion is the dual of the starting one.
Thus a marked core with all these coordinates multiples
of three poisons the entire one-apex ordinal chart even
if the left and right tournaments and their contacts are
chosen arbitrarily and are of unbounded order. This is a
scope restriction on a physical family, not on all ways
to glue tournaments.

The two cores differ by reversing all three arcs of Y.
That preserves degree sequences, automorphism structure,
and even the unrestricted cover polynomial, yet changes
whether the interface can ever reach a prime greater
than three. This is a concrete instance of **relative
phase**, not simply an analogy with symmetry size. The
natural next question is to find structural criteria for
which relative phases open a chart, and which such
charts cover prescribed integers rather than merely
admitting some primes.

Run `python replay_winding.py` for the literal kernel and
all mask checks.
