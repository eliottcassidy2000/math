# The carry word is a weighted cut graph

**Scope.** This section concerns the actual fixed-core two-tape formula
developed in the earlier [absolute steering packet](../chanpak_absolute_steering_20261115/TWO_TAPE_STEERING.md).
It converts the next inverse problem into an exact weighted binary cut.
It is not an equality or universal placement theorem.

## Fixed carry means a fixed cut geometry

Set q(x)=p(x)/2. Write x with r binary digits b_i from most significant
to least; leading zeros are allowed. By the recursion for p,

```text
q(x) = sum_{i<j} w_{j-i} b_i (1-b_j),
w_1=1,             w_d=2^(d-2) for d>=2.
```

Indeed appending a zero adds floor(prefix/2) plus the last prefix bit
(if any). Appending a one adds zero. The lowest old bit contributes
one; each earlier bit at distance d>=2 contributes 2^(d-2).

Now fix the sum of two r-bit words x,y and their binary carry word
(c_j is the incoming carry at the j-th position from the right).
The bit total a_j at that position is
`a_j = bit_j(x+y)+2 c_(j+1)-c_j`. Values 0 and 2 are fixed to 00
or 11; value 1 is free to choose 10 or 01. Conversely every assignment
of these free choices gives the fixed sum and carry, by column
addition. Let I be the positions with bit total 1. Then

```text
q(x)+q(y) = C_carry + sum_{i<j in I} w_{j-i} (b_i XOR b_j).
```

The identity for a free pair follows by adding the two complementary
terms. For fixed/fixed positions the contribution is deterministic,
and for fixed/free positions the two entries of a free pair sum to one,
so the contribution is again constant. The exact constant C_carry
is the checksum with all free positions assigned the same orientation.
Indices in this formula may be taken left-to-right or the reverse;
the nonconstant part uses only their distance.

The target now has two separate types of choice: choose a carry word
to choose the induced weighted complete graph on I, then pick an
orientation of its vertices to choose a cut. First choosing a full
dyadic residue bank is not necessary to ask the inverse question.

For the earlier two-tape shell, the high words of x and y can coincide.
For x=Q0*u+v, y=Q0*u+w, Q0=2^r, the exact split of p has cross terms
depending on rev(v)+rev(w), and boundary terms depending on the two
top low digits. Both totals are constant in the carry cell because
both are weighted sums of the a_j. Thus the only varying surcharge is
the cut formula on the low words. The smooth term remains the
quadratic product of the two factors; on a fixed sum this depends
on their difference. In particular, global complement of the free
assignments leaves the surcharge invariant, but generally does not
fix the smooth scalar when the two linear coefficients differ.

## Exact small-response exclusion in the all-complement shell

If every low bit is free, the minimum nonzero surcharge q(x)+q(y)
is as follows:

```text
r=2t+1, t>=2:  2^t
r=2t, t>=3:    3*2^(t-2)
```

These are the weighted degrees of a single central vertex (either
one if there are two). Indeed the sum of the weights at distances
1 through k is 2^(k-1) for k>=1. The minimum degree is therefore
the stated central degree m_r. If the first and last positions have
different colors, their mutual edge of weight 2^(r-3) alone gives at
least m_r for r>=5. If their colors agree, let A be the other color.
A one-element A has cut cost its degree. If |A|>=2, every vertex in A
has at least m_r/2 in edges to the two endpoints: each internal
distance weight w_d is at least one half of the sum of the weights
from 1 to d. Thus again the cut is at least m_r.

This gap grows as sqrt(Q0). It constrains the *isolated checksum*
at a frozen all-complement shell; it does not exclude those path
counts when the product and linear terms vary or the carry changes.

## A linear-length carry-selected differential (not an absolute phase fix)

There is a constructive reason not to insist on all bits free. Put an
anchor free position at an index h>=2a+3 from the right. Given any
integer R in [0,2^a), make the positions at distances j+2 below it free
precisely for the 1 bits j of R. Make all other low positions 00.
The two assignments "all free positions oriented alike" and "only
the anchor reversed" have the same exact sum and binary carry, and
their surcharge difference is exactly R, since those cut edges have
weights 2^j. All chosen positions are at indices >=a+2. Therefore
for the physical fixed-core formula the polynomial part has the same
residue mod 2^(a+1) at both assignments, and their count difference is
90R mod 2^(a+1). This uses only O(a) low positions, as opposed to
constructing a single fixed aligned block serving every syndrome.

The catch is essential: the baseline of that chosen cell still depends
on R. This gives any desired **relative** correction with an adaptive
carry, not an arbitrary absolute residue around a frozen baseline.
The two numerical counts can also differ greatly in absolute value.
A possible next inverse step is to force a small upper/lower product
envelope while selecting the free-position graph; on its own the
differential does not solve this.

## A finite local feedback interpretation

At dyadic precision 2^a, edges at distance >=a+2 vanish. In a fixed
carry cell, scan bit positions high-to-low. Keep the last a+1
assignment bits and remaining cut debt mod 2^a. Exact backward
reachable sets decide whether a given debt can be completed;
recursively choose the first viable bit. This makes the choices a
directed graph with one prescribed outgoing edge from every viable
nonterminal state, and length drops at every move. For the actual
count residue mod 2^(a+1), one must also keep or enumerate the bottom
a+1 bits of both addresses for the product and linear parts; these
are local and can be included in the terminal condition. The debt
states are finite for fixed precision. No claim is made that the
initial state is always viable or that the matching word is within
an absolute count window.

This is the precise use of the Collatz analogy: state determines one
next edge, but a distance or termination proof must be separate.
Here termination comes from remaining tape length rather than an
unproved property of the Collatz graph.

## One exact carry decision at the first output bit

The carrier coefficients satisfy Q=0 mod 4 and T=S=2 mod 4.
Consequently at fixed exact sum, only the parity of the cut can change
F mod 4. Weights w_d are odd precisely for d=1,2, so write G_I
for the graph of free positions with an edge at distance 1 or 2. The
parity of a cut is the sum of degrees at its selected vertices mod 2.
Thus the fixed carry cell supplies both odd classes mod 4 if some vertex
of G_I has odd degree, and supplies just its baseline class otherwise.

**Classification.** Each component of G_I is Eulerian if and only if it
is an isolated position or a run of three consecutive positions. Indeed,
a distance-two gap between consecutive free positions leaves exactly one
edge across that gap, a bridge that cannot occur in an Eulerian component.
Thus a nonsingleton Eulerian component must be a consecutive run.
A two-vertex run has odd degrees, a three-vertex run is a triangle, and
every longer run has degree three at its second vertex. The parity
classifier is therefore a local, radius-two obstruction to deciding a
carry without remembering which positions are free.

For the differential construction above, the two addresses may have the
same fixed high word u. Set A=2^h and
L=2^(h-a-1) rev_a(R). The two low pairs are (A+L,0) and (L,A).
Their exact count difference, from the first to the second, is
`Q*A*L + (S-T)*A + 90R`. This follows from the fixed-sum
product and linear terms and the cut formula. It is independent of u;
the same carry-selected differential occurs at all common heights.

## Sources and what the analogies actually buy

**Quasi-Riemann.** The [11/12 alternate proof, §2.2](https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex)
does not infer a hard individual row directly from a mean-square bound.
It first supplies many near-replicas of that row at sixth powers of
primes; their only differences come from the divisibility exclusions.
The [7/8 version](https://github.com/openai/math/tree/main/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026)
is the stronger stated zero-free half-plane. The repo used it for the
[idoneal/class-group exceptional-set comparison](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/binary_information_class_groups_largest_anomalies_20261017.md),
with a separate local Paley ranking detector. The transfer here is
methodological: assign a hard residue to concrete invariant copies, not
to an unconditioned average of carry cells. One exact replication is
visible in the cut law. At precision 2^a, an arbitrary additional group
of free positions separated by distance >=a+2 from the existing free
positions contributes a new cut without interacting with the old cut
modulo 2^a. The constant C_carry must still be recalculated; it cannot
be erased from the physical target equation. There is currently no
uniform real-window estimate with enough such copies to force equality,
and neither Quasi-Riemann result implies one.

**Planar Jacobian.** [El Hilany and Tsigaridas, §§3–5](https://arxiv.org/html/2101.05245)
decompose the non-properness set of planar polynomial maps by boundary
faces related to their Newton polytopes. In a different counting problem,
the repo's [Boolean Newton carry boundary](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-2502-endpoint-boolean-newton-carry-tournament-and-dipole-boundary.md)
exhibits outer projections that lose the actual carry and a top atom.
For the two-tape inverse, all genuine exact positive-integer fibers are
bounded: since p(x),p(y)>=0, a target N bounds x<=(N-B)/T,
y<=(N-B)/S. This is a discrete properness test, not an application of
the Jacobian conjecture. It is exactly what is lost when passing to
the modular charts. For example 2857 has no witness in this one-core
family: x=y=0 gives 2855 and any nonzero address gives at least 6769.
Yet the [earlier unrestricted modular compiler](../chanpak_absolute_steering_20261115/PREFIX_PROGRAMMING.md)
supplies congruent counts for 2857 at all fixed dyadic precisions.
Those witnesses must escape every fixed coordinate bound, since a
bounded subsequence would repeat one pair and give exact equality.
This is a failure of local-to-global reasoning in a literal family,
not an obstruction to other tournament constructions at 2857.

An effective real-side branch test is to retain both the interval for
the factor difference and the cut debt. For fixed sum P, let z=x-y,
so xy=(P^2-z^2)/4. On any partial digit assignment take a containing
interval [z_l,z_u] for z; the polynomial smooth part has exact extrema
on that real interval, since it is a concave quadratic. A lower bound
for the cut is the weight of already assigned cut edges; an upper
bound adds every remaining edge weight. These give a safe enclosing
count interval. Any branch excluding N can be discarded *before*
modular reachability. This substitutes an explicit digit-cell
boundary decomposition for a claim that a modular lift is proper.

## Replay

`python replay_carry_cut.py` (also `python -O replay_carry_cut.py`).
The replay checks the cut formula and mod-4 graph criterion on literal
small addresses; verifies the complementary min-cut values; checks
the adaptive differential at assorted common high words; and checks
the canonical feedback selector against brute-force viable carry cells.
The infinite conclusions follow from the identities, not these tests.
