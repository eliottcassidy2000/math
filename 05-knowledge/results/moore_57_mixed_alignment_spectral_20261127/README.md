# The Moore spectrum forbids mixed equal-or-converse alignment on nine branches

2026-11-27. This strengthens the [outside-route obstruction](../moore_57_deficient_incidence_external_routes_20261127/README.md):
the sign of each aligned edge can now be either equal or converse,
independently. No score bound is needed for the main obstruction.
There is no change to the unrestricted entrywise window
160/162/164/166 and no assertion about Moore graph existence.

## A proper orbit is a test vector

Let A be the adjacency of a hypothetical degree-57, diameter-two
Moore graph. The 57 root branches each have 56 slots in the second
shell, with no internal slot edges and a perfect matching to every
other branch. The Moore identity A²+A=56I+J gives eigenvalues
7 and -8 on the orthogonal complement of the all-ones vector.

Call a pair exactly aligned if its matching transports its whole
interior tournament to the other interior or its converse. Take
any connected graph of such pairs and transport along its paths.
The closed-walk permutations at a reference branch form a group
G whose elements preserve or reverse the reference tournament.
The choice is unique for a particular permutation. The preserving
kernel has odd order (an involution would flip an oriented pair);
the index of the kernel is at most two. No G-orbit can have 56
slots, or even 28, because orbit sizes divide |G|. Transport an
orbit to each branch in the component. The resulting sets O_i
are path-independent, with a common size s, and respected by all
aligned pair matchings.

Use the centered indicator 1_{O_i} - s/56 at slots of each chosen
branch and zero on all other vertices, including the root and
its neighbors. The resulting vector is orthogonal to constants.
On a clique of k aligned pairs its Rayleigh quotient for A is
exactly k-1, since every matched pair of centered vectors agrees
and each branch has zero internal slot edges. Hence **k<=8**.
In particular the entire equal-or-converse template cannot exist
on any of the four Paley defect cliques, even if one chooses which
pairs are equal and which are converse. The previous conditional
count of equal blocks did not guarantee feasibility of that template.

For an aligned K8 the vector attains 7, so it is a 7-eigenvector.
At each outside slot the eight neighbors in the K8 therefore
include exactly s/7 marked slots. Every orbit size is divisible
by seven. In the same alignment component, an outside branch
can have at most s/7 aligned links into the K8, by testing an
outside slot in the transported orbit. Use the smallest orbit:
s<=27 gives at most three links. If G contains a reverser, every
orbit splits into equal odd score-sign halves (a reverser negates
all signed row scores, which are odd). Thus orbit sizes are
2 mod 4, so the smallest is at most 26 and its only multiple of
seven congruent to 2 mod 4 is 14. In this case the limit is two links.

## A count of the necessarily unaligned pairs

On a connected alignment component of k chosen branches, take
its smallest orbit. Each pair of centered indicators across an
unaligned matching has normalized inner product at least
-s/(56-s) (its uncentered overlap is nonnegative). An aligned
matching gives +1. With m unaligned pairs the Rayleigh bound gives

    m >= max(0, ceil(k(k-8)(56-s)/112)).

Always s<=27. If G has a reverser, s<=26. If it does not, all
orbit sizes are odd. With |row score|<=3 there must be at least
four orbits: suppose there were two of odd sizes a,b, a+b=56.
Transitivity on each orbit gives zero signed score internally
and a constant cross-score x or y, respectively. These are odd
and ax=-by. Their maximum magnitude is at least
max(a,b)/gcd(a,b)>=29/7>3. The gcd is an odd divisor of 56.
Thus in the positive-only case the score cap gives s<=13,
and s<=26 works for both cases combined.

A stronger count in the positive-only case uses all r orbit
indicators simultaneously. Normalize the r indicators by
square roots of their sizes and restrict their span to vectors
orthogonal to the branch constant. The matrix across a matched
pair in these coordinates is the normalized orbit-overlap matrix.
It has nonnegative diagonal and is the identity for aligned
pairs. The constant direction is fixed by its action and its
transpose, so on the remaining r-1 coordinates the trace is at
least -1, or exactly r-1 on aligned pairs. Sum the Moore Rayleigh
bound over a common orthonormal basis for these r-1 directions,
using the same coefficient vector in every branch. This gives

    2[(choose(k,2)-m)(r-1)-m] <= 7k(r-1),
    m >= max(0, ceil((r-1) k(k-8)/(2r))).

In particular r>=4 in a score-capped component with no reversing
closed walk, so its connected count is at least the r=4 bound.

If the alignment graph is disconnected, every pair across two
unions of components is unaligned. For each applicable connected
bound C(k), use the recurrence

    B(1)=0,
    B(k)=min(C(k), min_{1<=a<k} [a(k-a)+B(a)+B(k-a)]).

For sizes 14,14,14,15 this forces the following totals:

| Assumptions | Each 14 | The 15 | Total |
|---|---:|---:|---:|
| Arbitrary interiors | 22 | 28 | 94 |
| Absolute row scores at most 3 | 23 | 29 | 98 |
| Score cap and no reversing closed walk within each aligned component | 32 | 40 | 136 |

These count pairs that are neither equal nor converse under the
actual matching. They do not yet force any particular entry in
an unaligned block to reach 164. To turn this into a bound on the
matrix entries, the natural next step is to charge large deficits
in the normalized overlap to the individual arc inequalities at
caps 160 and 162. The present argument does not make that step.

## Small positive control

For degree seven the nonprincipal spectral upper bound is 2.
In the Hoffman–Singleton graph each of the 35 triangles of
root branches has closed transport a product of three disjoint
transpositions. Make a transitive reference tournament whose
order is reversed by that permutation; transport it to the two
other branches. All three pairs then align (one sign may be
converse). Each transported two-slot orbit gives a Rayleigh
quotient of 2, attaining the bound. The replay checks the Moore
identity, transports and exact eigenvector equation, as well as
the split arithmetic above. Run node replay.js.
