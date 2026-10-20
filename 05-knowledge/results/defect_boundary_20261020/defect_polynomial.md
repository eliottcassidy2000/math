# The ordered-defect polynomial of a fixed permutation

**PROVED:** the product formula, minimum and sparse-support bound below.
**FINITE-EXACT:** [defect_check.py](defect_check.py) and its recorded output.
The polynomial depends on the reference total order. It is **not** an
isomorphism-invariant weight for a Burnside sum.

Let a vertex permutation consist of consecutive cycles of odd lengths
`lambda=(l_1,...,l_t)`, each written in increasing cyclic order within its
block. Let `F_lambda(z)` count the tournaments fixed by this permutation,
weighted by their number of backward arcs in the displayed order. Then

\[
F_\lambda(z)=
\prod_i\prod_{d=1}^{(l_i-1)/2}(z^d+z^{l_i-d})
\prod_{i<j}(1+z^{\operatorname{lcm}(l_i,l_j)})^{\gcd(l_i,l_j)}.
\tag{1}
\]

Within a cycle of odd length `l`, the orbit of cyclic difference `d` contains
`l` pairs. Its two orientations have exactly `d` or `l-d` backward arcs.
Between two separate blocks, every cross pair has the same reference
direction. There are `gcd(a,b)` orbits of length `lcm(a,b)`, with costs zero
or the entire orbit length. Independent choices over unordered-pair orbits
give (1). An even vertex cycle instead fixes no tournament.

Three immediate consequences are

\[
F_\lambda(1)=2^{\,\sum_i(l_i-1)/2+\sum_{i<j}\gcd(l_i,l_j)},\quad
\min\deg F_\lambda=\sum_i(l_i^2-1)/8,\quad
z^{\binom n2}F_\lambda(z^{-1})=F_\lambda(z).
\]

The first is the inherited Davis fixed-permutation count. The second realizes
the universal lower bound proved in [sparse_theory.md](sparse_theory.md);
the third is tournament converse. This is a cost refinement of the same
edge-orbit mechanism, not an independent replacement for it.

For example,
`F_(3)(z)=z+z^2`,
`F_(5)(z)=(z+z^4)(z^2+z^3)`, and
`F_(3,3)(z)=(z+z^2)^2(1+z^3)^3`.
At `z=1`, the last gives 32, the first two-three-cycle rootless contribution
before the cycle-index denominator is applied.

## A minimal hostile against weighted Burnside

On four vertices, the consecutive permutation `(0 1 2)(3)` has polynomial
`z+z^2+z^4+z^5`. The conjugate permutation `(0 1 3)(2)`, relative to the
**same** reference order, has `z^2+2z^3+z^4`. A scalar conjugacy-class
substitution therefore cannot retain this ordered cost. The missing coordinate
is the placement of the cycles within the order.

The intrinsic minimum feedback statistic avoids this defect, but then counting
fixed tournaments of a given minimum is a different problem. The capped-core
compiler addresses its bounded-budget layers instead of applying an invalid
weighted Burnside formula.

Checks: exhaustive labelled orientations for every odd cycle type through
order six; signed pair-orbit walks for all 5,913 permutations of orders one
through seven, including even-cycle negative controls; all 2,035 odd cycle
types through total order 30 for the three identities; sharp triangle-join
controls through 30 triangles. The general proof, rather than these finite
checks, establishes the all-order statement.
