# One enumeration, several growth coordinates

**PROVED:** the abstract response lemma below. The accompanying four grades
organize separately proved counting results; their combination is not an
unproved assertion of a complete tournament-isomorphism carrier.

## 1. The common mechanism: conserve a grade, then quotient its response

Suppose a partial construction has nonzero weight w and response h in R^r.
Its allowed future multiplicity vectors form a nonempty finite set C of
nonnegative integer vectors. Continuing by k in C contributes

    w * q^(h dot k) * f(k),   q>1, f(k)>0.

Two responses h,h' are equivalent up to one absorbable weight factor iff

    h'-h belongs to span_R(C-C)^perp.                         (1)

Proof: their ratios are q^((h'-h) dot k), constant on C precisely when
(h'-h) dot (k-k')=0 for every pair k,k'. The resulting constant can be
absorbed into w. QED.

If all feasible continuations have the same grades Ak=b, every vector in
the row space of A is therefore removable from the response. Equality with
the full removable space requires the additional spanning guard
span(C-C)=ker A; it must not be assumed just from Ak=b.
For example a singleton feasible set makes every response removable even
when the displayed grade has a nontrivial kernel.

For tournament cycle insertion the conserved grade is total vertex mass,
k_1+3k_3+5k_5+...=n. Unit cycles supply the necessary test continuations,
so h_j-jh_1 is the exact projective quotient. This is the proved smaller
recurrence in [PROJECTIVE_RESPONSE.md](recurrence/PROJECTIVE_RESPONSE.md).
It answers the 'bit for free' question in a precise way: a coordinate that
acts uniformly on every allowed future contributes only a scalar weight.
It does not mean named port alignments are dispensable for a different test.

## 2. Four grades, four different kinds of partial structure

### Vertex mass: exact coefficient recursion

The cycle-insertion state remembers future gcd responses, modulo the
conserved-mass direction. Remaining mass shrinks, so responses at cycle
lengths exceeding that mass disappear. This counts all tournament classes
without retaining individual representatives. It does not recover the
metagraph adjacency of those representatives.

### Module surplus: composition refinement

For a prime root of size k expanded only to k+h vertices, the relevant
cycle-index data is its weighted Taylor jet through h. Increasing h adds
moments indexed by odd partitions of the new degree. The child-palette
observer can identify further jet directions; independent formal palettes
and actual bounded tournament palettes are different observers.
See [base jets](composition/RESULT.md) and
[observable quotients](composition/QUOTIENT_EXTENSION.md).
Nonuniform attachment inside an old module leaves this grammar and requires
its incidence data. Prime-count input cannot silently be omitted.

### Missing comparisons: bounded even symmetry

A d-hole partial tournament allows even-cycle symmetry supported on at most
2d vertices in each fixed-permutation contribution. Its automorphism group's
2-primary order is at most 2^d. The budget is sharp, but the missing-edge
shape matters: a star with at least two holes fixes its center and forces
odd automorphism order. See [defect grammar](defects/MISSING_COMPARISON_GRAMMAR.md)
and [star exposure](partials/PARTIAL_STAR_CLASSIFICATION.md).
Filling a hole and forgetting a mark are distinct operations. Their exact
incidence counts need stabilizer/flag information; the C3 exposure hostile
and the four-vertex unmarked reversal loop are explicit tests.

### Symmetry cycle deficit: analytic bands

Strip fixed vertices from an automorphism and let its remaining cycle core
have m vertices and c cycles. The deficit m-c=2k controls how each additional
fixed vertex changes its normalized weight. Summing those fixed vertices
gives one rational pole bank at 4^k. The core types in that bank are indexed
by ordinary partitions of k, with maximal support3k attained by3^k.
Thus pole order3k+1 is a structural growth invariant, not numerology.
See [global meromorphic expansion](limits/ALL_POLES_CORE_EXPANSION.md).
The exact raw coefficients are recovered by multiplying normalized counts
by 2^(n choose2)/n!; approximation alone does not justify integer rounding.

## 3. What a full metagraph theory still owes

These grades are exact for their declared operations. They do not yet supply
a finite closed classification of every partial attachment under arbitrary
vertex insertion, re-marking, and edge reversal. A proposed common state
must pass three tests:

1. It must distinguish histories that a permitted future operation separates.
2. Its transition weights must include stabilizers, multiplicity and any
   forgotten-root or forgotten-exposure correction.
3. Its observer must be stated: one coefficient, all palettes through a
   surplus, named completions, or actual isomorphism-level neighbors.

The present session supplies new exact quotients and classifications for
these tests, not a proof that they all collapse into one small scalar state.
The next precise frontier is to combine module splitting under a new vertex
with the proven surplus observer, while retaining just enough aligned ports
to determine when the old prime decomposition ceases to apply.
