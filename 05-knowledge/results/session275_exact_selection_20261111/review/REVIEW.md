# Independent focused review of Session275

**Final status: PASS WITH MINOR CORRECTIONS (applied).** The review is scoped to new bridges, not a re-audit
of the repository's inherited theorem library. Chan–Pak remains OPEN.

## Landscape

The fixed-path bit encoding enumerates all framed objects. The canonicalization
is sound: ordered degree blocks plus all within-block permutations give a
complete canonical representative; equality of upper-triangle encodings also
fixes the degree data. Joining equal canonical codes implements actual
reframing of an isomorphic tournament, not accidental identification by H.
All cut masks are generated, and simultaneous two-arc flips require only their
endpoints to satisfy the scalar filter. The component interval lemma is the
ordinary integer intermediate-value argument and remains valid for nonmonotone
walks. Padding supplies size, not the missing unbounded-count hypothesis.

Independent verifier enumerates permutations of vertices directly to count
paths, canonicalizes over ALL permutations rather than degree-block reduction,
and computes components by explicit adjacency/BFS rather than union-find.
Orders3,4,5 reproduce the reported components and class counts. The literal
six-vertex H=23 seed is checked by permutations. At orders6,7 the producer
algorithm was inspected and the stored complete-spectrum runs compared with
component intervals; this review did NOT independently rerun those full
orders. Their finite exact claim rests on the exhaustive producer run.

## One-way switches

The square-zero identity is valid against every upper-triangular spacer:
D_i=e_1 r_i, r_i e_1=0, and C e_1=C_11 e_1, hence D_i C D_j=0.
Product expansion therefore eliminates every mixed term, with the entire
physical boundary retained. The identity is mathematical algebra, not a
proposal to implement subtraction as a tournament operation.

Independent permutation enumeration of every induced subset reproduces all
seven rooted packet entries of both literal pairs. Independent endpoint-recursive
DP checks four eleven-vertex pair composites and four spacer cubes, without
using producer carrier helpers. The explicit weight formula retains r=0
correctly; baseline and weights at m=1,2 agree. The gcd4 argument, greedy
superincreasing inverse, complete-sequence criterion and unit-deletion binary
incidence reconstruction are valid in the stated scopes.

Minor proof repair requested and applied: strengthen the displayed inequality
for r≥2 to g_r>64g_(r−1), rather than merely>(64/5)g_(r−1). The stronger
inequality follows from the displayed recurrence and is what yields the stated
ratio of successive weights. No formula or output changed.

## Lift-sensitive descent

The normalized-cofactor21 classification follows from the inherited graft law
and bounds. For n≥5, h+q≥7h forces h=1 after the divisibility-by4 control;
small n have the listed complete packets. Physical replacements dispose of
q=4,8,16 without treating arbitrary normalization as realizable. At q=12,20,
D≥3 is excluded by the inherited nonsink endpoint bound, then f=0 forces
q=20 and the source's penultimate vertex. The converse and the order≤19
bound follow immediately. The ordinal-product formula preserves both H and Q
and the strongness argument uses the actual cyclic quotient.

Independent permutation counting gives H(B)=29, Q(B)=580, and H(graft)=609.
Minor definition clarification requested and applied: the first rooted path
in U is nonempty, and the second is allowed empty. This is the convention
needed for U=2^t. The mathematical formula was already using that convention.
The classification relies on inherited ordinary exclusion H=21 and proper-graft
inequalities, which this review did not separately reprove.

## Explicit first corridor stabilization

The three moves in PORTAL63.json / BRIDGE63.json are independently checked
by permutation counting, including actual flipped-arc sets, canonical relabel
maps, and common Hamilton paths. A common path ensures the flipped arcs are
free in a lawful chart, not merely arbitrary arcs of an unframed tournament.
The lower path goes63→61→59 and reaches a source-padded order-seven graph;
the upper move goes63→65 and reaches another source-padded graph. The source
vertices are7 and6 respectively. Together with the finite n7 component input,
this genuinely stitches23..61 and65..105 through63 in the variable-order
corridor. This is a finite structural stabilization, not its all-order recurrence.

## Two-sided seams and the moving-root bridge

The inherited one-sided recurrences expand to the stated bifan polynomial
with all zero-fan boundaries kept separate. The four-measurement formulas
recover the polynomial and its determinant invariant correctly. Completing
the product square gives exactly the shifted divisor equation, with positive
factors from the coefficient bounds. The finite inverse is therefore both
necessary and sufficient for its fixed seed. The sparse-family bound is scoped
to that family and a fixed finite bank.

Independent path enumeration reproduces all seven explicit root packets in the
note. Independent endpoint-recursive DP counts their four(1/2,1/2) bifans.
The balanced-fan pair has constant gap2 because its leading and middle
coefficients agree; this is not assumed to connect to other count corridors.
The cut at the two old roots flips the intended fan incidences, and repairing
ten old arcs restores T. The note correctly does not claim that the repair
steps preserve a frame or scalar corridor.

The subsequent balanced-seam translation theorem is also valid: for fixed U,
substituting c=R−L and z=e−s into the bifan polynomial gives the stated
common additive shift. A balanced port has zero contributions to both charges;
converse-invariant orbit unions and the middle of an odd transitive core
supply actual examples. The61/63 portal packets at root4 are independently
recomputed, and three balanced fan scales retain exactly the same single
changed edge and exact gap2. Preservation of an existing common frame follows
from the explicitly shared insertion rule. This does not imply successive
translated corridors overlap.

The final nonempty vertex-transitive-core corollary is sound. Vertex
transitivity makes the endpoint and two-cover marginals constant; their sums
are h and q−h respectively, giving every vertex zero imbalance. Thus every
apex cut has c=z=0. Independent controls check all eight cuts of C3 and their
common translation by24 at the first balanced fan. This transports changed
arcs and count differences; it does not manufacture a Hamilton frame shared
by the entire cut cube. Shared-frame preservation retains its stated premise.

## Mixed-bank and product-gate expansion

The N/P bank's packet is independently reproduced. The gcd2 proof uses the
last two weights to restrict the gcd to2 times a power of5, and the first
weight removes that power. This is valid but does not imply interval supply.
The forward/converse trade multiplication gives precisely
(36 Q_S+4 k_S)E14. Reversed order vanishes, and three differences vanish,
so the positive causal quadratic compiler is exact. For a nontrivial spacer
k_S=0 leaves the asserted Q_S u(x)v(y) multiplication term. The explicit
three initial gate tables are independently counted as literal tournaments.
This still retains baseline/linear costs; no physical subtraction is used.

The later additive-energy comparison is independently elementary. The strong
weight separation excludes all nonzero coefficient relations with digits
−2..2, so signed differences have unique digit vectors, support3^m and
multiplicity2^(number of zero digits). Squaring and summing gives6^m.
The note correctly treats the concurrent S93 draft only as inspiration and
does not inherit its unreviewed exponent claim.

## Re-marked escape and valuation schedule

The new root7 packet is independently checked by full and subset permutation
counts. Its first three actual fan counts are independently counted by DP.
The all-length full packet follows by deleting the retained root, counting
paths ending there, and comparing two future fan sizes; the k=0 boundary is
properly excluded. Order2 modulo29 is28, and v_29(2^28−1)=1. The binomial
lift therefore gives exactly one surviving exponent class at every level,
and all exact finite valuations occur infinitely often. This is an ordinary
valuation-control theorem for the one exponential progression, not exact
integer target coverage. No issue was found in that distinction.

## Reproduction and boundary

`python3 -B 05-knowledge/results/session275_exact_selection_20261111/review/verify.py`

Required checks use explicit failures, not removable assertions. Normal and optimized runs agree. RECEIPT.json
records their scope and the source hashes seen on replay. Sparse exact digit
selection, a normalized forbidden cofactor, and a finite corridor census do not
individually or jointly establish cofinite coverage. The packet keeps these
three missing implications visible.
