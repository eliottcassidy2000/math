# Exact children with mixed binary, ternary and quaternary modes

2026-11-13, Codex. **PROVED IN NOTES** for each explicit frame identity;
**FINITE-EXACT** for completeness of the stated order-eight catalogs.
The universal expansion and notation are in [CORE_REPLACEMENTS](CORE_REPLACEMENTS.md).

## A one-arc, one-vertex +4 family

For every m>=1 there is a strong tournament of count

    h_m = 6*3^m + 6*2^m - 5

whose core can be changed to give exactly 3h_m+4, reversing one old core
arc and adding one vertex. Every old vertex and every arc touching the
old transitive module is retained. In particular:

    25 -> 79,  73 -> 223,  205 -> 619,  577 -> 1735.

The rooted source and target are respectively

    A = (2,12,41,49,39,3), root2,       a=(25,24,6,0),
    B = (2,36,25,19,67,93,15), root6,  b=(79,72,18,0).

Their full response vectors satisfy b=3a+4(1,0,0,0). Consequently
H(B[M])=3H(A[M])+4H(M) for every nonempty tournament M. Taking M=TT_m
and K_3=3^m-3*2^m+3 gives the formula above. Root insertion-slot
enumeration certifies every coefficient, including the zero fourth one;
the order bounds forbid five visits.

The embedding is source0->target3, 1->4, 3->1, 4->5, 5->2, with the
whole source module2 mapped to target module6. Append target0. Comparing
every pair shows that the only reversed old arc is source3->0;
`renewal.json` stores the
embedding and literal finite checks. This family is a physical mixed-mode
lift of the previously known scalar instance25->79, not a claim that
79 or the subsequent scalar counts were previously unknown.

## Other certified families

Each row holds for every m>=1. Cost is (arcs per module vertex, fixed core
arcs, new vertices), measured for the displayed core injection. These
are counts of reversals, not a statement about endpoint neutrality.

| source h_m | target | cost | first instance |
|---|---|---|---|
| 14*3^m+6*2^m-13 | 3h_m+2 | (1,1,1) | 41->125 |
| 6*4^m+36*3^m-17*2^m-11 | 3h_m+2 | (0,3,1) | 87->263 |
| 12*4^m+50*3^m-30*2^m-15 | 3h_m+4 | (0,4,0) | 123->373 |
| 6*4^m+52*3^m-5*2^m-15 | 3h_m+2 | (0,4,1) | 155->467 |
| 12*4^m+40*3^m+3*2^m-15 | 3h_m+4 | (0,2,1) | 159->481 |

The third row is a **same-order** transformation: four fixed old core
reversals, with no inserted vertices and no changes to module incidences.
The target count can grow by three while the vertex count stays fixed
for that individual move; the rule is not asserted to be renewable.
Every source and target in this table is strong. Their exact masks,
ports, full run vectors, injections, and literal checks m=1..5 are in
`renewal.json`. The universal certificate is the complete run vector,
not interpolation in m. Cyclic-module controls explicitly check the
scaled-digit boundary in every row.

## Finite catalogs and an actual chain

The universe is all ordinary tournament isomorphism classes of orders
3..8, restricted to strong frames, and every vertex occurrence as a
possible module port. The input bank has the known ordinary counts;
converse tournaments are retained according to the inherited ordinary
bank, not silently identified. Duplicate response vectors are then
merged, retaining a least-order witness. This quotient forgets the
particular frame and its other endpoint statistics.

- Two-run catalog: 2,477 distinct affine responses and 2,074 rays after
  factoring powers of two from their amplitude. There are149 +2 and115
  +4 normalized-exponent transitions. Their lower exponent thresholds
  are part of the data. Different normalizing shifts can mean different
  module lengths, so these264 transitions are not all arbitrary-M laws.
- Full catalog: 22,817 distinct response vectors, with312 exact
  arbitrary-module transitions b=3a+d e_1. There are228 two-visit,
  80 three-visit and4 four-visit transitions. Their digits are respectively
  (133,95), (44,36), and (2,2) for (+2,+4).
- The full-vector matches give an actual three-step chain for every m>=1:

      3*2^m-1 -> 9*2^m-1 -> 27*2^m-1 -> 81*2^m+1,
                     +2             +2             +4.

  At m=1 this is5->17->53->163. Intermediate frames may be replaced by
  the stored equal-response representative; that is an explicit valid
  construction, but no count-neutral arc path between those
  representatives is claimed. The longest chain in this finite
  full-vector catalog has three edges. Larger frames are outside that
  absence statement.

The independent check compares the visit DP with permutation-and-slots
enumeration at all85 rooted occurrences of ordinary frames through5.
The full catalog additionally has91,268 literal transitive-clone DP
checks; the two-run catalog has9,908. The run expansion supplies the
all-m and all-module interpretation of each certified finite vector.

## Composition boundary

Neither source count nor outgoing digit determines the next full vector.
For example, a successful child of the mixed-mode source above is not
automatically a source of any listed rule. [FINITE_PHASE.md](FINITE_PHASE.md)
proves that repeatedly alternating finitely many one-module frames cannot
renew a +4 step under uniform length increments. This identifies a
specific lost resource, rather than attributing failure to an arithmetic
residue class. Exact target-specific replacement remains viable because
it need not use a finite renewable family grammar.
