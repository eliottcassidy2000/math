# Ten-vertex halving repair: unrestricted quotients cost 3..8; the canonical regular quotient costs 7..15

**Status.** PROVED: the cost formula for a prescribed perfect matching,
the universal lower bound of three reversals, and the three-choice macro
formula below. FINITE-EXACT: the complete minimum
repair census for all 64 labeled cores `Q` in `F_Q(C3)=Q[C3,C3,C3,1]`,
including all 945 perfect matchings of ten vertices; the constrained
census requiring quotient isomorphic to the regular five-vertex tournament
`H5`, over all 24 labeled orientations of `H5`. Both complete enumerations
were independently reproduced by the coordinating session's Mahler/Catalan
lane. No general odd-size extension or Collatz descent implication is
claimed. No novelty claim.

Session: collatz decoder, 2026-09-25, bounded pair-repair follow-up.
The coordinating lane asked whether the failure of exact module halving
can be measured by a small, explicit loss budget.

## Inheritance and target

The inherited construction has three directed triangles and a singleton,
with uniform interblock directions prescribed by the four-vertex core `Q`.
The coordinating lane's exact module obstruction says this construction
does not already possess a halving into pairs. This experiment quantifies
the nearest repair. The canonical hostile is an unmodified cyclic triangle:
no two of its vertices have the same relation to the third. The positive
control is the transitive ten-vertex tournament, whose consecutive pairs
are already modules. The corrected near miss is inferring a structural
halving from even cardinality alone. The least-used sidecar is the list
of arcs changed during a proposed contraction.

Anchor: find the smallest number of arc reversals making every pair of
some perfect matching a module. Niche: reduce the 945 matching choices to
a small formula on the core. Wildcard: dependence on the location of the
singleton inside an otherwise isomorphic core.

The live board is pair modules, core orientation, cyclic defect, singleton
location, and inverse reconstruction. The census changes the board:
singleton location matters even when the unrooted core is unchanged;
cyclic defects have different costs; an optimal contraction still needs a
record of changed arcs to be reversible.

## 1. Exact cost for any prescribed matching

Fix a perfect matching `M` of a tournament on `2m` vertices. All pairs are
modules iff every `2 by 2` block of arcs between distinct pairs is uniform.
Indeed, this condition gives each outside vertex an identical relation to
both members of a pair; conversely two pair modules force the four
crossing directions to agree. Internal pair arcs are unconstrained.

For pair blocks `P,R`, let `k(P,R)` be the number of the four arcs directed
from `P` to `R`. The minimum number of reversals is exactly

    c(M) = sum_{unordered pair blocks {P,R}} min(k(P,R),4-k(P,R)).     (1)

Each block can be made uniform in its majority direction, with either
direction allowed at a tie. The blocks have disjoint edge sets, so the
costs add with no compatibility condition left over. This proves both
necessity and sufficiency of (1).

The script enumerates every matching and also constructs actual repaired
tournaments for every optimum. It then checks the module predicate at
each individual outside vertex, independently of the block-sum formula.

### A universal lower bound of three, independent of the census

For a pair `{x,y}`, let its discrepancy be the number of outside vertices
`z` for which `x->z` and `y->z` disagree. An internal pair in one triangle
has discrepancy exactly one. A pair crossing two triangle blocks has
discrepancy at least two: among the other two vertices in each original
triangle, exactly one distinguishes the pair. A pair involving the
singleton has discrepancy at least one by the same argument in its one
triangle block.

If a perfect matching contains `p` internal triangle pairs, then `p<=3`.
It contains one pair involving the singleton and `4-p` other crossing
pairs. The total discrepancy is therefore at least

    p + 1 + 2(4-p) = 9-p >= 6.

One arc reversal changes at most two discrepancy indicators: the comparison
for each endpoint's matching pair against the other endpoint. Making all
pairs modules requires discrepancy zero, so at least three reversals are
necessary. If there are fewer than three internal pairs, the bound is at
least four. Hence every cost-three optimizer must belong to the natural
family even without enumeration. The stronger all-optima claim at costs
six through eight remains FINITE-EXACT.

## 2. Complete ten-vertex result

Label the core vertices `0,1,2,3`, replace the first three with the
directed triangles on `0,1,2`, `3,4,5`, and `6,7,8`, and replace core
vertex `3` with singleton `9`. Within each triangle, use the cyclic order
displayed. A core mask uses the bit order
`(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)`; bit1 means the lower index points
to the higher index.

| Core and singleton location | Labeled cores | Minimum flips | Optimal matchings per core |
|---|---:|---:|---:|
| Transitive core, any singleton location | 24 | 3 | 27 |
| Source or sink over a triangle, singleton inside that triangle | 12 | 6 | 54 |
| Strong core, deleting singleton leaves a transitive triple | 12 | 6 | 27 |
| Strong core, deleting singleton leaves a cyclic triple | 12 | 7 | 27 |
| Source or sink over a triangle, singleton is the source or sink | 4 | 8 | 81 |

Here "strong core" denotes the four-vertex strongly connected tournament,
whose sorted outdegrees are `(1,1,2,2)`. Source-over-triangle and
sink-under-triangle have scores `(1,1,1,3)` and `(0,2,2,2)`, respectively.
The table is FINITE-EXACT over all 64 cores and all 945 matchings per core:
60,480 candidates, not a random sample.

Permuting the three triangle blocks and reversing every arc leaves the
optimization value unchanged. The 64 cores form six classes under these
operations. The table combines two transitive classes because they have
the same result. The root-deletion bit distinguishes costs `6` and `7`
inside the single unrooted strong-core type.

### Every unrestricted optimum retains one internal pair per triangle

The exhaustive census finds that every optimal matching has exactly one
pair internal to each of the three triangles. Three vertices remain, one
from each triangle, together with the singleton; these four representatives
are paired in one of three ways. Thus all optima lie in a natural family
of just `3^3 * 3 = 81` matchings. The fact that none of the remaining
864 matchings improves or ties the optimum is FINITE-EXACT, not asserted
as a theorem for larger odd blocks.

For each optimal pairing of representatives, the three cyclic rotations
of each original triangle yield 27 optimal matchings. The census verifies
exactly `27`, `54`, or `81` optima as indicated above.

## 3. The cost can be read directly from a three-choice core formula

Choose one internal pair from each triangle, and let `M0` be a perfect
matching of the four representatives, viewed as the core vertices. For
each triangle index `i=0,1,2`, let `{j,k}` be the pair of `M0` not
containing `i`. Define

    d_Q(M0) = # {i in {0,1,2} : Q(i,j) != Q(i,k)}.

Let `t_Q(M0)` be the number of core arcs directed from one pair of `M0`
to the other; its replacement by `4-t` does not affect the expression.
Then the lifted five-pair matching has cost

    3 + 2 d_Q(M0) + min(t_Q(M0),4-t_Q(M0)).                         (2)

**Proof.** Between two internal triangle pairs the core directions are
already uniform, contributing zero. Between the internal pair from
triangle `i` and the representative pair containing its remaining
triangle vertex, the triangle contributes one arc in each direction;
the other representative contributes two identical directions. The block
has one or three forward arcs and costs one. These three blocks give
the initial `3`. Against the other representative pair `{j,k}`, the
internal pair contributes two copies of each of `Q(i,j)` and `Q(i,k)`;
it costs two precisely when those values differ. Finally, between the two
representative pairs the block is exactly the corresponding four core
arcs, costing `min(t,4-t)`. These are all ten cross-pair blocks. QED.

The script checks (2) against (1) for all 64 cores and all three choices.
Together with the complete optimal-family census, the ten-vertex minimum
is therefore computed by taking the minimum of three explicit numbers.

## 4. What the repaired quotient loses

The unrestricted optimization above requires only that a five-vertex quotient exists.
It does not require that quotient to be the coordinating lane's canonical
five-vertex arithmetic tournament, or even to be regular. For example,
the transitive-core case repairs in three flips to a tournament whose
quotient has scores `(0,1,2,3,4)`, rather than `(2,2,2,2,2)`. Representative
optimal quotients for costs `6,7,8` are likewise recorded; tied blocks
can allow more than one repaired quotient.

The minimum costs measure disagreement between the initial tournament
and some pair-module tournament. They are not Collatz edge distances,
nor a decreasing rank along arithmetic steps.

A reversible structural decoder can retain the matching, the repaired
five-vertex quotient, the five internal pair orientations, and the list
of flipped arcs. Those data reconstruct the initial ten-vertex tournament
exactly. Forgetting the flip list destroys information. The unresolved
arithmetic obligation is to show how such retained data certify an allowed
Collatz route, and why its complexity is controlled under iteration.

## 5. Requiring the canonical regular quotient changes the optimizer

The coordinating lane's canonical target is `E10=H5[TT2]`, where `H5` is
the regular tournament on five vertices and each vertex is replaced by
a two-vertex transitive block. There are exactly 24 labeled regular
five-vertex tournaments, all isomorphic. The script verifies this by
enumerating all `2^10` orientations and checking outdegree two at every
vertex, then independently comparing with the orbit of the circulant
orientation `i->j` iff `j-i mod5` is `1` or `2` under all 120 permutations.

For a prescribed matching and a prescribed quotient orientation, each
cross-pair block with `k` forward arcs costs `4-k` if the target directs
it forward and `k` otherwise. Internal pair arcs can remain unchanged:
either orientation is a copy of `TT2`, and swapping the pair's two
vertices supplies the isomorphism. Thus all `945*24` choices exhaust
the cost of reaching `E10` up to isomorphism.

The complete census covers `64*945*24=1,451,520` matching/quotient choices:

| Core and singleton location | Cores | Unrestricted | Canonical `E10` | Best if restricted to three internal pairs |
|---|---:|---:|---:|---:|
| Transitive | 24 | 3 | 15 | 15 |
| Source/sink, singleton inside triangle | 12 | 6 | 13 | 14 |
| Strong, deleted singleton leaves transitive triple | 12 | 6 | 12 | 12 |
| Strong, deleted singleton leaves cyclic triple | 12 | 7 | 7 | 7 |
| Source/sink, singleton at source/sink | 4 | 8 | 10 | 10 |

**The natural-family reduction cannot be reused without its hypothesis.**
In the second row every optimal matching has exactly two internal
triangle pairs; the previously sufficient three-pair family misses the
optimum by one reversal. In the transitive row there are 135 optimal
matchings: 81 with no internal triangle pairs and 54 with three. In the
remaining rows every canonical optimum has three internal pairs.

Counts of optimal matchings and matching/quotient combinations are,
respectively, `135/135`, `27/81`, `27/27`, `27/27`, and `81/162` in the
table's order. The extra quotient choices count distinct target block
orientations, not additional matching choices.

For an explicit counterexample to reusing the natural restriction, take
core mask `5`, with bits indexed as in section 2. A canonical optimum is

    M=((0,1),(2,3),(4,9),(5,6),(7,8)),
    quotient bits=(1,1,0,0,1,0,1,1,1,0),
    flips={(0,2),(0,5),(1,5),(2,5),(2,7),(2,8),(3,7),
           (3,8),(3,9),(4,6),(4,7),(4,8),(6,7)}.

The quotient bits use pair-index order `(0,1),(0,2),(0,3),(0,4),(1,2),
(1,3),(1,4),(2,3),(2,4),(3,4)`. The 13 listed edges are reversed from
their original orientations. Only pairs `(0,1)` and `(7,8)` are internal
to original triangles. The resulting five quotient vertices each have
outdegree two. The full census proves that 12 reversals cannot suffice.

This changes which core appears cheapest: a transitive core costs least
for unrestricted halving and most when the quotient must be canonical.
The strong core whose singleton-deleted triple is cyclic needs seven
reversals in either problem. Therefore the target predicate must be
retained while optimizing a proposed decoder; a low cost to some quotient
does not measure a low cost to the required quotient.

Every canonical optimum is realized in the script, checked for pair
modules and regular quotient degrees, and checked to preserve internal
pair orientations. Root-preserving relabeling and global reversal retain
the counts. The canonical `E10` itself is the zero-cost positive control.

## 6. Reproduction and stopping boundary

    python 04-computation/experiments/decoder_pair_repair_20260925.py

Output: [decoder_pair_repair_20260925.out](decoder_pair_repair_20260925.out).
The script is pure Python, uses exact integer arithmetic, and validates
with explicit exceptions. It prints all 64 results, representative optimal
pairings and flip sets, quotient scores, all-optimum structure checks,
and rooted relabeling/global-reversal checks. It then prints the canonical
quotient results and explicit optimal repairs. The positive controls have
minimum zero, and all 64 inherited hostile constructions have unrestricted
minimum at least three and canonical minimum at least seven. Runtime is
about two seconds in this session.

The bounded probe produces complete local costs for unrestricted and
canonical halving, and exposes a failed transfer of the optimizer between
those targets. It stops here because larger regular odd blocks and
iteration of the repaired map require new arguments. Neither a bounded
global repair budget nor a Collatz certificate follows from this
ten-vertex census.
