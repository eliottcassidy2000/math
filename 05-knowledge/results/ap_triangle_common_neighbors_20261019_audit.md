# Independent audit: AP triangle common neighbors

**Verdict: PASS. No mathematical or certificate blocker found.**

Audited 2026-10-19, including the added divisor bound in Section 6.
Reviewed the note, checker, JSON certificate, and the cited local non-torsion
proof. No main-checkout or GitHub changes were made.

## Proof audit

- **Quantifiers:** Correct: for each specified integer AP triangle and each
  finite target count, choose rational neighbors first and then one common
  denominator scale `L`. This does not yield one integer scaling for all
  target counts. Root denominators alone suffice for integer anchors:
  `L²x=(Lp)²−L²a` is integral.
- **Three-anchor criterion:** Correct for every positive rational `a,d`.
  The maps between three AP squares and a positive rational right triangle
  of area `d` are inverse at the level needed. The impossible rational
  isosceles-right case rules out degeneracy. The congruent class depends on
  `d`, while positivity and avoidance of anchors follow from unboundedness.
- **Infinite order and recurrence:** All six valuation cases in the cited
  squarefree-curve argument are correct. First doubling gives negative
  `v₂(X)`, and subsequent doublings lower it by two. Scaling transfers this
  to rational `d`. Compactness and continuity give nonzero rational
  differences of positive multiples tending to the identity. A finite set
  of nonzero points cannot converge to the identity, so a distinct
  subsequence exists. On the real curve this forces `X→+∞`; the triangle
  formulas then force `p²−a→+∞`. No density, rank, or doubling-density claim
  is required. The linked [Elkies course notes](https://people.math.harvard.edu/~elkies/M223.24/index.html)
  do contain the identified September 18 group-law/real-curve discussion.
- **Simple graph and graph type:** The initial witness is the middle anchor
  and must be excluded; the proof does so. The universal conclusion is a
  non-induced `K3 join I_k` subgraph. The induced conclusion is made only
  for the explicitly checked six-vertex certificate. The `K4`, `K3,3`, and
  nonplanarity consequences follow without implying unbounded chromatic
  number or `K5`.
- **Four-anchor obstruction:** Correct for nonzero rational gap, including
  signed anchors and a signed proposed neighbor. Four equations would give
  four distinct rational squares in AP. The cited [van der Poorten paper](https://arxiv.org/html/0712.3850v1)
  explicitly states and proves that impossibility. No finite test is used
  as its proof.
- **Added divisor bound:** Correct. With positive anchors and neighbor,
  `0<p<q`, so `(q−p,q+p)` is a strict positive factor pair of `D` with
  matching parity. The map is injective and gives at most
  `floor(tau(D)/2)` possibilities. For the scaled AP core, the middle
  anchor consumes one equation solution but is forbidden in the graph.
  Therefore external neighbors number at most
  `floor(tau(L²d)/2)−1`, and `k` leaves require `tau(L²d)≥2(k+1)`.
  The strictness caveat for signed/zero-root variants is appropriate.

## Independent exact checks

`ap_triangle_common_neighbors_20261019_audit.py` imports no submitted code.
It checks certificate values directly with fractions and integer square
roots, and certifies elliptic multiples by curve membership plus chord
collinearity rather than the submitted addition routine.

Result: **81 exact checks pass**, also under Python `-O`.

- All six labels are positive and distinct and use the claimed single scale.
- All 15 pairs are checked: the 12 required edges are squares and the three
  leaf pairs are nonsquares. Every published edge root is verified.
- All rational neighbors, root denominators, curve multiples, and the
  non-clique control are verified.
- The submitted normal and optimized executions are byte-identical to its
  stored JSON, with **6,615** submitted checks.

### Complete independent unscaled hostile

Enumerating all **24** unordered positive integer factor pairs `AB=5760`
finds exactly one integer right triangle of area `2880`:

`(A,B,C)=(36,160,164)`.

Any integer common-neighbor solution yields the integer triangle
`(r−p,r+p,2q)` of this area, so this enumeration is exhaustive. Recovering
`p=(B−A)/2`, `q=C/2`, `r=(A+B)/2` gives only

`(x,p,q,r)=(3362,62,82,98)`.

The value `x=3362` is the middle anchor. Thus the original core has **zero
external integer common neighbors**, independently of the submitted
factor-difference and bounded-root scans.

The full 24-pair output is saved in
`ap_triangle_common_neighbors_20261019_audit.json`.
