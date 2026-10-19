# Exact exclusion of seven-thin mixed-(2,2) roots with three common neighbors

## Theorem and scope

Let P be a finite poset of maximum incomparability degree seven, and suppose
it has no balanced pair. Orient each pair in the direction whose linear-
extension probability exceeds 2/3. Choose an incomparable pair x,y of maximum
index span in that strict majority order, with x preceding y.

The root cannot have two lower slaves, two upper slaves, and exactly three
common incomparable neighbors.

This is a finite morphology theorem, not a proof for all seven-thin posets.
In particular, it makes no claim for mixed-(2,2) with four common neighbors.
The proof below uses exact integer certificates and injection lemmas. Its
short-span screen is vacuous, so this c=3 theorem does not depend on the
published six-thin theorem used by the separate c<=2 packet. Together with
that packet, it covers mixed-(2,2) for c=0,1,2,3. No historical novelty is
asserted.

## 1. Canonical core and majority-order restrictions

The strict majority orientations form a total order. A directed triangle
would require three event probabilities whose sum exceeds 2, whereas every
linear extension satisfies at most two of the three cyclic events. A cyclic
tournament contains a triangle. Comparable pairs have their poset direction.

The root core Q=I(x) union I(y), including the roots themselves, has nine
vertices

`x,y,b1,b2,z1,z2,z3,a1,a2`, labeled `0,1,...,8`.

Each bi is below y and incomparable with x; each aj is above x and incomparable
with y; each zi is incomparable with both roots. The maximum-span property
forces all three zi between x and y in majority order. A common neighbor
outside that interval would have larger span to the opposite root. In
particular the projected root span is at least four.

The possible nonroot comparabilities run within the lower, common, and upper
groups or forward from lower to common/upper and from common to upper.
Reverse comparabilities would make a declared incomparable vertex comparable
with a root by transitivity. Label within each group by a linear extension.
Every possible induced core is therefore represented by a naturally labeled
poset on seven nonroot vertices, adding bi<y and x<aj.

There are 96,428 such naturally labeled nonroot posets. Exactly 13,441 violate
the core degree bound seven, leaving 82,987 cores. Only a common vertex can
have degree eight, and this happens exactly when it is isolated in the
nonroot poset. Inclusion-exclusion gives

`13,441 = 3*4,824 - 3*357 + 40`.

The counts on six, five, and four naturally labeled vertices are independently
reproduced in the preceding c<=2 packet. More importantly, the verifier here
reconstructs all 96,428 seven-vertex nonroot posets directly and checks equality
of the entire degree-valid core set, not just its cardinality. Search generated
predecessor ideals forward; independent replay constructs successor-closed
subsets backward.

## 2. Exterior decomposition and exact generator reduction

Every vertex outside Q is comparable with both roots. It lies below both or
above both, since a mixed direction would make x and y comparable. Denote
these two exterior parts by A and B. Necessarily A<B. Thus no exterior vertex
appears between the roots, and the projected root span equals its full span.
Every incomparable core pair must have projected span at most that root span.

Fix a linear extension w of Q. The A vertices interleave only with the core
prefix before the first root; the B vertices interleave only with the suffix
after the last root. Their choices are independent. Write their interleaving
counts as L(prefix) and R(suffix).

Appending an available core vertex to a prefix cannot decrease L: append that
vertex to any valid completed left interleaving. Hence the increments

`U(empty)=L(empty)`, `U(qv)=L(qv)-L(q)`

are nonnegative. The analogous suffix increments V are also nonnegative.
Telescoping the product L*R writes the ambient extension measure as a
nonnegative combination of generators G(u,v), where G(u,v) counts all core
extensions beginning with ordered prefix u and ending with ordered suffix v.
The generator coefficients may obey additional realizability restrictions;
those restrictions are not needed for a certificate valid on every generator.

A vertex already having seven incomparable neighbors inside Q is comparable
with all exterior vertices. Necessarily all A lie below it and all B above it.
The reverse direction would make it comparable with both roots. Every prefix
increment at or after the first such saturated vertex is zero, and every
suffix increment at or before the last such vertex is zero. Thus a generator
whose prescribed prefix or suffix contains a saturated vertex has zero
coefficient and may be omitted. The included saturation lemma proves the
vanishing for all subsequent increments, not only at the saturated vertex.

## 3. Analytic exclusions with propagated side budgets

For incomparable u,v define

`delta_P(u,v)=|D_P(v)\D_P(u)|+|U_P(u)\U_P(v)|`.

The included directed-defect injection lemma proves that delta<=1 implies
`p(u<v)<=2/3`, with the stronger bound 1/2 when delta=0. Both directed bounds
<=1 imply that the pair is balanced. A single such bound forbids that majority
direction.

Use the roots and core-saturated vertices as seeds. Let F_A be their upward
closure in Q and F_B their downward closure. All A lie below F_A and all B
above F_B. Set alpha(u)=0 on F_A and otherwise `7-degree_Q(u)`; set beta(v)=0
on F_B and otherwise `7-degree_Q(v)`. Then

`delta_P(u,v) <= delta_Q(u,v)+alpha(u)+beta(v)`.

Only lower exterior vertices missing from D_P(u), and upper exterior vertices
missing from U_P(v), can add to these directed set differences. Forced side
relations make the corresponding term zero; otherwise the remaining
incomparability slack bounds it. This proves the displayed inequality and
includes the useful rules alpha=0 for upper slaves and beta=0 for lower slaves.
The full proof is in the included side-specific note.

The analytic phase excludes 47,298 cores with a balanced pair from these
bounds, and 3,800 more by a directed cycle of forced majority arcs. These arcs
include poset comparabilities, x<y, x<zi<y, and the reverse of every direction
whose defect bound is at most one. All 51,098 exclusions are independently
replayed from their stored witnesses.

## 4. Exact residual certificates

There remain 31,889 cores. For each, enumerate every core extension that could
be the projected majority order: the root is forward, every incomparable core
pair has index span at most the root span, and no forward comparison is
forbidden by its directed defect bound.

For each such order the stored certificate consists of forward incomparable
comparisons E_i with positive integer weights w_i. For each retained generator,
let E count its extensions and N_i count the extensions realizing comparison
E_i. The verifier checks exactly

`3*sum_i(w_i*N_i) <= 2*E*sum_i(w_i)`.

The inequality therefore holds for every nonnegative mixture of generators,
including the actual ambient extension measure. It contradicts the claim that
all selected comparison probabilities exceed 2/3. A non-strict certificate
inequality is sufficient because majority probabilities are strictly above
2/3.

Search sometimes used floating-point linear programming to suggest weights.
Each suggestion was converted to positive integers and checked on every
column. Independent standard-library replay recomputes all extension counts
and uses Python's unbounded integers. Numerical solver status is never taken
as a proof.

The completed residual coverage is:

- 9,037,424 projected orders excluded by directed-defect bounds;
- 248,199 projected orders covered by positive generator certificates;
- of the latter, 189,488 use one-row certificates, 53,793 use two-row
  certificates, and 4,918 use rationally exactified LP-derived certificates;
- 2,978 distinct certificates were discovered through LP suggestions;
- zero short-span exclusions.

These categories describe the recorded proof pipeline. The independent replay
regenerates every admissible order and checks its assignment to a forward
positive certificate, including an exact order-assignment digest. It checks
every generator and every integer-weight inequality anew. There are 27,326 distinct certificates: 14,926 one-row, 9,422 two-row, and
2,978 LP-derived. The maximum support is ten, and the maximum integer weight
is 19,986,611,404,658,120,028,280,320. These values are recorded and independently
checked in the included replay receipt.

## 5. Exhaustiveness and reproduction

Every degree-valid canonical core belongs either to one of the 51,098 analytic
cases or to the complete list of 31,889 residual cases. The verifier checks the
core set, exact residual index set, and correspondence of every residual record
to its core. Thus no degree-valid core or admissible projected order is omitted.
The contradiction in Sections 3 or 4 proves the theorem.

The earlier bounded result is documented in
[the mixed-(2,2), c<=2 packet](../mixed22_c02_exact_exclusion_20261019/mixed22_c02_exact_exclusion_20261019.md).
Its published-six-thin dependency is explicit there; none is required for the
present c=3 proof because the three common neighbors force span at least four.

From this packet directory, run with Python 3 and no external packages:

```sh
python -O verify_mixed22_c3_exact_exclusion_20261019.py
```

This writes `independent_replay_receipt_20261019.json`. Every check remains
enabled with `-O`. Both compressed data files are hash-checked first. The
manifest records the self-contained verifier, proof dependencies, data, and
receipts. The original discovery records remain preserved separately; the
repack receipt links them to the compact frozen data without rerunning search.
