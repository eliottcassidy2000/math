# Exact exclusion of seven-thin mixed-(2,2) roots with four common neighbors

## Theorem and scope

Let P be a finite poset of maximum incomparability degree at most seven with no
balanced pair. Orient every pair in the direction whose linear-extension
probability exceeds 2/3, and choose an incomparable pair x,y of maximum index
span in that strict majority order, with x preceding y.

The root cannot have two lower slaves, two upper slaves, and exactly four
common incomparable neighbors.

The proof is a complete finite morphology certificate. Together with the
previous c<=2 and c=3 packets, it excludes the entire mixed-(2,2) branch for
seven-thin maximum-span roots: the root degree is c+3, so c cannot exceed four.
The present c=4 proof does not use the published six-thin theorem. Four common
neighbors force the root span to be at least five. No assertion for every
seven-thin poset and no historical novelty claim is made here.

## 1. Canonical core and exhaustive catalogue

The strict majority orientations form a total order. A directed triangle
would require three event probabilities whose sum exceeds two, while any
linear extension satisfies at most two cyclic comparisons. Every cyclic
tournament has a triangle. Comparable pairs retain their poset direction.

The maximum-span root is critical: `D(x) subseteq D(y)` and
`U(y) subseteq U(x)`. Indeed, if t<x but t is incomparable with y, then t
precedes x in majority order and the incomparable pair t,y has greater span.
If t is comparable with y, its direction must be t<y, since y<t<x would
contradict x being incomparable with y. The upset inclusion is the order-dual
argument. These inclusions force every nonroot vertex incomparable with x alone to
be below y, and every nonroot vertex incomparable with y alone to be above x.

The core Q=I(x) union I(y), including the roots, has ten vertices

`x,y,b1,b2,z1,z2,z3,z4,a1,a2`, labeled `0,1,...,9`.

Each bi is below y and incomparable with x; each aj is above x and incomparable
with y; each zi is incomparable with both roots. All four zi lie strictly
between x and y in majority order, or their span to the opposite root would
exceed the maximum root span.

Nonroot comparabilities can run within the three groups or forward from lower
to common/upper and from common to upper. A reverse relation would contradict
a declared root incomparability by transitivity. Label each group by a linear
extension. Thus every possible induced core occurs among naturally labeled
posets on the eight nonroots, after adding bi<y and x<aj.

There are 2,800,472 naturally labeled eight-vertex posets. For a nonroot vertex
let C be its number of comparable nonroots. A common vertex has core degree
9-C and each lower or upper slave has degree 8-C. The roots each have degree
seven. Therefore degree validity is exactly the condition C>=2 at every
common vertex and C>=1 at every slave. It leaves 1,165,395 cores and rejects
1,635,077.

The verifier checks equality of the entire core set, not just the count.
Discovery added predecessor ideals in forward label order; independent
reconstruction adds successor-closed subsets in reverse label order. A
28-bit relation code records each nonroot core uniquely. Every degree-valid
code must occur exactly once in the analytic classification file.

## 2. Exterior reduction and generator cone

Every vertex outside Q is comparable with both roots and hence lies below
both or above both; a mixed direction would make the roots comparable. Call
these exterior parts A and B. Necessarily A<B. In particular no exterior
vertex lies between the roots, so projected root span equals full root span.
Every incomparable core pair has projected span at most this span.

For a fixed core linear extension w, exterior vertices interleave only with
the prefix before its first root and the suffix after its last root. The two
interleavings are independent. Write their counts L(prefix) and R(suffix).
Appending an available core vertex cannot decrease L, since appending that
vertex to a completed left interleaving gives an injection. The corresponding
prefix and suffix increments U,V are nonnegative. We use unnormalized
counts: `U(empty)=L(empty)=e(A)` and `V(empty)=R(empty)=e(B)`, where e counts
linear extensions and e(empty)=1. Nonempty increments are the differences of
successive prefix or suffix interleaving counts.

Telescoping L*R expresses the ambient extension measure as a nonnegative
combination of generator measures G(u,v), where G(u,v) counts core extensions
beginning with ordered prefix u and ending with ordered suffix v, with both
roots in the unprescribed middle between these boundary words. Equivalently,
u and v contain neither root, and are disjoint boundary words of the counted
extension. These
coefficients can have further realizability restrictions. A certificate valid
on every generator remains valid with all such restrictions.

If a core vertex already has degree seven, it is comparable with every
exterior vertex. The direction is forced: all A precede it and it precedes
all B. Consequently every prefix increment at or after the first saturated
vertex vanishes, and every suffix increment at or before the last saturated
vertex vanishes. The included saturation proof establishes the vanishing of
all subsequent increments, not only the increment at the saturated vertex.
Generators whose prefix or suffix contains a saturated vertex therefore have
zero coefficient and may be omitted.

## 3. Analytic exclusions

For incomparable u,v set

`delta_P(u,v)=|D_P(v)\D_P(u)|+|U_P(u)\U_P(v)|`.

The included exchange-injection proof shows that delta<=1 implies
`p(u<v)<=2/3`; delta=0 gives the stronger bound 1/2. Two opposite bounds <=1
prove balance. One bound forbids its direction in the strict majority order.

Use roots and saturated core vertices as seeds. Their upward closure F_A has
all A below it, and their downward closure F_B has all B above it. Set
alpha(u)=0 on F_A and otherwise `7-degree_Q(u)`; set beta(v)=0 on F_B and
otherwise `7-degree_Q(v)`. Then

`delta_P(u,v)<=delta_Q(u,v)+alpha(u)+beta(v)`.

Only missing lower-exterior downset relations at u and missing upper-exterior
upset relations at v can increase the directed defect. Forced directions make
the corresponding contribution zero, and residual incomparability slack
bounds every other contribution. The included side-specific proof supplies
the details.

These bounds certify 733,492 cores with a balanced pair. Another 216,636 have
a cycle of forced majority arcs. Allowed arcs are poset relations, x<y,
x<zi<y, and the reverse of every direction with defect bound at most one.
Every stored balanced witness and every arc of every stored cycle is checked
independently. The remaining core count is exactly 215,267.

## 4. Exact certificates for every residual order

For each residual core, enumerate all its linear extensions that could be
the projected majority order. Require the root to be forward and every
incomparable core pair to have span at most the root span. Reject any order
containing a forward direction forbidden by the defect bound.

Every surviving order is assigned a certificate of forward incomparable
comparisons E_i with positive integer weights w_i. For each retained generator
let E be its number of extensions and N_i the number realizing comparison
E_i. Independent replay verifies the integer inequality

`3*sum_i(w_i*N_i)<=2*E*sum_i(w_i)`.

It holds for every nonnegative mixture of retained generators, including the
actual ambient extension measure. It contradicts the requirement that all
selected forward probabilities exceed 2/3. The non-strict inequality suffices
because the majority probabilities are strictly greater than 2/3.

The search first tries short certificates, then may use numerical linear
programming to suggest weights. Suggested weights are converted to positive
integers and checked exactly on every generator. Independent replay uses only
Python's standard library and unbounded integers. Solver feasibility is not
accepted as mathematical evidence.

Each construction record was committed only after independent exact replay.
The compact final packet uses four independently hash-checked residual streams,
partitioned by residual index modulo four. It is replayed again. Every core extension, retained
generator, event count, integer inequality, and assignment to a forward
certificate is recomputed. A deterministic order-assignment digest provides
an additional check. The replay receipt reports the exact projected-order
and certificate counts, maximum support, and maximum integer weight.

The completed residual totals are 113,423,189 orders excluded by
directed bounds and 2,902,976 orders covered by generator certificates.
There are 184,897 stored per-core certificates: 72,131 one-row,
112,616 two-row, and 150 larger certificates. The latter
were obtained from LP suggestions and then checked exactly. A term pattern
reused on different cores is counted separately. Maximum support
is 12, and the maximum integer weight is
164,426,573,285,821,101,279,024,000. There are zero short-span exclusions.

## 5. Coverage and reproduction

The degree-valid core set is partitioned into 950,128 analytic exclusions and
215,267 residual cores. The verifier checks the exact residual index set and
its correspondence to the independently reconstructed catalogue. No core or
admissible projected order is omitted. Sections 3 and 4 contradict the
assumption of no balanced pair for every listed morphology, proving the
stated theorem.

The earlier results are documented in
[the c<=2 packet](../mixed22_c02_exact_exclusion_20261019/mixed22_c02_exact_exclusion_20261019.md)
and [the c=3 packet](../mixed22_c3_exact_exclusion_20261019/mixed22_c3_exact_exclusion_20261019.md).
Only the former uses the published six-thin theorem; the current c=4 result
has no short-span exclusions.

From this directory, reproduce all independent checks without external
packages (tested with Python 3.12.14):

```sh
python -O verify_mixed22_c4_exact_exclusion_20261019.py
```

Alternatively, run four independent residual parts, with the analytic census
included in part zero, and then check their union:

```sh
for k in 0 1 2 3; do
  python -O verify_mixed22_c4_exact_exclusion_20261019.py --part "$k" > "replay_part${k}.log" 2>&1 &
done
wait
python -O verify_mixed22_c4_exact_exclusion_20261019.py --assemble
```

Every check remains enabled under `-O`. Inputs and proof-code dependencies are
hash-checked, and the assembly verifies all four receipt identities plus the
complete residual set. The manifest records the theorem, self-contained proof
dependencies, compressed certificates, scripts, and receipts. The analytic
census note is a preserved earlier-stage description; the completed theorem
and final replay receipt above govern the full packet.
