# Exact exclusion of the seven-thin mixed-(1,2), c=4 maximum-span branch

Date: 2026-10-19. This is a scoped computer-assisted exclusion with explicit integer certificates. It makes no claim of historical novelty or of resolving the general 1/3–2/3 conjecture.

## 1. The statement and its scope

A finite poset is seven-thin here when every vertex is incomparable with at most seven other vertices. A balanced pair has ordering probability in [1/3,2/3] under the uniform distribution on linear extensions.

**Theorem.** Suppose a seven-thin finite poset has no balanced pair. Form its strict >2/3 majority order, and select an incomparable pair x,y of maximum index span, with x before y. Its core Q=I(x) union I(y) cannot have exactly one lower slave b, exactly two upper slaves a1,a2, and four common incomparable neighbors Z={z1,z2,z3,z4}. The order-dual mixed-(2,1) branch is excluded as well.

Here b is incomparable with x and satisfies b<y, each ai satisfies x<ai and is incomparable with y, and every z is incomparable with both roots. The theorem concerns this selected maximum-span root and this specific local morphology. It does not assert that every arbitrary incomparable pair has a balanced core, nor does it settle other seven-thin morphologies.

Combined with the [separate frozen c<=3 mixed-(1,2) exclusion](../mixed12_witness_migration_20261019/mixed12_witness_migration_20261019.md), this removes the entire mixed-(1,2) maximum-span branch in a seven-thin poset: y already has c+3 incomparable core neighbors, so c<=4. The same conclusion holds for the dual branch. The present packet verifies the c=4 step independently of that earlier c<=3 computation.

## 2. Why the finite core family is complete

If there is no balanced pair, every pair has exactly one ordering probability greater than 2/3. These directions form a total tournament. A directed triangle is impossible because its three cyclic comparison indicators sum to at most two in every linear extension. Thus the majority tournament is a transitive total order M extending the poset.

The maximum-span condition implies that every common neighbor z lies strictly between x and y in M. If z were before x, the incomparable pair z,y would have larger span; if z were after y, the pair x,z would. In particular,

    p(x<z)>2/3 and p(z<y)>2/3.

The same maximum-span argument forces every slave incomparable with x but comparable with y below y, and every slave incomparable with y but comparable with x above x. More generally it gives D(x) contained in D(y) and U(y) contained in U(x).

Label x=0,y=1,b=2,Z={3,4,5,6},a1=7,a2=8. The fixed relations are b<y and x<a1,a2. Among nonroots, the only permitted directions between blocks are b<Z, b<ai, and Z<ai; the opposite directions would compare a common neighbor with a root or compare x,y. Relations within Z and between the upper slaves are arbitrary subject to transitivity. Relabel Z in a topological order and, if necessary, exchange the upper-slave names.

After deleting the roots, this is exactly a naturally labeled poset on the seven labels (b,Z,a1,a2): every forward pair is an allowed optional relation. The canonical ideal/upset constructor exhausts these possibilities without duplicates. It produces 96,428 cores. All core vertices must satisfy incomparability degree at most seven. The roots have degrees six and seven; each slave is already comparable with a root. The only degree violation is therefore a common neighbor isolated from every other core vertex. Exactly 79,121 cores survive, and these are the exhaustive target of the certificate packet. Repeated isomorphic posets under different natural labels cause no completeness problem.

## 3. The exact analytic portion

For incomparable u,v define

    delta(u,v)=|D(v)\D(u)|+|U(u)\U(v)|.

The included exchange proof shows that delta=0 gives p(u<v)<=1/2 and delta<=1 gives p(u<v)<=2/3. It uses explicit swaps and three-label injections on actual ambient linear extensions, without a cone relaxation.

The exterior is A below both roots and B above both roots. A core vertex with seven incomparable core neighbors is comparable with every exterior element; the only possible directions are A<vertex<B. For a nonroot v, at most epsilon(v)=7-d_Q(v) exterior elements can be incomparable with v, while epsilon is zero for each root. Therefore

    delta_P(u,v)<=delta_Q(u,v)+epsilon(u)+epsilon(v).

Two opposite bounds at most one give an ambient balanced pair. For a root/common comparison, the sharper one-sided formulas are delta_P(x,z)=delta_Q(x,z)+q_B(z) and delta_P(z,y)=delta_Q(z,y)+q_A(z), where q_A+q_B<=7-d_Q(z). Either upper bound at most one contradicts the maximum-span majority orientation.

The compressed analytic manifest records an explicit witness for every analytically excluded core, tagged by the exact logical scope:

- 49,592 cores: an actual ambient balanced pair;
- a further 7,038 cores: a directed maximum-span majority contradiction;
- union: 56,630 cores;
- residual family: 22,491 cores.

The directed test by itself covers 54,235 cores; its overlap with the balanced-pair test is retained in the count receipts but is not double-counted in the union. The core catalogue ordering has SHA-256 `4bd17845e2d931c28633c6b7733bc9846fd92817756e0721da198fb9e3cf7148`.

## 4. Residual coefficient certificates and arbitrary exteriors

Fix a core extension. Let u be a prefix lying before the first root, and v a suffix after the last root. Let G(u,v) count core extensions with that exact prefix and suffix, recording both their number E and the number N(a<b) with each oriented incomparable comparison.

For an actual exterior, the extension-count vector is a nonnegative combination of these generator columns. Indeed, the normalized left and right interleaving counts have nonnegative boundary-tree increments U(u),V(v), and the vector is

    sum_(u,v) U(u)V(v) G(u,v).

No independence or realizability assumption on arbitrary nonnegative coefficients is required. The actual coefficients form a subset of the cone used for the certificate screen.

If a prefix contains a saturated core vertex, its increment is zero: all A elements must occur before that vertex, so adjoining it and every subsequent prefix letter does not change the left interleaving count. The suffix statement is dual. Hence only generators whose prefix and suffix avoid saturated vertices need be retained. The included saturation note proves the forced directions and every later zero increment explicitly.

For each admissible projected majority order M, the certificate lists positive integer weights w(a,b) on comparisons forward in M, and verifies on every retained generator

    sum w(a,b) [3N(a<b)-2E] <= 0.

The same inequality holds for every actual ambient exterior by nonnegative combination. If all listed comparisons had probability greater than 2/3, the left-hand side after normalization would be strictly positive. This contradicts the majority order. Thus every certified M is impossible.

The residual packet covers every one of the 22,491 cores and every admissible projected order: a core linear extension with x before y and no incomparable-pair index span larger than the x,y span. There are 4,253,947 such orders in the residual family. An actual maximum-span majority order projects to one of these: all exterior A precedes x, all B follows y, and deleting exterior vertices cannot create a larger core span relative to the unchanged x,y interval.

All residual certificates use only one or two rows. No three-row, weight-nine, or linear-programming certificate was needed. The construction made zero LP calls. The standalone verification uses Python integers throughout and has no third-party dependencies.

## 5. Reproduction and certificate normalization

The construction reused already committed checks rather than repeating them. Its raw receipts therefore refer to some prior sampled cases. The independent standard-library replay resolves every such omission by materializing an exact positive one/two-row certificate when needed. The normalized residual manifest in this packet is self-contained: its verification does not depend on earlier sample counts, NumPy, SciPy, the original scan, or any unstored search decision.

Run from the extracted packet directory:

```sh
python -O verify_frozen_c4_packet_20261019.py
```

This independently rebuilds the complete core catalogue; checks every analytic witness; recomputes every residual weighted generator margin directly from core linear extensions; checks every admissible majority order against the stored exact coefficients; and verifies the case/certificate stream hashes. Explicit exceptions keep all checks active under optimization. It writes `replayed_c4_controls_20261019.json`. A `--limit` option performs only a labeled prefix control and must not be treated as the full theorem replay.

The compact whole-core view is optional explanatory data. It extracts already-verified certificates whose comparisons are all forced by maximum span, such as p(x<zi)+p(zj<y)<=4/3. Each is a separate exact per-core inequality; the data do not assert one universal untested inequality for all morphologies.

The cone is only a relaxation of actual exterior distributions. A hypothetical failure would have shown a limitation of that screen, not a counterexample to the balanced-pair conjecture. Here the complete listed family passes. No claim is made about other root counts or the general conjecture.
