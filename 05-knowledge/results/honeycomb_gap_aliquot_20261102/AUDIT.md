# Audit and reproduction

**VERIFIED same-agent mathematical review and exact replay.** No independent
agent or proof-assistant audit is claimed. Baseline:fb7c88537a4a, including the
latest shared-tail normalization corrections. No theorem ID is reserved.

## Exact checks

`verify.py` uses active exceptions, Python integers and Fractions.

- Hyperbolic coordinates: all127 vertices through depth6 at each of six
  rational bases2,9/4,5/2,3,4,5. All48,006 vertex pairs check sharp separation
  and exact edge identification;40,005 pairs also check the rough-distance
  bounds through rational cosh comparisons. The b2 cross-root hostile and
  general geodesic circle equations are explicit. The all-depth and all-real-b
  conclusions rely on the proof, not interpolation of these samples.
- The incoming two-parameter comparison adds44,919 exact pair checks in23
  configurations through depth5, including three rational equality cases.
  The full minimum-distance graph is checked, not just its number of edges.
  The sharp arcosh(9/7) infimum is restricted to this family and the
  minimum-distance relation; no universal induced-unit-distance bound follows.
- Independent native predecessor enumeration through radius26 and BFS in
  both lattices reproduce704/694 and2263/2107, including their first radius
  violations for root5. The information-theoretic equality case is a direct
  balanced-fibre calculation; it need not be realizable as an adjacency map.
- All132,540 combinations a>b>=1 through a6 and odd primes p,q<=500 check
  the sigma-collision equivalence using independently enumerated divisors.
  The prime guards, repeated-prime9 hostile, supplied numerical roles and
  initial families are retained. Infinitude of the linear prime pairs is open.
- Exact small square graphs at15,20,24,25,26,33 have checked connectivity,
  edge counts and three-colourings. Leaf stripping at15 leaves the explicit
  odd cycle1,3,6,10,15. Full graph enumeration through46 verifies the first
  triangle at30 and first four-cycle at46. No geometric embedding of G25 or
  G33 is inferred from the colouring or absent four-cycles.
- The inherited square25 path enumerator is imported by its explicit path;
  its20 paths, four compatible pairs and26 union arcs are replayed. The
  existing33-state minimality theorem is cited, not rediscovered. All60 even
  permutations on five labels independently check the A5 order histogram.

## Source discipline and corrections

The primary unit-distance source is
[Alexeev, Mixon and Parshall, version2, Theorem1, Table1 and section4](https://arxiv.org/html/2412.11914v2).
Its extrema are edge counts, not chromatic numbers or observer states. Its
rhombus rule is applied only as a necessary initial linear constraint.
The source supplies bounds, not exact values, at25 and26.

[Kuehn's Theorem2](https://web.mat.bham.ac.uk/D.Kuehn/unvacc3.pdf) and the
explicit Session140 construction justify the inherited honeycomb subdivision.
The latter retains its unbounded routes, original markers and unknown
components. The hyperbolic metric formula is the standard one cited by
Session140 from Felikson's university notes; all gap-specific claims are
proved locally. [De Grey's primary paper](https://arxiv.org/abs/1804.02385)
supports the historical five-chromatic example; no assertion about a new
six- or seven-colour lower bound is made here.

The aliquot20 overlap and two-prime-power excess-two family already occur in
`aliquot_divisor_overlap_20261028/THEORY.md`; they are inherited. The new
joint-row collision parametrization is a simple multiplicativity deduction,
with no claim of literature priority. The June15/20/24 reflection is used
through the corrected cardinality atlas, not its retracted uniform-fibre claim.
The old phrase about a3-4-5 triangle at square size12 concerns component
cardinalities; it is not a graph triangle. The literal first C3 is at30.

The prompt's deficit20 is corrected to10. A separate entropy calculation
explains the appearance20/704 under a uniform-source observation, without
identifying it with the abundant integer20. No prior correct704/694 receipt
was altered. The81/27 hostile remains the prior unequal joint-port responses;
this packet does not replace them by scalar geometry or a colouring.

## Reproduction

From the repository root:

```text
python -B -X utf8 05-knowledge/results/honeycomb_gap_aliquot_20261102/verify.py
python -B -O -X utf8 05-knowledge/results/honeycomb_gap_aliquot_20261102/verify.py
python -B -X utf8 05-knowledge/results/honeycomb_gap_aliquot_20261102/plot.py
```

The first two commands need no third-party libraries and regenerate
verification.json. The illustration uses matplotlib and numpy; its floating
coordinates are visual only, not proof data. Dependencies for plotting in
this session were installed into a separate Temp directory, leaving the
global Python environment unchanged. The PNG preview was visually inspected;
the SVG export records the same figure.

Normal and optimized runs produced byte-identical receipts. All local
Markdown links resolve, `git diff --check` passes, and the repository's
agent-facing documentation checker passes its startup, correction and router
checks. These are packaging checks, separate from the mathematical audit.

## Incoming mathematical signal

The clean rebase onto9e1ae5aaf7 brought Session141 and fc2d85f9a8. Session141's
`curvature_triangle/PROOF.md`, section5, independently supplies our gap repair
in the broader (q,W) family, with a sharper upper constant and the near-log2
edge-scale construction. We credit it explicitly. Its warning about close
siblings led to the new exact matching/triangles/tree phase classification
and the restricted arcosh(9/7) infimum. This is additional adjacency control,
not an independent claim to invent Cantor-boundary geometry.

Session141 also proves the complete sigma42 fibre {20,26,41}, the excess-two
classification through tau17, the order26 joint-path polynomial and exact
binary coefficient content. Our row-family parametrization and finite entropy
bound are distinct; the order25 polynomial was already inherited. The
arithmetic/count note links these incoming results without redoing their
censuses. The separate eleven-wedge packet proves sharp order11 stability,
so the current stability gap is12..41. That certificate is inherited here,
not represented as a fresh full census or independent replay in this packet.

The meta-patterns used were **correct the object before sharpening the
technique**, **retain local incidence through counting**, and **keep observer,
recurrence class and finite head**. The boundary-gap theorem is a candidate
mechanism, not a promoted universal meta-pattern from one geometry example.
