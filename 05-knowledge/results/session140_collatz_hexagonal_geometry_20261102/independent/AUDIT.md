# Session140 independent geometry audit — frozen checkpoint 1

**Outcome:** PASS after one corrected display. Independently checked written
proofs plus separate exact finite controls; this is not formal verification.
Repository baseline: `00f55b6c70a43b10ef9e4763a8f5b9aefd97c3c3`.
The exact reviewed local bytes are pinned in `SOURCE_PINS.json` using paths
relative to the common packet parent, not machine-specific paths.

## Found defect, correction, and regression

The original Collatz proof (`8384004d7e04695b4534dc4eac99988187e6403a8946dedb29c80ded43ee8b65`)
printed `T^(-3)(7)={56,9}`. The independent forward-only complete-fibre
control failed on this display. Correctly,

```
T^(-1)(7)={14},  T^(-2)(7)={28,9},  T^(-3)(7)={56,18}.
```

The prior unit-filtered check was insufficient to validate the displayed
full fibre: both 9 and 18 are excluded by that filter. The producer corrected
the display and added full-fibre regression checks. Corrected proof SHA256:
`e742e475d53cd2bd8532dc3aaef7ace817ede5fed9803b46477260f9f8bfc07c`.
The sharp four-step shortcut clock and all main theorems are unchanged.

## Collatz arithmetic and geometric obstruction

For a positive 3-unit n, the independently derived first safe inverse gate
has k values `2,1,2,3,4,1` at residues `1,2,4,5,7,8` modulo nine. Taking k
doublings, the legal odd inverse, and `4-k` doublings produces
`F1(n)=(16n-2^(4-k))/3`; five doublings produce `F0(n)=32n`.
Both are strictly larger 3-units and return to n in exactly five standard
steps. For the shortcut map the first alternative is `16n`, and both return
in exactly four steps. Same-depth fibres are disjoint because the iterated
map is a function. Strict increase excludes selected macro-level collisions.
Ordinary walks inside a macro-edge can repeat vertices near the known cycle:
the proof properly does not identify them with a disjoint binary subdivision.

Hence every 3-unit root has at least `2^r` ancestors at distance at most 5r.
Any positive integer reaches a 3-unit within its dyadic valuation plus one
steps, without assuming eventual convergence. Both target lattices have
at most `(2R+1)^2` sites in graph radius R, giving
`D_r >= (2^(r/2)-1)/(10r)` for injective edge-endpoint coding. Uniformly
bounded fibre multiplicity changes only the constant, not impossibility.
The Euclidean packing variant needs positive uniform separation; it does
not prohibit arbitrary unit-distance drawings or topological plane drawings.

## Honeycomb subdivision: external hypothesis and finite conversion

[Kühn's author-hosted primary paper](https://web.mat.bham.ac.uk/D.Kuehn/unvacc3.pdf),
*Minor-universal planar graphs without accumulation points*, Theorem 2,
printed page 4, was independently read. Its disconnected-graph condition
requires a VAP-free drawing where any component containing infinitely many
nested cycles has only finite components beside it. Mere VAP-freeness is
not enough in general. The Collatz proof uses the stronger sufficient case
that **no** pseudoforest component contains an infinite nested-cycle family.
Kühn's theorem is CITED, not reproved or formally audited here.

The remaining deductions pass. Functional outdegree one implies at most
one cycle per weak component; local finiteness makes rooted levels finite.
Nested horizontal intervals and integer heights give accumulation-free
component drawings, and separated component strips preserve this property.
The brick-wall vertical perfect matching contracts to the square grid:
its neighboring lower endpoints `(x±1,y±1)` become the four axial neighbors
under `((x+y)/2,(x-y)/2)`.

An infinite connected minor branch set causes no infinitary edge path:
join its at most three attachment vertices by finitely many finite paths,
then trim to a finite tree. A three-attachment tree has a median with
pairwise internally disjoint arms. Repeated attachments yield zero-length
arms, not an obstruction. Disjoint branch sets give a genuine subdivision.
Orient each resulting path according to the original successor. Every used
vertex has one outgoing edge and at most two incoming edges. Routing lengths
are necessarily unbounded. The theorem supplies a subgraph, not an induced
or spanning subgraph, an effective canonical layout, or a convergence proof.

## Unit-distance and chromatic audit

The HN packet's induced realization theorem is sound: seed each finite cycle
by a regular side-one polygon (nonadjacent chords have length greater than
one), then add tree vertices on their parent's unit circle. At any finite
stage, old vertices and the other old unit circles exclude only finitely
many points. Generic translation seeds new components; fair scheduling keeps
each stage finite and reaches all vertices. This proves existence, allowing
crossings; it does not discover the unknown Collatz components or cycles.
Fixing the labeled equilateral positions of 1,4,2 leaves continuum many
extendible positions for 8; an isometry fixing those three anchors is the
identity, so these labeled realizations are noncongruent.

Every finite subgraph of a functional pseudoforest is 2-degenerate; countable
compactness yields a 3-coloring. One must not say the whole infinite graph
is necessarily 2-degenerate (an infinite cubic tree is a counterexample).
The standard Collatz triangle gives the matching lower bound.

A stronger explicit coloring, independently checked here, is
`c(n)=(bit_length(n)-1) mod3`. If `k=floor(log2 n)`, an even step changes k
by -1. For odd n, `2^(k+1) <= 3n+1 < 2^(k+3)`, so k changes by +1 or +2.
Thus every edge changes color and the standard graph has chromatic number
exactly three, with no convergence premise. This is a theorem by inequalities,
not by the finite sample. Subdivision loses this chromatic information:
routing each edge of the 1,4,2 triangle by two steps yields a bipartite hexagon.

For the exact cycle-lift criterion, necessity is telescoping; sufficiency
uses spanning-tree path sums, checking the remaining edges on fundamental
cycles. Distinct quotient labels guarantee injectivity. In `Z[omega]`,
`a+b omega -> a+2b mod7` maps the six units bijectively to nonzero residues.
The kernel is `(2-omega)Z[omega]` (multiplication matrix determinant seven).
The quotient triangle `0->1->2->0` uniquely lifts to `1,1,-omega`, summing
to nonzero `2-omega`. Thus quotient K7 does not contradict a three-colorable
triangular source and does not give a Euclidean HN lower bound.

## Receipts and scope limits

Run from any directory, preserving sibling packet names:

```
python3 independent/verify.py > independent/normal.json
python3 -O independent/verify.py > independent/optimized.json
cmp independent/normal.json independent/optimized.json
```

The separate verifier uses explicit exceptions, never optimization-removable
asserts. It exhausts complete fibres by forward enumeration (predecessors at
depth d are at most `2^d` times the root), independently of the printed inverse
formula. Other controls include 103,024 arithmetic/color samples, power-of-two
boundaries through exponent 1024, eight macro-tree tests to depth 12, all
1,441 labeled trees on two through six vertices with 77,305 repeated-attachment
triples, all 30 rooted K7 quotient triangles, and 21,504 finite canonical
suspensions. The suspension controls are **not** acceptance of any pending
synthesis text; that audit belongs to a separate checkpoint.

The external HN six lower bound remains externally CITED, not locally
Lean-built, dependency-audited, or independently proved. Inherited rigidity,
Selmer, and modularity results are not promoted by this audit. Graph planarity
is not crossing-free Euclidean unit-distance realizability. No assertion of
Collatz convergence, all-site local filling, or unique geometric coding is
made. There is no single forward Collatz orbit visiting all positive integers:
an orbit either omits 1 or eventually repeats its finite 1,4,2 cycle.
