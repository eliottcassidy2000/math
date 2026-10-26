# Session71: marked Penrose guards and exact dissection boundary carriers

**Status:** elementary proofs below are **PROVED** in their declared protocols;
controls are **FINITE-EXACT**. Penrose's matching/inflation facts and the
Dudeney paper's optimality result are **CITED**, not independently reproved.
No claim of literature priority, new infinite tiling theorem, or equivalence
between Penrose aperiodicity and Dudeney optimality is made.

## Results first

1. The Dudeney paper's *whole-side* edge matching diagram is exactly a pair
   of involutions on named piece sides, with fixed points for degree-one
   nodes and node-type sidecars. Keeping only their product is insufficient.
2. Marked Penrose kites and darts themselves realize arbitrary Boolean
   conjunction guards in a finite placement protocol, using only translations
   and rotations by multiples of 36 degrees, and with exact nonoverlap.
   Thus the Session70 correlation obstruction transfers to a native geometric
   operation, not merely to a picture or a shared integer five.
3. Seam-isometry cycle consistency reconstructs placements up to a common
   rigid motion, but does not ensure nonoverlap. Five unit squares on a
   locally consistent open seam path already give a counterexample.
4. Gaussian norms, golden lengths and equilateral directions remain distinct
   carriers. Equal-area square/triangle realizability is not a consequence of
   the kite/dart count matrix or Pythagorean norm identities.

The reproducible package is `check_carriers.py`, normal/optimized JSON,
optimization-safe injected failures, `SOURCE_PINS.json`, and external-source
metadata. It does not enumerate all Dudeney dissection possibilities.

## 1. Inheritance and the actual source map

Repository input pin: `c2039c32f9726b09ed0e09d5f5c126fc72975897`.
Full paths, Git blob identities, and SHA256 values are in `SOURCE_PINS.json`.

- **Closest proved mechanism:** Session70 `session70-root/SYNTHESIS.md`,
  sections1,3: continuation kernels depend on the named operation protocol;
  native conjunction guards separate the entire alignment law.
- **Canonical hostile:** even/odd parity laws have identical proper-subset
  marginals. The new contribution here is a geometric marked-tile realization
  of those same test functions with disjoint supports.
- **Corrected near miss:** `fivefold_boundary_20261018g.md`, section3, proves
  that the physical fivefold projection is injective but nondiscrete. Its
  retained conjugate coordinate and bounded window restore discreteness;
  a window alone is explicitly not a Penrose tiling classifier.
- **Least-used relevant sidecar:** whole-side endpoint order, local metric,
  and permitted orientation group. They determine which matching combinatorics
  can be developed geometrically and which later seams can be attached.
- **Pythagorean inheritance:** THM-3333,
  `01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md`:
  `(m+in)^2=(m²−n²)+2mni`, norm `m²+n²`, and the actual Farey/light-cone
  isomorphism preserving determinant squared. This is an exact map, unlike
  identifying every appearance of five with golden geometry.

Bounded targeted history search recovered the fivefold/window packet and
Gaussian-square theorem. Older Penrose-related `.out` files are inspiration,
not a recovered proof of a kite/dart boundary-carrier theorem. The unrelated
Riemannian *Penrose inequality* transport notes must not be substituted for
Penrose tiling results merely because their filenames contain Penrose.

Concept board: named guards; relative frame; boundary pairing; geometric
holonomy; nonoverlap; arithmetic field. Each theorem below says exactly which
coordinates it preserves and what it does not establish.

## 2. Literal edge matching diagrams are two involutions

**CITED definition.** Demaine--Kamata--Uehara, *Dudeney's Dissection is Optimal*,
arXiv:2412.03865v3, Definition2 and section3.1, takes **whole intrinsic piece
sides** as links. Nodes have degree1 or2. A side meeting flat subdivisions can
receive a virtual degree-one node. The construction is not a universal
atomization of all overlapping subsegments.

**Theorem P1 (PROVED; named combinatorial isomorphism).** Let `S` be the set of
whole piece-side occurrences. For assembly `X=A,B`, define `α_X` on `S` by
swapping the two links incident to each degree-two matching-diagram node and
fixing the link incident to a degree-one node. Then the edge-labelled,
bipartition-labelled diagram is isomorphic to

    vertex sets: {A} × (S / <α_A>), {B} × (S / <α_B>),
    edge s:       (A,[s]_A) -- (B,[s]_B).

Conversely these two involutions are recovered uniquely from this labelled
incidence multigraph. Preserve boundary/flat/internal node types as labels.

*Proof.* Each matching-diagram node has one or two incident links. Therefore
its incident set is exactly one orbit of the corresponding involution. Every
whole side has one endpoint in each assembly's node set. Sending a node to
its incident-link orbit and a link to its side name gives mutually inverse
incidence-preserving maps. Parallel edges are retained. No geometric existence
is asserted for an arbitrary involution pair. ∎

Alternating continuation is governed by the two named involutions, with
`α_B α_A` giving a two-step motion. The product is **not** the complete carrier:
on `S={0,1}`, both `α_A=α_B=id` and `α_A=α_B=(01)` have product `id`, while their
diagrams have respectively two one-link components and one two-link cycle.
The checker includes this hostile, credited to the independent carrier audit.

**Flat-side guard.** A long whole side opposite two shorter collinear sides
is not silently replaced by two links. Such a common segment refinement is a
different object and requires a map back to the paper's side/virtual-node
structure. The paper's alternating length relations at split sides would be
lost by pretending that all its links are unit seam atoms.

**Endpoint-order sidecar.** Add each piece's cyclic side order and its side
endpoint incidence. Crossing a shared seam between orientation-preserving
placements reverses the boundary direction, hence swaps its two ordered
endpoint slots. A chain of `r` such seam crossings swaps slots exactly when
`r` is odd. This is the direct permutation mechanism behind the no-flip
endpoint-parity rule of the paper's Lemma10, using its path convention.
Reflections change this bookkeeping and must not be admitted without a
recorded orientation sign. Vertex identifications and angle budgets remain
additional data; the unlabelled dual graph of pieces cannot reconstruct them.

**Relative-frame corollary (PROVED).** Let `G` be the declared group of
admissible boundary relabelings, preserving all fixed local piece data. If the
two assembly involutions separately have orbit representatives `a,b`, their
joint `G`-orbits are classified by

    C_G(a) \ G / C_G(b),

with representative `(a, r b r^{-1})`. Indeed write the two presentations as
`g_a a g_a^{-1}` and `g_b b g_b^{-1}`. Common conjugation removes `g_a`; the
remaining frame is `r=g_a^{-1}g_b`. Choices of representatives change it by
left/right stabilizers, and precisely these changes produce the same joint
orbit. If more assembly labels vary, use stabilizers of those **full labelled
states**, not merely centralizers of their involutions. With fully named
operations, the group may be trivial.

The finite `S4` control has two perfect matchings of the same individual type:
`a=(01)(23)` with `b=a` gives two components, whereas `b=(02)(13)` gives one.
The centralizer has size8 and there are two double cosets. Separate orbit
types therefore fail to determine joint continuation. This is the exact
relative-frame bridge to Session70, not an isomorphism of all its objects.

## 3. A native marked Penrose guard with exact geometry

**CITED primary input (visually checked).** Penrose, *Pentaplexity*, printed
pages33--35, Fig9, specifies kites/darts and the same-color vertex matching
rule. Its kite's two flank vertices are black, its 72-degree tip and
144-degree vertex white; the dart's 72-degree tip and 216-degree reflex vertex
are black, its two 36-degree tips white. The formulas below are our coordinate
realization of those shapes and markings.

Put

    φ=(1+√5)/2,  q=2−φ=φ^{-2},
    U=e^{iπ/5},  V=e^{-iπ/5},
    dot(U,V)=(φ−1)/2.

Write coordinates in the **oblique** basis `(U,V)`; these are not ordinary
Cartesian coordinates. Their squared Euclidean norm is

    N(x,y)=x²+y²+(φ−1)xy.

The dart and kite, in boundary order, are

    D=[(0,0),(1,0),(q,q),(0,1)], colors=[black,white,black,white],
    K=[(1,1),(1,0),(q,q),(0,1)], colors=[white,black,white,black].

Their long edge lengths are1, short edge lengths `φ^{-1}`, and area ratio
`area(K)/area(D)=φ`. To see the angles directly, use physical
`U=(cos36°,sin36°)`, `V=(cos36°,−sin36°)`; `q(U+V)=(φ^{-1},0)`.
These are the marked shapes of Fig9, not rhombi misnamed kites.

Define two candidate placements of the **same** marked kite:

    K_good=K−V,
    K_bad=R K+t,
    R=[[1−φ,−1],[1,0]],  t=(φ,−1).

`det(R)=1` and `N(Rv)=N(v)`. In complex coordinates `R` is multiplication by
`−U/V=−e^{2πi/5}`, a rotation of252 degrees (−108 degrees), so no reflection
is used. Both candidate kites lie in `y≤0` and meet the dart only in the full
segment `[0,U]`; the dart lies in `y≥0`. Their interiors are disjoint.
All three polygons lie in the oblique strip `0≤x≤1`.

For `K_good`, vertex `W=U+V` (white) lands at `U` (white) and vertex `V`
(black) lands at `0` (black). For `K_bad`, vertex `W` (white) lands at `0`
(black) and vertex `U` (black) lands at `U` (white). Thus one guard succeeds
and the other fails, despite the same piece counts, the same simple adjacency
`K2`, the same shared length1, and nonoverlap in both configurations.

**Theorem P2 (PROVED; finite marked-placement protocol).** Tile counts and
simple piece-incidence graphs do not preserve legal Penrose attachment.
This statement concerns the named same-color matching guard, not the
existence of an infinite extension of either patch.

*Proof.* The two displayed placements satisfy all listed geometric
constraints. Their endpoint-color predicates differ at both seam endpoints. ∎

## 4. Arbitrary conjunctions are real translated-subset operations

For `k≥1`, put a dart at `3jU`, for `j=0,...,k−1`. For each hidden bit `b_j`,
put a kite at

    3jU−10V + (K_good if b_j=1 else K_bad).

These are legal initial **disconnected finite marked-tile placements**:
there are no contacting pieces. Every piece of gadget `j` has oblique
`x` coordinate in `[3j,3j+1]`, so different gadgets cannot intersect. Staged
kites have `y∈[−11,−10]`; their darts have `y∈[0,1]`.

The named operation `move(j)` translates that kite by **the same vector
`+10V` for either bit** and accepts exactly when all resulting contacts obey
the same-color rule and interiors are disjoint. The straight translation
path itself also remains nonoverlapping. For subset `J`, perform these named
moves and require all to succeed. Disjoint columns make the acceptance test

    g_J(b)=∏_{j∈J} b_j.

**Theorem P3 (PROVED).** This finite marked-Penrose protocol realizes the
entire Boolean zeta family. Individual states need all `k` hidden bits for
exact named-guard preservation. For arbitrary signed laws on the states, an
exact linear carrier evaluating all these guards has dimension `2^k`.
Keeping only subset guards of order at most`t` has rank
`Σ_{r=0}^t binom(k,r)`.

*Proof.* Singletons distinguish all individual bits. Index rows by subsets
`J` and columns by supports `B`; the guard matrix is `1[J⊆B]`. Ordered by
bitmask it is triangular with unit diagonal. Consequently its full rank is
`2^k`, and every selected row family remains independent. ∎

Uniform even- and odd-parity laws agree on every proper-subset guard: fixing
fewer than `k` coordinates leaves a free coordinate that pairs the two
parities. At `J={0,...,k−1}`, exactly one law assigns mass `2^{1−k}` to the
all-ones state. Hence bounded-order summaries lose correlations observed by
actual later geometric operations. This is a **native realization** of the
Session70 algebraic obstruction. It does not import infinite Penrose
extendibility, recurrence, or undecidability into the finite statement.

## 5. Euclidean cocycles develop frames; they do not certify tilings

Fix finitely many labelled polygon templates `P_i` and a connected declared
seam graph. On directed seam `i→j`, let `h_ij` be its specified local isometry
from frame `j` to frame `i`, with `h_ji=h_ij^{-1}`. Lengths, endpoint ordering,
matching labels, and local opposite-interior-side compatibility are checked.
For no-flip problems take all transitions in `SE(2)`; allowing `E(2)` instead
is a different problem.

**Theorem P4 (PROVED; development criterion).** Global placements `g_i` with
`g_j=g_i h_ij` exist iff every closed seam walk has transition product1.
When they exist they are unique up to common left multiplication by an
allowed Euclidean isometry.

*Proof.* Necessity follows by telescoping `h_ij=g_i^{-1}g_j`. For sufficiency,
choose a root, put its placement equal to identity, and define each `g_i` by
the product along a root-to-`i` path. The closed-walk condition makes this
independent of the path and ensures every edge equation. Any two choices
have the same relative placements and differ by their root placements. ∎

**Hostile to geometric sufficiency.** Take five unit squares with origins

    (0,0), (1,0), (1,1), (0,1), (0,0).

Declare the four consecutive seams in the path, using steps `E,N,W,S`.
Each seam is a full unit edge; every interior square uses two distinct ports;
every local seam has opposite-side interiors. The graph is a tree, so its
holonomy condition is vacuous. Yet square5 is exactly square1: piece-area
sum5, union area4. The same path scheme `E,N,W,W` develops without interior
overlap, with origins `(0,0),(1,0),(1,1),(0,1),(−1,1)`.
Here `P5` is the **declared seam graph**, not the full graph of all accidental
contacts. Such accidental contacts are exactly among the global geometric
conditions the coarse carrier omitted.

Thus a valid reconstruction carrier still needs pairwise interior
disjointness and the target-union predicate. For a fixed regular closed polygon
target, containment plus disjoint interiors plus equal total area suffices:
a nonempty omitted part would contain relative interior of positive area,
since the finite union is closed and the target is the closure of its interior.
Keeping a geometric cocycle is not permission to drop this final guard.

## 6. Arithmetic fields: exact links and prohibited shortcuts

| Carrier | Natural field / algebra | Retained operation | What does not follow |
|---|---|---|---|
| Gaussian/Pythagorean | `Q(i)`, norm `a²+b²` | squaring, rational rotations, determinant pairing | golden inflation, equilateral geometry, minimal pieces |
| Golden lengths | `Q(√5)`, `φ²=φ+1` | conjugation and inflation length scale | ordinary planar rotations all lie in this real field |
| Fivefold planar directions | `Q(ζ5)=Q(ζ10)` | rotations by36-degree multiples and conjugation | window/matching/nonoverlap conditions |
| Equilateral directions | `Q(√−3)=Q(ζ6)` | 60-degree rotations | a rational-area or Gaussian square assembly |
| Dudeney no-flip geometry | `SE(2)` with unrestricted real translations/rotations | actual rigid placements | a restriction to one of the preceding fields |

Our oblique coordinate coefficients lie in `Q(√5)`, but their physical complex
coordinates lie in `Q(ζ5)` because `U=−ζ5³`, `V=−ζ5²`. This distinction matters.
Complex conjugation is a reflection, not an allowed move in a no-flip model;
only when reflections are explicitly admitted may it be used geometrically.

Penrose's cited half-dart regrouping gives the inventory matrix
`M=[[2,1],[1,1]]`: an inflated kite uses two kites and two half-darts; an
inflated dart one kite and two half-darts. Its area eigenvector is `(φ,1)`
with eigenvalue `φ²`. This does **not** make `M` a faithful adjacency or
matching-rule carrier, nor a universal intact-tile subdivision of every
finite patch without boundary corrections.

The prior fivefold packet's missing conjugate coordinate, the present
missing boundary frame, and Session70's missing joint moments are all
precisely specified quotient losses. They are not identical quotients:
each requires its own observer and its own retained sidecar.

### The genuinely golden norm obstruction

**Theorem P5 (PROVED; equal-area similar-copy field obstruction).** No two
nondegenerate polygons, one similar to this kite and one similar to this dart,
can simultaneously have all vertices in a finite abelian extension of `Q`
inside `C` **and have equal area**. Reflected similarities do not evade it.
The ordinary unequal-area native kite/dart pair is not obstructed.

*Proof.* Suppose their vertices lie in such a field `F`; take the compositum
`L=F Q(ζ5)`, again finite abelian. The native templates have complex vertices
in `Q(ζ5)`. Write the two similarities as `aK+t` and `bD+s`, using conjugate
templates if needed. Differences of two corresponding vertices show that
`a,b∈L`; the templates have a nonzero edge and conjugated templates are still
in `Q(ζ5)`. Equal areas and `area(K)=φ area(D)` give

    (b/a) conjugate(b/a)=φ.

In any finite abelian field, complex conjugation `c` commutes with every
Galois automorphism `σ`. Therefore

    σ(w c(w))=σ(w)c(σ(w))=|σ(w)|²≥0.

But the automorphism of `Q(√5)` sending `φ` to `1−φ<0` extends to `L`, since
restriction from the Galois group of a finite Galois extension onto that of
a Galois subextension is surjective. This is a contradiction. ∎

This parallels the norm lane's square/triangle obstruction with `√3/4`
replaced by `φ`. The common theorem concerns **total nonnegativity of norms
in abelian fields**, not an isomorphism of `Q(i)`, `Q(√−3)`, and `Q(√5)`.
Area ratios of unlike native tiles need not themselves be complex norms;
only imposing equal-area rescaling creates the displayed norm equation.
The hypothesis is sharp enough to allow escape outside the abelian class:
`K` and `√φ D` are equal-area copies with vertices in
`Q(ζ5,√φ)`, of degree8 over `Q`. Indeed `[Q(ζ5):Q]=4`, and `φ` cannot be a
square in that field by the same norm obstruction (a real square root would
have norm `φ`). The extension by its positive square root has degree2.
No minimum-piece conclusion is implied.

The root/norm lane separately proves the particularly sharp arithmetic
boundary: equal-area square and equilateral triangle cannot both have all
vertices in any finite abelian number field. This uses total nonnegativity
of complex norms under commuting Galois conjugations, not Fermat's theorem
or the golden numeral5. Its sharp degree-eight `D4` realization is independent
of the no-three-polygon-pieces lower bound.

## 7. What the connection does and does not establish

Demaine--Kamata--Uehara prove the four-piece optimum for polygonal pieces
moved without flipping. Their matching diagrams retain enough necessary
metric/angle/topological constraints to exclude all candidate three-piece
cases. We provide a formal boundary-carrier re-expression and independent
geometric guard hostiles; we do not replace their proof by a tile count,
Fibonacci recurrence, Fermat norm representation, or shared numeral.

The exact bridge to Session70 is the operation kernel: preserve every
named legal attachment, including failure. In the Penrose protocol this
forces the missing Boolean alignment data; in the paper's diagram it forces
whole-side pairing and endpoint continuation; in actual assembly it still
requires development, nonoverlap, and target coverage.

**Open scope:** minimal carriers for restricted connected Penrose patches
that are required to extend to full-plane legal tilings; whether a useful
smaller carrier exists when motions/observations are anonymous or resampled;
and a finite computational verifier for the full Dudeney case analysis.
None is claimed solved by these finite controls.

## Sources

- [Demaine--Kamata--Uehara, arXiv:2412.03865v3](https://arxiv.org/html/2412.03865v3),
  Theorem1, Definition2, section3.4 (version dated31 July2025).
- [Roger Penrose, *Pentaplexity*, Mathematical Intelligencer2,32--37](https://web.ma.utexas.edu/users/radin/Pentaplexity.pdf),
  especially Fig9 and printed pages33--36;
  [DOI](https://doi.org/10.1007/BF03024384).
- Repository source pins and exact source snapshots listed in `SOURCE_PINS.json`.
