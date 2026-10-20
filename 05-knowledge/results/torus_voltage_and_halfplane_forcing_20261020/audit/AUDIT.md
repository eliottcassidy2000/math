# Independent audit: half-planes, K7 quotients, and voltage-sensitive returns

2026-10-20. **Mathematical verdict: PASS for G1–G4 and the explicitly scoped
finite controls. Source verdict: PASS for pinned-source provenance and claim
description only; NOT BUILD-VERIFIED and not a complete proof audit.**

The audit found two presentation/provenance defects, not a counterexample to
the four bridge results: an ambiguous machine-readable lattice label, and
source line totals counting a terminal newline as an extra line. Both were
repaired before release. The reviewed final file hashes
and repair status are recorded in `audit_receipt.json`.

## 1. Scope and independent methods

The target is the geometric-forcing packet, not a proof that the Euclidean
plane needs six colours. No new plane chromatic bound is established here.
The audit used independent mathematics, read-only primary-source retrieval,
and small exact calculations. It did not execute third-party source code or
install Lean. It did not change mathematical payloads; authorized metadata-only
corrections are documented below. It made no repository or remote writes.

The producer's standard and optimized Python runs both pass **160 checks** and
reproduce its frozen JSON byte-for-byte. The audit's separate implementation
passes **272 checks**, identically with and without `-O`. It does not import
the producer's code. Its different computational paths are:

- Derive unit steps from the norm equation, rather than assume their list.
- Count walks by multinomial multiplicities, rather than displacement DP.
- Independently retain period coordinates in a group-ring adjacency DP.
- Compute characteristic polynomials by the Leibniz determinant, rather than
  traces/Newton identities.
- Extend Hamiltonian paths recursively, rather than use the producer's full
  permutation enumeration or subset DP.
- Change section representatives by an explicit nonconstant gauge and check
  the closed-triangle voltage sums are unchanged.

These finite tests corroborate, but do not replace, the all-length proofs.

## 2. G1: half-plane/plane colourability equivalence

**PASS.** For a positive finite palette size k, restriction gives the easy
direction. Every finite point configuration can be translated into either an
open or closed nonempty half-plane without changing its unit edges. A proper
colouring of that half-plane therefore supplies a k-colouring of every finite
plane unit-distance subgraph. The de Bruijn–Erdős compactness theorem supplies
the converse global colouring. The same proof works for any subset of the
plane containing a translate of every finite configuration.

The cited original is de Bruijn and Erdős, *A colour problem for infinite
graphs and a problem in the theory of relations* (1951), Theorem 1. The
[publisher/university record](https://research.tue.nl/en/publications/a-colour-problem-for-infinite-graphs-and-a-problem-in-the-theory-/)
and the [original-paper archive](https://combinatorica.hu/~p_erdos/1951-01.pdf)
support the finite-palette compactness attribution. The latter's theorem text
was available in the indexed primary-paper result; direct PDF opening timed
out. This access limitation does not affect the elementary application.

Necessary boundaries are preserved:

- This uses unrestricted colourings and compactness. It is not a measurable
  colouring transfer, effective finite bound, or explicit construction.
- Translation, not scaling, preserves the relevant exact distance.
- Arbitrarily large disks are sufficient. The packet does not assert that
  every fixed-width strip has the plane's chromatic number.
- Seam gluing is exact only for complete labelled boundary relations and
  every crossing-edge constraint. Separate side colourability or palette
  sizes alone do not supply compatible labels. The seam is not bounded-size.

## 3. G2: the triangular-lattice quotient really is K7

**PASS.** With Q(a,b)=a²+ab+b², the six solutions of Q=1 are the stated unit
steps. Their values under q(a,b)=a+3b mod 7 are distinct and exhaust the
nonzero residues. This makes the quotient K7 and makes the projection a
locally bijective graph homomorphism, not a seven-vertex Euclidean subgraph.

Both vectors u=(-3,1) and v=(1,2) lie in the kernel. Their determinant is -7,
so their span is the full index-seven kernel. An additional exact observation
is that q(a,b)=0 implies Q(a,b)=0 mod 7. Thus every nonzero period has squared
length at least 7, and u,v attain it. The usual three-colouring a-b mod 3 is
proper and changes along both period generators.

The categorical direction is important: a quotient colouring pulls back to
a periodic lattice colouring, giving χ(lattice)≤χ(quotient). The reverse
inference is false. A non-five-colourable quotient rules out only the imposed
periodicity, not all lattice or plane colourings. K7 cannot be simultaneously
realized by seven distinct pairwise-unit points in the plane; already K4
fails, because the two possible third vertices on a unit edge are separated
by √3.

The spectral identity and phase restriction are correct. The infinite symbol
has minimum -3; the seven periodic phases give the six nontrivial quotient
eigenvalues -1. The operators have different domains, so a quotient bound
cannot be silently applied to the unrestricted lift.

## 4. G3: the missing coordinate is full displacement

**PASS.** Each quotient edge has one designated unit step. For representatives
t(r)=(r,0), its voltage is t(r)+d-t(s). In the basis (u,v), the coordinates of
a period displacement (x,y) are

    ((y-2x)/7, (x+3y)/7).

### Coordinate bridge to the Euler-torus companion

The geometry packet uses u=(-3,1), v=(1,2). To distinguish notation, write
U=(7,0), V=(-3,1) for the knot packet's basis (called u,v there). Then
u=V and v=U+2V. Consequently geometry period-basis coefficients (a,b)
convert to knot period-basis coefficients by

    (a,b) → (b,a+2b).

Both pairs represent the same ambient lattice displacement. Zero winding is
unchanged, but marked slope pairs must be converted before comparison. The
Hamiltonian-cycle census below records ambient lattice displacements, not
coefficients in either period basis. See the
[Euler-torus companion](../../euler_torus_hamilton_transport_20261020/EULER_TORUS_BRIDGES.md).

For a closed quotient walk the representative terms telescope, so its total
voltage is exactly its actual lattice displacement. Changing representatives
adds a coboundary and has no effect on this sum. A return to the same quotient
vertex merely says displacement lies in the period lattice; return to the
same lifted vertex requires displacement zero.

The example 0→1→2→0 has displacement (3,-1)=-u, hence is not a closed lattice
walk. The complete unordered-triangle census is 14 zero-winding and 21
nonzero-winding triangles. The third trace is 210 in the quotient and 84 in
the zero-winding coefficient. The factor seven is essential: the coefficient
of the trace counts walks based at the seven selected representatives, not
all infinitely many lattice vertices. Per lattice starting vertex the values
are 12 and 30 respectively.

All undirected moments through length eight match the packet:

    lift:     1, 0, 6, 12, 90, 360, 2040, 10080, 54810
    quotient: 1, 0, 6, 30, 186, 1110, 6666, 39990, 239946.

Extracting the identity coefficient of the group-ring trace is the correct
repair. Evaluating every monomial at 1 discards the repair.

## 5. G4: orientation gauge, path counts, and return formulae

**PASS.** The distinction is geometric, not a relabelling convention.

- Carousel steps (1,0),(0,1),(-1,1) map to residues {1,2,3}. The height
  a+2b strictly increases, proving absence of every positive-length directed
  closed walk in the lift. The quotient first has spurious returns at length 3.
- Paley steps (1,0),(0,-1),(-1,1) map to {1,2,4}. A zero displacement forces
  the three multiplicities to be equal. Hence the number of length-3m returns
  per lifted starting vertex is (3m)!/(m!)³; other lengths give zero. The
  quotient first acquires extra returns at length 4.

The independent enumeration yields 175 carousel Hamiltonian paths and 189
Paley paths. There are respectively 25 and 27 paths from each fixed starting
vertex. Directed Hamiltonian cycles modulo rotation of the starting vertex
number 17 and 24. Reversal is not being divided out.

Both stated characteristic polynomials are exact. The Paley polynomial is
(t-3)(t²+t+2)³. If α,β are the roots of t²+t+2, then
u_n=α^n+β^n obeys u_0=2, u_1=-1, u_n=-u_(n-1)-2u_(n-2).
Its trace is 3^n+3u_n; vertex transitivity divides this by seven for a rooted
return count. This proves the recurrence for all n, independently of the
finite control through n=12.

### Additional audited boundary: neither Hamiltonian cycle lifts closed

All 24 Paley quotient Hamiltonian cycles also have nonzero winding: their
length is 7, not a multiple of 3. Thus the distinction from carousel is the
existence of some closed directed walks, not closed Hamiltonian lifts.
The exact rooted-at-zero displacement censuses, with each directed cycle
counted once modulo starting vertex, are:

    Carousel: (-7,7):1; (-1,5):14; (0,7):1; (7,0):1.
    Paley: (-7,7):1; (-2,3):7; (-1,-2):7;
           (0,-7):1; (3,-1):7; (7,0):1.

This does not identify Hamiltonian paths/cycles, Euler tours, or embedded
knots. In particular, it does not conflict with a possible zero-flux Paley
Euler tour of length 21. The separate Euler-torus packet owns that calculation.

## 6. Amenability and the six-colour source claim

**PASS as scoped description, not as an independent proof of the cited
six-colour lower bound.** The plane-isometry group is amenable even when
viewed discretely: translations and planar rotations are abelian, and
reflections add a finite extension. The distinction from the three-dimensional
paradoxical construction is legitimate. Wagon's
[Chapter 10 summary](https://www.cambridge.org/core/books/abs/banachtarski-paradox/measures-in-groups/6F9AE5EC8448E16EEFBD9B9796F4FADF)
explicitly identifies solvable-group amenability and the nonexistence of the
corresponding line/plane paradoxes. None of this supplies a chromatic ceiling;
the integer Cayley graph and free-group tree examples in the packet correctly
separate amenability from high chromatic number.

The audit independently fetched all ten listed files from `openai/math` at
commit `1581d4cbdd69a0faf34b2089614fe2a814316a6a`, recomputed each Git blob
hash from its bytes, and checked every recorded declaration position.
The read-only main-ref endpoint independently returned this same commit.
`source_verification.json` records the resulting evidence without copying the
primary sources into the proposed publication.

The relevant source boundaries are as the producer states:

- The [actual definition](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/lean/OAI/Geometry/PlaneColoring/Coloring.lean#L7-L8)
  is the ordinary exact unit-distance proper-colouring predicate. It does not
  bake in the desired impossibility. The comparator uses the same predicate.
- The comparator's `sorry` is in a challenge statement, not the actual
  solution module. Validation configuration declares permitted axioms and
  disables Nanoda; configuration is not a successful validation receipt.
- The [final solution theorem](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/lean/OAI/Geometry/PlaneColoring/Five.lean#L124-L128)
  visibly invokes proper-to-Borel-weak transfer and the measurable-weak
  obstruction. `lean/OAI.lean` imports that solution module at line 479.
- The [introduction](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/build/sections/introduction.tex)
  claims an unrestricted lower bound and an existential finite obstruction.
  This is not an explicit coordinate realization or independently checked
  five-colour UNSAT certificate.
- The [transfer section](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/build/sections/transfer.tex)
  does not simply assert that an invariant law has measurable samples.
  Continuous-factor identification and Haar rigidity are additional
  load-bearing assertions. The converse's density-point and finite-translate
  mechanism is distinct from this difficult forward direction.

No `lean` or `lake` was on PATH. No solution build, comparator execution,
transitive axiom audit, independent kernel check, or complete trusted-code
review was performed. The audit therefore cannot promote the source claim
to an independently verified theorem. Nor does the audit establish a
refutation of that claim or validate a historical definition-baking allegation.

## 7. Corrections, release gate, and reproduction

The producer replaced `plane_chromatic_number` with
`triangular_lattice_chromatic_number`, and `plane_returns_per_vertex` with
`lattice_returns_per_vertex`. This is a type/scope correction; no numeric
result changed. The auditor then made the coordinated metadata-only repair:
count source lines with `splitlines()` instead of treating a trailing newline
as an extra line. Exactly ten numeric line-total fields changed. The original
manifest was retained locally, the exact diff is included, and the producer's
final manifest hashes were regenerated. Git blob hashes and all declaration
line anchors were already correct.

The publication overlay resolves the staged Euler-torus crosslink to
`../euler_torus_hamilton_transport_20261020/EULER_TORUS_BRIDGES.md`; the
publisher must verify that accompanying target before writing. Preserve the
corrected manifest and limited source status. Do not include downloaded
primary-source copies or interpret the
audit's finite checks as a new Euclidean lower bound.

From the audit directory, with the producer receipt one directory above:

    python3 independent_checks.py ../exact_checks.json
    python3 -O independent_checks.py ../exact_checks.json

In the current staging layout the input path may instead be
`../session16-geometric-forcing/exact_checks.json`; omitting the argument finds
either layout. `independent_checks.json` is the deterministic audit result.
The producer can independently be rerun normally and with `-O` in a copied
directory. Only explicit exception checks are used, so optimization does not
disable validation.

To recheck primary-source byte provenance, obtain read-only copies of the ten
files at the pinned commit, under their repository-relative paths, then run:

    python3 verify_source_manifest.py ../source_manifest.json SNAPSHOT_ROOT

This verifies all ten Git blob hashes, physical line totals, and 22 recorded
line anchors. It never imports or executes the fetched source code.

## 8. V3 addition: classic lattice nonduplication

The separate [classic proof](../CLASSIC_LATTICE_NONDUPLICATION.md) is **PASS**
for square-lattice translations and triangular-lattice translations semidirect
the order-twelve group D6. It makes no novelty claim and imposes no regularity
on the pieces. A finite partition cannot duplicate the lattice: target counts
would force quadratic box volume to be bounded by linear boundary defects.
Exact inequalities and the explicit cutoff R=3D are proved, including D=0.

This proof was separately checked for its quantifiers, counting direction,
boundary formulas, symmetry invariance, and affine inverse. The independent
finite controls pass 13,885 explicit checks, with identical normal/optimized
outputs; those finite checks do not replace the all-R proof. The
[dedicated audit receipt](classic_lattice_nonduplication_audit.json) records the
scope and exact bounds. This addition leaves the earlier 160 producer and
272 geometric audit checks unchanged and does not promote the external
six-colour source claim.
