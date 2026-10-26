# Gaussian norms, the 44/5 ledger, and the Dudeney boundary

Session 71 norm lane. Base: `c2039c32f9726b09ed0e09d5f5c126fc72975897`.
Statuses: **PROVED elementary statements**, **INHERITED** canon, **CITED** paper,
and **FINITE-EXACT** controls are separated below. No new theorem ID is reserved.

## Findings first

1. Fermat's two-square theorem belongs to the arithmetic of Gaussian norms;
   it does **not** prove the real-geometry, no-flip polygonal four-piece optimum.
2. The inherited PS-D2 boundary-length argument is the direct geometric
   bridge: an equilateral triangle plus its two uncut altitude halves cannot
   rigidly tile a square, in any orientations. Its edge-length field is too
   small. In a general normalized triangle–square dissection, each square
   side must expose at least one piece edge outside `Q(sqrt(3))`.
3. A precise norm obstruction connects the coordinate geometries: in
   `K=Q(sqrt(3))`, the standard equilateral triangle exists but its equal-area
   square has no coordinate realization, in any orientation. The obstruction
   is removed by `Q(3^(1/4))`. It is independent of piece count.
4. The repository already has a stronger-than-scalar Pythagorean carrier:
   raw Gaussian-square images retain represented norms, determinant magnitude,
   and Farey adjacency. A norm alone loses exactly the direction data that
   geometry and matching operations can later inspect.
5. The paper's 44 square cut-pattern types and 5 triangle types count different
   objects from earlier 44s and 5s. Several prior contexts are recovered below;
   no cross-census bijection has been established.

## 1. Paper scope and arithmetic input (CITED)

[Demaine–Kamata–Uehara, v3](https://arxiv.org/html/2412.03865v3), Theorem 1,
excludes at most three polygonal pieces when flips are forbidden. Its §4 uses
triangle side `tau=2`, square side `sigma=3^(1/4)`. Lemma 13 obtains
`sigma=tau` from matching-path length equations; Lemma 14 uses the diameter
inequality `sqrt(2sqrt(3))<2`. The proof retains typed edge/vertex matching
information rather than imposing lattice cuts.

The §5.1 cc-diff buckets for square patterns have sizes
`1,3,12,15,11,2`, summing to 44. Triangle types are `D1,D2,F1,D3,E1`.
Matching cc-diff leaves `12+2*15+2*11=64` pairs; tri-diff removes 18,
leaving 46. These are counts reconstructed from the displayed lists, not
an independent proof that the lists exhaust all dissections.

## 2. What Fermat gives, and what it does not

**CLASSICAL / INHERITED.** For a positive integer `n`,

    n=a^2+b^2 (a,b integers)
    iff every prime p=3 mod4 occurs in n to even exponent.

The prime case is Fermat's theorem; the general form follows from Gaussian
unique factorization and multiplicativity. For a positive rational `q`,
the same statement holds with integer valuations `v_p(q)`, possibly negative,
and `a,b` rational. For reduced `q=m/n`, this is equivalent to `mn` being an
integer two-square norm, since a representation `mn=x^2+y^2` gives
`q=(x/n)^2+(y/n)^2`.

Source checked: Keith Conrad, [The Gaussian Integers](https://kconrad.math.uconn.edu/blurbs/ugradnumthy/Zinotes.pdf),
Theorems 1.2, 8.3 and 9.6–9.9. Repository inheritance is stronger and already
proved: `normseven_arithmetic_f1_20261026/boundary/RESULT.md`, §C3, supplies
the rational norm criterion for the exceptional conic
`a^2+b^2=delta*c^2`, with the delta=2 versus delta=3 hostile. The eight coarse
rational nodes do not by themselves guarantee rational lifts to a resolution.

Here `5=1^2+2^2`, but `44=2^2*11` has odd inert valuation `v_11=1` and is
not even a rational two-square norm. This is a legitimate arithmetic contrast,
not a reason for the paper's type counts.

**PROVED scaling hostile.** Let `k_min(P,Q)` denote minimum number of polygonal
pieces using translations and rotations. For every positive real `lambda`,

    k_min(lambda*P,lambda*Q)=k_min(P,Q).

Indeed, scale a dissection and conjugate each rigid motion by dilation;
the inverse dilation supplies the reverse inequality. Thus equal-area
triangle–square pairs of area 5 and 44 have the same minimum, despite the
different two-square norm statuses. Likewise equal-area square–square pairs
need one piece. Area/norm status cannot classify minimum piece count. This
hostile was independently highlighted by the norm-geometry audit lane.

## 3. A precise field obstruction (PROVED, elementary)

### 3.1 Direct recovered bridge: the boundary-length field

**INHERITED PROVED:** [PS-D2, prime-shells dissection note](sources/10-prime_shells_20260921_dissections.md),
§2, already proves an orientation-free rigid obstruction. An equilateral
triangle of side `a` together with its two altitude-half triangles, left
uncut, cannot tile a square. All piece-edge lengths lie in `a*Q(sqrt(3))`,
so the square side would too, contradicting the required
`(s/a)^2=sqrt(3)/2`. The two altitude halves alone have the analogous
`sqrt(3)/4` obstruction. PS-D1 separately proves that a determinant-one
affine repair preserves areas but changes the equilateral metric; an affine
map is not a legal rigid-piece motion. These results are inherited, not new.

**Immediate general boundary-length lemma.** Let finitely many polygons,
with pairwise disjoint interiors, have union a square. Each square-side
length is a sum of lengths of complete exposed piece edges. Here polygons
are regular polygonal regions (closures of their interiors); zero-area spurs
are discarded before forming the exposed-edge ledger. To prove this, every
piece is contained in the square. Any piece edge whose intersection with a
square side has positive length is collinear with its supporting line,
and hence that complete edge lies on the square side. Such edges cover
the side except possibly finitely many endpoints. Two of them cannot overlap
in positive length: both polygons would locally have interior on the same
inward side, contradicting interior-disjointness. Lengths therefore add.

**Corollary for the normalized Dudeney targets.** For any polygonal dissection
from triangle side2 to its equal-area square, on each of the four square
sides some exposed piece edge has length outside `Q(sqrt(3))`. Otherwise
the boundary-length lemma would put `3^(1/4)` in that field, impossible
by the rational-part argument below. This allows arbitrary translations,
rotations and even reflections. It applies to every finite piece count,
not only three, so does not prove optimality. It explains exactly what
additional cuts must supply when escaping PS-D2: a larger edge-length module.
It does not assert that the cut vertices must leave `K^2`: distances between
`K`-points can themselves lie outside `K`. For example, the rational edge
vector `(1,1)` has length `sqrt(2)`, which is not in `Q(sqrt(3))`.
Nor do four distinct exposed edge occurrences imply four distinct pieces:
one square of side theta, treated as one piece, supplies four such occurrences.
This is a hostile to the counting inference, not a one-piece triangle–square
dissection. The piece-owner labels cannot be discarded.

### 3.2 Related coordinate-field formulation

Let `K=Q(sqrt(3))`, and let

    T=conv{(0,0),(2,0),(1,sqrt(3))}.

Then `area(T)=sqrt(3)`. There is no square of this area with all its vertices
in `K^2`, regardless of position or orientation.

**Proof 1: norm and the other real embedding.** An edge vector `(a,b)` of
such a square would satisfy `a^2+b^2=sqrt(3)`. Apply the field embedding
`j(sqrt(3))=-sqrt(3)`. It produces

    j(a)^2+j(b)^2=-sqrt(3),

impossible in the real numbers. Equivalently, the Gaussian norm
`N_{K(i)/K}(a+ib)` is nonnegative under each real embedding, whereas `sqrt(3)`
is not. More generally the same obstruction works in every totally real
number field containing `K`, because the negative embedding of `K` extends
to a real embedding of that field.

**Proof 2: coefficient audit.** Write `a=r+s*sqrt(3)` and
`b=t+u*sqrt(3)`, with rational coefficients. Equating rational parts gives

    r^2+3s^2+t^2+3u^2=0.

It forces `r=s=t=u=0`, contradicting the required sqrt(3)-coefficient 1.
This proof is specific to `K`; it is a useful independent path.

**Removal and degree boundary.** Put `theta=3^(1/4)`. The square
`(0,0),(theta,0),(theta,theta),(0,theta)` has the required area, and the
triangle also lies in `Q(theta)^2`. Eisenstein at 3 proves `x^4-3`
irreducible. Thus among number fields containing the normalized triangle's
coordinate field, degree four over `Q` is the least possible degree allowing
both coordinate realizations. `Q(theta)` is a real field but not a totally
real field: two conjugates of theta are nonreal. These statements concern
target coordinates, not the coordinate field of a dissection.

**Exact dissection consequence.** There cannot be a finite dissection of this
triangle into polygons with `K`-coordinates followed by rigid motions whose
matrix entries and translation vectors all lie in `K`, assembling the square.
Every moved polygon then has `K`-coordinates. Intersections of nonparallel
boundary lines are solutions of linear equations over `K`, so the polygonal
union has `K`-vertices, contradicting the preceding theorem. This restriction
would exclude **every** finite piece count. Arbitrary real rotations and cut
coordinates are not restricted in the cited optimality theorem. In particular,
the result does not show that rational/K-coordinate cuts alone are impossible
when the motions are unrestricted.

The Penrose lane independently develops a wider abelian/cyclotomic field
obstruction. The elementary `K^2` proof here is its easy-to-audit special case,
not a competing novelty claim.

## 4. Exact inherited Pythagorean connection

**INHERITED PROVED:** [THM-3333](sources/01-THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md)
defines the raw Gaussian-square map

    Phi(m,n)=(m^2-n^2,2mn,m^2+n^2).

For the Lorentz form `<(A,B,C),(A',B',C')>_L=CC'-AA'-BB'`,

    <Phi(u),Phi(v)>_L=2*det(u,v)^2.

On primitive projective integer rays and their compatible raw images,
Farey adjacency is exactly Lorentz pairing 2. This is a genuine graph
isomorphism onto the defined image, not an isomorphism to arbitrary
Pythagorean triangles or to the Dudeney cut-pattern graph.

The minimal inherited norm-collision hostile is

    u=(8,1), u'=(7,4), v=(1,0);
    N(u)=N(u')=65;
    |det(u,v)|=1, |det(u',v)|=4;
    Lorentz pairings=2 and32.

Thus represented norms retain information that the scalar norm discards.
Similarly, raw `(8,6,10)` and `(3,4,5)` form a Farey edge but become the same
unordered primitive triangle after normalization and forgetting leg order.
The hypotenuse5 alone is particularly insufficient.

**INHERITED PROVED:** [THM-3336](sources/02-THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md)
retains primitive Gaussian content under multiplication; raw multiplication
is a similarity, primitive reduction is not an innocent coordinate deletion.
[THM-3358](sources/03-THM-3358-admissible-composite-parabolic-compiler-and-hensel-normal-offset-atlas.md)
uses the stricter primitive norm condition: an odd integer admits a primitive
two-square representation precisely when every prime factor is 1 mod4
(with the trivial norm1 boundary treated separately).

## 5. Historical 44/5 recovery, not numerology

| Located context | Precisely counted object / status | Comparison with dissection |
|---|---|---|
| Paper §5.1 | 44 typed square cut-pattern classes; 5 triangle classes, modulo its typed subdivision equivalence; CITED lists | The directly relevant 44/5 |
| Fresh re-marking note §2, HYP-9416 | Strong-atom rooted-complement distinct target counts5 at n=4 and44 at n=5; source says FINITE-EXACT through7, not independently audited, all-orders OPEN | A genuine recent pair of 5/44 counts, but different objects and equivalence |
| THM-644 small census | 44 pure-black/non-transpose-self classes at order6; mixed classes5 at order5; source is a bounded verified census with conjectural general parity material | Neither 44 cut patterns nor a five-object source category |
| THM-1600 §A | `L((v-1)^m)=!m`; `!5=44`, PROVED exact identity | A particularly clean genuine 44–5 formula, but no dissection map |
| Gaussian norms | 5 represented;44 not represented, PROVED classical criterion | Explicitly not a pattern-count explanation |

The exact intended earlier 44/5 reference is not uniquely determined by the
available prompt. These are recovered candidates with semantics and source
pins, not an assertion that one is the user's intended referent. A set-theoretic
bijection between two 44-element sets always exists but preserves no named
operation automatically. A structural bijection would require corresponding
objects, morphisms/operations, and preserved predicates; none is supplied here.

## 6. Transfer contracts and failure anatomy

| Source -> target / map | Preserved | Lost / required sidecar | Decisive hostile |
|---|---|---|---|
| Gaussian integer -> integer norm | squared length and multiplication | factor allocation, orientation, representation | two norm65 rays with determinant1 versus4 |
| Primitive ray -> raw Phi image | exact Farey adjacency, Lorentz form, norm | already-identified spin sign only | arbitrary leg swap fails the fixed pairing |
| Equal-area targets -> scalar area | area necessary for dissection | shape, piece incidence, direction, connectivity | square/square needs1, triangle/square4 |
| K-coordinate square edge -> Gaussian K-norm | squared edge length | not a loss here; target norm is forbidden | conjugate sqrt(3) is negative |
| Real dissection -> K-restricted model | only instances actually defined over K | unrestricted cuts/motions excluded | extension Q(theta) already removes coordinate obstruction |
| Pattern family -> cardinality44 or5 | size | typed subdivision, incidences and matching constraints | no recovered operation-preserving bijection |

Closest proved mechanism: PS-D2's orientation-free boundary-length field,
represented Gaussian norms and rational conic
solubility. Canonical hostile: norm65 directional collision. Corrected near
miss: general number-field geometry is not lattice geometry, and a scalar
area is not a piece-count invariant. Least-used sidecar: the second real
embedding (total positivity), rather than just the positive physical length.

Concept board: integer norm / represented Pythagorean ray / rational conic /
real-embedding sign / typed dissection matching / the44-and5 census. The
field obstruction connects the first four to geometry, but it simultaneously
proves why that connection cannot distinguish three from four pieces.

## 7. Reproduction and honest boundary

`python3 -B check.py` and `python3 -O -B check.py` produce identical receipts.
The checker verifies source hashes, independent bounded norm enumeration,
Gaussian/Farey identities, the representation-loss hostiles, exact field
coefficient identities, positive extension controls, derangements by explicit
permutations, and the paper-list arithmetic. Bounded checks corroborate proofs;
they do not prove the universal norm criterion, source census exhaustiveness,
or the paper's theorem. The minimum-four result is cited, not re-proved here.

No shared repository file was changed. The startup packet is read-only; message
peek was unavailable because this checkout has no `.machine-id`. No identity
was fabricated. Publication and full cross-lane integration belong to the
parent session.
