# Independent audit of the finite-field surface addendum

20 October 2026, UTC.

## Verdict and exact artifact

**APPROVED.** No substantive mathematical or reproducibility defect remains.
The approval is for the six exact files in `exact_manifest.json`, whose
SHA-256 is

```
bf5919e41985f34c815c7358c59c42fb1337048ee1200813ac84c0cf0c41ba8c
```

The six files are README.md, FINITE_FIELD_SURFACES.md, LITERATURE.json,
verify_surfaces.py, receipt.json, and SHA256SUMS. All source byte counts
and hashes match. They are independently frozen under `frozen/`.
The checksum file covers the five other files; the manifest hashes all six.

The two preliminary README wording suggestions were incorporated before
this final freeze: its checksum-file count is now exact, and its F9/F25
summary says tetrahedron boundaries and icosahedron surfaces. No proof
repair was necessary. The added no-torus corollary is correct.

The reviewer changed no producer file, repository file, remote state, or
member of the separately frozen eight-file prime/field base packet.
This audit does not assert that publication has occurred.

## 1. What was verified independently

The supplied checker passes 4,544 active checks. Its normal and optimized
runs are byte-identical to each other and to `receipt.json`.

The separate `independent_surfaces.py` imports no producer code. It passes
1,673 active checks in both modes, with byte-identical results. It uses
tuple-valued F9 arithmetic, direct graph incidence, oriented integer
boundaries reduced modulo three, face-dual components, and literal matrix
cosets. Unlike the producer, it also tests p=11.

Reproduce from this audit directory:

```
python3 -B frozen/verify_surfaces.py > producer-normal.json
python3 -B -O frozen/verify_surfaces.py > producer-optimized.json
cmp producer-normal.json frozen/receipt.json
cmp producer-normal.json producer-optimized.json
python3 -B independent_surfaces.py > independent-normal.json
python3 -B -O independent_surfaces.py > independent-optimized.json
cmp independent-normal.json independent-optimized.json
```

Both scripts use only the Python standard library. The independent checker
is portable: it does not read the working repository, base packet, or any
machine-specific path. Its q=9 model is F3[t]/(t^2+1), represented as
coefficient tuples rather than the producer's a+3b integer encoding.

If this report is published under the payload's `independent_audit/`
directory, the producer commands instead use `../verify_surfaces.py`
and compare with `../receipt.json`. The independent command is unchanged.
The frozen directory and exact approval metadata above refer to the
retained review-stage snapshot, not additional mathematical dependencies.

All counts and controls are exact. The largest tested graph has 60 vertices
(p=11), 330 edges, and 220 faces. The matrix-coset audit enumerates only
720 matrices in SL2(F9), not a large field or tournament census.

## 2. Correct carrier, characteristic, and local calculation

The vertices are nonzero vectors modulo only the two signs, not modulo all
nonzero scalars. Adjacency is determinant in {1,-1}. The complex contains
all triangles but deliberately no simplex above dimension two. These
choices are necessary and are explicit in the theorem.

At N=[(1,0)] every neighbour has one unique representative (x,1). Its link
has edges x~y iff x-y=+/-1, so its connected components are the additive
cosets of Fp inside Fq. Each is a p-cycle, including C3; there are q/p.
SL2 transitivity transports this description to every vertex.

At any edge choose det(v,w)=1. A common neighbour av+bw requires a,b in
{+1,-1}. Modulo simultaneous sign there are exactly the two distinct
classes [v+w] and [v-w]. Both are distinct from the endpoints. This proves
the exact two-faces-per-edge claim, including p=3.

The two third vertices are adjacent to one another iff +/-2 equals +/-1,
which, in the odd-characteristic scope, is exactly characteristic three.
Thus the filled-clique interpretation would genuinely change the object:
at p=3 it would turn a sphere boundary into a three-simplex. The producer's
explicit F2 hostile is also correct: three vertices, three edges, one face,
with each edge incident to only one face.

## 3. Connectivity and coherent orientation

The connectivity proof does not assume the prime-field surface conclusion.
Every upper unipotent fixes N and hence preserves its graph component.
The displayed S matrix takes N to a neighbour, so it also preserves that
component. Upper and lower elementary matrices generate SL2 over every
field; both printed diagonal-factorization identities were checked by
direct multiplication. Transitivity therefore gives graph connectedness.

There is also a short independent connectivity check: a vector (a,b) with
b nonzero has a common neighbour (x,1) with e1, where x is chosen to make
a-bx=1. A nonzero vector (a,0) is adjacent to (0,a^-1), which is of the
previous type. Thus every class is at distance at most three from N.

The orientation rule a+b+c=0, det(a,b)=1 is well-defined as a cyclic
ordering of classes. For a fixed ordering of the triangle's three classes,
the sign solution to the zero-sum equation is unique up to global sign;
reversal changes the determinant sign and cyclic rotation does not.
The two explicit oriented triples in the proof are correct:

```
(v,w,-v-w),       (-w,v,w-v).
```

They induce opposite directions along the common edge. Hence orientations
cancel across every edge, for primes and extensions alike. The independent
checker verifies each triangle's six signed representatives, its three
cyclic orders, and cancellation along every tested edge.

SL2 preserves this orientation. A determinant-minus-one linear map
preserves adjacency and reverses it, exactly as stated.

## 4. Prime surfaces, counts, genus, and no-torus corollary

When q=p every vertex link is one circle, so the finite complex is a
closed two-manifold; connectivity and orientability were already proved.
When r>1 a punctured small vertex neighbourhood has q/p>1 components,
excluding a surface point. The claimed if-and-only-if is valid.

The counts follow without any finite extrapolation:

```
V=(q^2-1)/2,
E=q(q^2-1)/4,
F=q(q^2-1)/6,
chi=(q^2-1)(6-q)/12.
```

Only for the prime-field surface may chi=2-2g be applied directly,
giving g=1+(p^2-1)(p-6)/24. The values at p=3,5,7 are 0,0,3.
For p>=7 the formula is at least three, so no prime surface or normalized
extension-field component is a torus. Nonplanarity of the one-skeleton
for p>=7 follows from E>3V-6. Neither observation implies a knot model
or a Euclidean polyhedral realization; the note does not claim one.

## 5. Prime-field planes and the normalization proof

Each edge determines W=Fp v+Fp w. The two endpoints are Fq-independent
because their determinant is nonzero, hence also Fp-independent. Signs
and endpoint exchange leave W unchanged. Both incident faces remain in W,
and every edge of either face spans the same W.

A determinant-one Fp-basis of W identifies its induced complex with K_p,
preserving adjacency, triangle incidence, and orientation. Distinct such
planes share no edge: its two endpoints would already span both planes.
They may share several isolated vertices, as the proof explicitly notes.

The assertion that these are the exact face-adjacency components is
justified, not assumed. In a K_p component, faces around each vertex form
a connected cycle in the link. The one-skeleton is connected, so these
local face cycles connect all faces. Therefore no further splitting of
a W-subcomplex occurs.

The planes are exactly g(Fp^2) for g in SL2(Fq). Stabilizing the standard
Fp^2 plane forces every matrix entry to lie in Fp, while determinant one
forces membership in SL2(Fp). The quotient is the left-coset space
g SL2(Fp), not a quotient group; the subgroup need not be normal.
The component count is consequently

```
|SL2(Fq)| / |SL2(Fp)| = q(q^2-1)/(p(p^2-1)).
```

At e1 the incident planes are W_x=Fp e1+Fp(x,1). Two coincide exactly
when x-y belongs to Fp, so they correspond bijectively to the link cycles.
For completeness, the scale used to normalize a second vector belongs
to Fp because the determinant form on W takes values in Fp, and it is
nonzero by independence. This supplies the implicit justification of
the normalization step in section 5.3.

Thus splitting link components separates precisely the W-subcomplexes.
Every original vertex has q/p distinct lifts in q/p distinct surfaces;
open edges and open faces have one lift each. The normalization theorem,
its coset interpretation, and its fibre counts are all valid.

## 6. F9, F25, and independent homology controls

F9 has f-vector (40,180,120), Euler characteristic -20, and links
C3 disjoint-union C3 disjoint-union C3. Its normalization consists of
30 tetrahedron boundaries, with three sheets at every original vertex.

The independent audit explicitly enumerates SL2(F9), its embedded
24-element SL2(F3), and all 30 left cosets. The 24 matrices producing
any one plane are exactly one coset g SL2(F3), checked as literal sets.
These planes equal the face-dual components, and every pair of distinct
planes is checked not to share an edge. The graph has exactly thirty
four-cliques, one per sphere sheet; their interiors are not filled.

Independent simplicial boundary ranks modulo three give:

| q | f-vector | Betti dimensions over F3 |
|---|---|---|
| 3 | (4,6,4) | (1,0,1) |
| 5 | (12,30,20) | (1,0,1) |
| 7 | (24,84,56) | (1,6,1) |
| 11 | (60,330,220) | (1,52,1) |
| 9 | (40,180,120) | (1,51,30) |

These are supplementary exact finite checks, not new universal homology
claims. In particular the F9 result independently rules out treating the
singular complex as a genus-11 surface. Its thirty fundamental classes
are compatible with the thirty oriented normalization sheets.

For F25 the stated (312,3900,2600), 130 normalization components, and
five sheets per original vertex are exact substitutions in the proved
formulas. They were not exhaustively enumerated, and the note correctly
labels that distinction.

## 7. Literature attribution and transport scope

The primary source was checked directly in both HTML and PDF:

- Bogdan Nica, *Unimodular graphs and Eisenstein sums*, arXiv:1505.05034v2,
  section 1.3, printed page 5, and section 8, printed page 15.
  https://arxiv.org/pdf/1505.05034
  https://arxiv.org/html/1505.05034v2

It gives precisely unimodular pairs modulo sign with determinant +/-1
adjacency. Over a field, unimodular means nonzero, so the carrier is exactly
the classical Platonic graph. Nica attributes the prime-ring construction
to Brooks, Perry, and Petersen. Its reference [11] has the DeDeo--Lanphier--
Minei bibliographic data recorded in the packet; these were also verified
against the author's university publication record:
https://scholars.unf.edu/en/publications/the-spectrum-of-platonic-graphs-over-finite-fields/

The graph attribution is supported. No priority claim for the normalization
is made, and no spectral result is needed by the proof. The field/ring
boundary is respected: a modular-curve statement about Z/(p^r) must not
be imported as a statement about F_(p^r). Their neighbour-step additive
orders differ, which is exactly the local-link mechanism here.

At q=5 the extension retains the previously audited marked determinant
graph and fills its triangles. It does not canonically recover field
operations from an unmarked path set, create a path family for every q,
or assert extra tournament automorphisms. No divisor/Collatz trijection,
knot identification, or open-problem resolution is implied.

## Closing boundary

The proof supports the stated all-odd-finite-field theorem. The additional
checks corroborate every delicate finite case and the literal coset model.
All six approved payload members are frozen separately from the untouched
eight-file base packet. This audit makes no claim of completed publication.
