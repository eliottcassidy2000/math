# Independent audit of the projective-eight bridge

2026-10-25. **Verdict: PASS** for P1--P4 and the stated marked transport
and loss claims in `projective_eight/RESULT.md`. The classical attribution
is retained; this audit independently checks the finite construction and
the coordinate proofs, not literature-wide priority. No substantive
mathematical discrepancy found.

## 1. The two projective designs

The seven prescribed U triples become the seven binary lines under lambda.
Adjoining the point mapped to zero gives exactly the seven linear planes;
their complements give all seven other cosets. Every three distinct binary
points determine one affine plane, so the S(3,4,8) assertion follows directly.
The independent checker generates PSL from translation and minus inversion,
then PGL by adding multiplication by a nonsquare. Their orbits of the seed
quadrangle have respectively 14 and 28 elements. The latter is independently
compared with the homogeneous equianharmonic polynomial test. The 28 blocks
are two disjoint Steiner systems, each triple occurring twice after their
union. Choosing a chirality is essential.

## 2. The actual conic marking and the missing cube basis

The three inherited differences are 7,3,5 and are independent. The dot-product
identification and reference conic are explicitly additional choices. Applying
them gives the displayed component order `(7,3,5,4,6,1,2,0)` and carries the
chosen design to the fourteen affine half-sets of the conic cube. This is a
bijection of finite labels, not a morphism from a finite-field projective line
to a complex surface.

Each complementary pair of affine half-sets gives one nonconstant Walsh
direction. For M the weights one, two, three have eigenvalues 0,-4,-8.
Thus the seven classes of the Steiner structure do not specify which three
are actual hyperplane relations. The stated shear sends 010 to 110 while
fixing 000; it preserves every affine plane and changes intersection 2 to 0.
This is an exact failure of curve-incidence transport, with a concrete
repair sidecar (the marked three directions), rather than a general no-go.

## 3. Affine PSL is not its usual linear action

The supplied formulas for T and S check on all eight points. Since these
generators are affine, their whole generated projective action is affine.
The independent generator closure has order 168 and 168 distinct linear
parts. It is transitive on the eight points. The minus-inversion map has
translation 001 and linear part interchanging the other two bits; its image
under I+L consists of vectors with first bit zero, so the translation cannot
be removed by any origin change. This single element already prevents a
common fixed origin.

The same linear part alone has four fixed binary points, while the affine
element has none. This additionally checks that the usual linear and the
displayed transitive affine actions are not being silently conjugated.
Exactly six elements of this projective group preserve the marked cube
edges, independently recounted. Preservation of all affine planes does not
imply preservation of the conic intersection graph or of surface geometry.

## 4. The two twenty-eight-element carriers

Multiplication by four has two fixed projective points and two three-cycles.
An invariant pair must consist of the two fixed points, giving exactly one.
The four invariant equianharmonic quadrangles consist of a three-cycle plus
one fixed point. The independent enumeration confirms all four satisfy the
equianharmonic test. Hence the two permutation characters differ at this
element, ruling out both an equivariant bijection and a rational permutation-
module isomorphism. The manuscript correctly calls this a formal action on
eight branch labels; it assumes no PSL automorphisms of a particular
eight-line arrangement. The inherited genus formula `31=28+3` is consistent
with this distinction.

## 5. Reproduction

The independent `check_projective_independent.py` imports no producer code.
It uses generator closure instead of the producer's matrix census, recovers
both designs and all triple multiplicities, tests all 70 four-subsets by a
homogeneous polynomial identity, recomputes affine data and the conic marking,
and checks all stated hostile witnesses. Run it normally and under `-O`.
Saved output is `projective_independent.json`; the source hash is embedded.

The original producer is also replayed separately. These exact controls
support the finite statements; the surface facts and conic equations remain
the explicitly cited inherited dependencies rather than being reproved by
a finite group calculation.
