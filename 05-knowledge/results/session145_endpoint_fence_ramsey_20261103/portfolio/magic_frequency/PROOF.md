# Nine decimal labels, a one-diagonal defect, and an exact Latin-square bridge

**Status.** The elementary statements below are **PROVED**. The fixed finite
checks and arrangement census are **FINITE-EXACT**. The interpretation of the
target numbers as elliptic-curve point counts is inherited from the separately
certified Session143 Frobenius packet, not reproved geometrically here. There
are no claims about acoustics, health, healing, or privileged physical meaning.
No global originality claim is made.

## 1. Inheritance and the objects

The repository's actual prior nine-number set is

    F = {174,285,396,417,528,639,741,852,963}.

Source: `05-knowledge/results/preprint_connections_solfeggio_digits_20261017.md`,
section 4, lemma L4, pinned at revision
`7b76044ca13e6684d1746923285b83dbf83828d2`. Its actual mechanism is decimal
encoding of the three fibers of Z/9 -> Z/3. Its claims about separate preprints
are not dependencies of this note. The later audit
`04-computation/experiments/first_arrival_20261018/audit_F2/f2_solfeggio.py`
checks the same set and distinguishes B(x)=4x+1 from multiplication by 4.
In particular B^3(x)=64x+21=x+3 mod 9. We use only this elementary residue
identity, not any heuristic about ROOT-energy shares or prime coincidences.

For a in {1,...,9}, let [z]_9 be the representative in {1,...,9} and put

    f(a) = 100a + 10[a-3]_9 + [a-6]_9.

Then f(1),...,f(9) is the increasing list above. The hundreds digit a is its
rank and is an injective address. Left rotation of the decimal digits sends
f(a) to f([a-3]_9). This is a C3 action whose three orbits are exactly the
fibers of a mod 3. On Z/9 this digit rotation is the inverse of the inherited
B^3 translation x -> x+3, so the orbit identification retains orientation
when that inverse is recorded. These are encodings of residue fibers, not nine special
real frequencies with a theorem about their physical behavior.

The requested target is

    M = 108 148 122
        140 126 112
        130 104 144.

The inherited Frobenius interpretation uses t=-4,...,4 and
N(t)=126-t^3+15t, with the trace placement

    T = -3  2  4
         1  0 -1
        -4 -2  3.

Session143's source is
`05-knowledge/results/session143_frobenius_magic_followup_20261103/producer/PROOF.md`.
It independently counts all 20 nonsingular short Weierstrass equations over
F5 and their extensions to F125. Here we independently check only its stated
nine-state numerical map and its eight line sums.

## 2. A decisive obstruction: raw labels cannot be magic

**Theorem 1.** No permutation of F is a real 3-by-3 magic square. The same is
true after any nonconstant affine relabeling x -> ax+b.

**Proof.** In a 3-by-3 magic square with line sum S, the center is S/3, hence
one ninth of the sum of all nine entries. Here the sum is 4995 and the required
center is 555, which is absent from F. An affine map with a nonzero a preserves
and reflects equality of line sums; alternatively its required center is
555a+b, not among the images. The excluded constant map gives a degenerate
magic square. QED.

This is stronger than failure of a particular layout. It also excludes an
affine equivalence between the raw frequency set and M.

## 3. The positive bridge: orthogonal Latin coordinates

Define rank coordinates for the original decimal labels by

    u(a) = ((a-1) mod 3)-1,
    v(a) = floor((a-1)/3)-1.

They range independently over {-1,0,1}. The increasing map

    Psi(f(a)) = 126 + 4u(a) + 18v(a)

has target list

    104,108,112,122,126,130,140,144,148.

Use the following Lo Shu rank placement:

    R = 2 9 4
        7 5 3
        6 1 8.

The two coordinate arrays are

    U =  0  1 -1       V = -1  1  0
        -1  0  1            1  0 -1
         1 -1  0            0 -1  1.

Each row and column of each array contains -1,0,1; their nine ordered pairs
are all distinct. They are therefore orthogonal Latin squares. More strongly,
each has all eight numerical line sums zero, including both diagonals.
Consequently

    M = 126 + 4U + 18V

has all eight line sums 378. This gives a constructive map with a preserved
predicate, not just a matching cardinality. It preserves strict order and
the three residue fibers as the three fixed-U fibers, but changes numerical
differences, ratios, actual frequency labels, and the original carry profile.
The source address a recovers every label; forgetting V alone loses a
three-element fiber.

### The one-diagonal precursor

The same rank placement of the raw decimal labels is

    F_R = 285 963 417
          741 528 396
          639 174 852.

Six rows/columns and the main diagonal sum to 1665. The anti-diagonal sums
to 1584, a deficit of exactly 81. This has a structural explanation:

    f(a) = 528 + 111u(a) + (567/2)v(a) + (81/2)v(a)^2.

Indeed the three centers of the consecutive rank blocks are 285,528,852.
Their first differences are 243,324; their second difference is 81. In the
U,V placement, every row/column and the main diagonal has V=-1,0,1, whereas
the anti-diagonal has V=0,0,0. Thus the quadratic term contributes 81 to
seven lines and zero to the eighth.

**Theorem 2.** For real c,p,q,d, the array

    c+pU+qV+dV^2

has seven line sums 3c+2d and anti-diagonal sum 3c. It is fully magic if and
only if d=0. In particular the raw decimal labels have exactly one quadratic
carry obstruction in these coordinates. Psi removes that obstruction and
rescales the two independent coordinates; it is not an affine transform of
the raw numeric label.

The 81 is a base-10 carry calculation, not evidence for a physical role of
that number. Likewise Psi is one mathematically natural choice among an
explicit family below, not a uniquely preferred transform of sound.

## 4. Complete classification for this prescribed placement

**Theorem 3.** Let g be any real-valued function on the nine labels F. Its
values in rank placement R form a magic square if and only if there are
unique real c,p,q such that

    g(f(a)) = c+p u(a)+q v(a)    for every a=1,...,9.

**Proof.** All such arrays are magic by the eight zero-sum checks above.
Conversely every real 3-by-3 magic square has the three-parameter form

    c+x      c-x-y      c+y
    c+y-x    c          c+x-y
    c-y      c+x+y      c-x.

To obtain it, the middle cell is c; opposite cells sum to 2c; the top corners
c+x,c+y fix the remaining top cell and the middle row by the row/column sums.
Thus the vector space of magic arrays has dimension three. The arrays
1,U,V are linearly independent because all nine pairs (u,v) occur, so they
are a basis. QED.

**Exact boundaries.**
- g is strictly increasing in the original frequency order if and only if
  p>0 and q>2p. The successive gaps are p,p,q-2p,p,p,q-2p,p,p.
- In the strictly increasing case, every output is positive if and only if
  c>p+q; in general the corresponding condition is c>|p|+|q|.
- g is injective, without monotonicity required, if and only if p and q are
  nonzero and p/q is not in {-2,-1,-1/2,1/2,1,2}. This follows by testing
  p*delta_u+q*delta_v=0 for coordinate differences in {-2,-1,0,1,2}.
- Constants and coordinate projections are legitimate magic arrays with
  repeated values; they are not bijections of nine states.
- Any polynomial or nonlinear description of a map on these nine points is
  subject to exactly this same three-dimensional classification. No global
  classification for all real inputs is claimed.

The original target is c=126,p=4,q=18. The raw frequency array contains the
independent V^2 term and is outside this space. Among the strictly increasing
maps there is a whole open region p>0,q>2p; the specific coefficients 126,4,18
come from the separately inherited Frobenius construction.

### A canonical closest magic array, if approximation is desired

The unique Euclidean least-squares projection of F_R onto the magic space is

    555 + 111U + (567/2)V.

Its residual is (81/2)(V^2-2/3), equal to 27/2 on the six cells with V=+/-1
and -27 on the three cells with V=0. The residual is perpendicular to 1,U,V,
which are mutually orthogonal as functions on the Cartesian product. Its
squared norm is 6561/2. This is a precise optimization statement, not a claim
that these modified labels are a recognized frequency system.

## 5. A second, distinct map through Frobenius traces

There is also a source-address map

    Phi(f(a)) = N(a-5) = 126-(a-5)^3+15(a-5).

Its ordered image is

    130,108,104,112,126,140,148,144,122.

It is bijective, but unlike Psi it is not order-preserving. The three inherited
frequency fibers map as follows:

    a=1,4,7  -> {130,112,148}, sorted {112,130,148}, mean 130;
    a=2,5,8  -> {108,126,144}, sorted {108,126,144}, mean 126;
    a=3,6,9  -> {104,140,122}, sorted {104,122,140}, mean 122.

Thus each fiber becomes an arithmetic progression of common step 18, and the
three fiber means form an arithmetic progression of step -4. This is the same
Cartesian decomposition {126+4u+18v}, with a fiberwise permutation and a sign
change of U relative to the rank-address map. One must not silently identify
Phi with Psi.

The traces themselves in placement T are not magic: their line sums are
3,0,-3,-6,0,6,0,0. Therefore the cubic transfers magic only after its specific
nonlinear evaluation, not because a magic trace grid already exists there.

### Complete polynomial classification on the fixed trace grid

Put

    U_*(t)=t(t^2-9)(3t^4-57t^2+124)/560,
    V_*(t)=-t(t^2-16)(t^4-12t^2+67)/840,
    Z(t)=t(t^2-1)(t^2-4)(t^2-9)(t^2-16).

Evaluation on T gives exactly U,V. Direct expansion gives the polynomial
identity 4U_*+18V_*=-t^3+15t.

**Theorem 4.** A real polynomial P is magic on this fixed trace placement if
and only if

    P(t)=c+pU_*(t)+qV_*(t)+Z(t)Q(t)

for some real c,p,q and real polynomial Q.

**Proof.** Theorem 3 classifies its nine values. Subtract the corresponding
c+pU_*+qV_*; the remainder vanishes at the nine distinct nodes -4,...,4, hence
is divisible by Z. The converse follows by evaluation. QED.

The degree-seven leading coefficient is (9p-2q)/1680. Consequently, among
polynomials of degree at most six, the only ones magic on this fixed grid are
c+k(-t^3+15t). The second independent nonconstant mode first appears at degree
seven. This does not contradict the earlier arbitrary-rearrangement cubic
classification: a nonconstant affine function of t can have a different
magic placement, but cannot be magic on T.

## 6. Loss ledger and hostile boundaries

| Map/claim | Preserved | Lost or not proved | Decisive check |
|---|---|---|---|
| Decimal encoding a -> f(a) | a, mod-3 fibers, digit-rotation C3 action | no magic property | mean 555 absent |
| Increasing chart Psi | order, fiber partition, bijection, constructed magic arrangement | metric gaps, original carry curvature | raw anti-diagonal defect 81 |
| Trace map Phi | nine-state address, fibers, cubic count values | order, curve equation/isomorphism/group/points | distinct map tables |
| Projection to U or V | one three-fiber coordinate, zero line sums | injectivity and other coordinate | each output occurs three times |
| Transported digit rotation | C3 orbit partition | numerical magic, elliptic-curve operation | Psi-rotated grid has seven sums 378, anti-diagonal 324 |
| All fixed-grid magic functions | eight equal line sums | no theorem on arbitrary placements or global real functions | exact three-dimensional basis |
| Center 126 vs order-five carrier dimension 126 | equality of integers only in this lane | no canonical carrier/map/operation correspondence | replacing 126 by any c keeps the same magic incidence |

The rotation hostile is important: even though a C3 action transports through
either bijection, it does not thereby become a magic-preserving operation.
In the increasing chart it rotates the three V symbols and moves the special
zero symbol on the anti-diagonal. A bare orbit partition is not enough.

For the same reason, the occurrence of 126 as q^3+1 at q=5 cannot by itself
identify a 126-dimensional labeled species carrier. The center can be shifted
arbitrarily without changing this Latin-square construction. A genuine bridge
to another 126 must supply an additional map and a preserved operation; no such
map is claimed here. This lane does not solve a Ramsey number or identify a
knot/fence operation.

## 7. Reproduction

Run from this directory:

    python3 verify.py --output result.json
    python3 -O verify.py --output result.optimized.json
    python3 test_controls.py --output hostile_results.json
    python3 -O test_controls.py --output hostile_results.optimized.json
    python3 verify_source_pins.py

The verifier uses exact integers/rationals, not floating-point tolerances or
Python assertions. It checks the formulas, both maps, all eight line sums,
Latin orthogonality, a rank-six magic constraint matrix, the degree filtration,
projection orthogonality, and two independent enumerations of all target magic
arrangements (constructive parameterization versus all 8! center-fixed
permutations). Both give exactly eight arrangements. Source blobs and excerpts
are pinned separately. Mutated-result rejection is exercised under ordinary
and optimized Python; hostile mathematical examples are also tested directly.
