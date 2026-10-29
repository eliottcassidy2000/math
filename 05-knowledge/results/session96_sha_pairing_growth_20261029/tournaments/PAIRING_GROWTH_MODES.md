# Tournament growth modes as integral alternating-pairing transports

**Status: PROVED in the stated lattice category, with exact finite controls.**
This is a bridge to the *finite pairing mechanism* behind square orders, not
an arithmetic construction of Sha or a proof of its finiteness. The inherited
source snapshot is `d92d6c0d742c7bab872d0e62769e70e9f0595d2c`.

## 1. Recover the actual operations before comparing their numbers

The repaired [four-state recursive tournament note](../../tournament_recursive_four_state_20261004.md)
already separates these operations:

* endpoint padding K1 -> T -> K1 has n+2 vertices and preserves H(T);
* self order join T -> T has 2n vertices and squares H(T);
* self substitution T[T] has n² vertices and does not generally square H(T);
* transitive cloning and skew-Sylvester doubling are different 2n operations.

For example C3 -> C3 has9 Hamilton paths, whereas C3[C3] has3159 in that
source. These counts are inherited, not needed in the proofs below. The
old [THM-447 skew-Sylvester identity](../../../../01-canon/theorems/THM-447-skew-sylvester-doubling.md)
provides the exact doubling matrix; its path-tiling corrections are retained
in the repaired note. The old
[THM-1475 odd Pfaffian theorem](../../../../01-canon/theorems/THM-1475-the-pfaffian-is-the-odd-function-and-its-spectrum-has-gaps.md)
already proves the odd-square determinant. We extend the comparison to the
whole finite pairing and its named port, not rediscover that determinant law.

We do not use historical claims equating skew-Hadamard normalization with
Hamilton-path normalization. Those are different operations unless an
explicit common transport is supplied. We also do not use the old spectrum
paragraph's phrase "doubly regular tournament on8 vertices": doubly regular
tournaments have odd order; the relevant even matrix is a bordered
skew-Hadamard object.

The arithmetic comparison uses [Poonen–Stoll, *The Cassels–Tate pairing on polarized abelian varieties*](https://math.mit.edu/~poonen/papers/sha.pdf):
a principal polarization gives an antisymmetric pairing, but alternation can
fail. Their characteristic obstruction belongs to the2-primary analysis.
The theorem does not say that every finite pairing is realized as Sha. The
lattice calculations below are proved directly and supply only the alternating
side of this comparison.

## 2. The finite alternating discriminant carrier

Let M be an invertible integral skew-symmetric matrix of even size. Define

    G(M)=Z^n/MZ^n,
    b_M([x],[y])=x^T M^(-1)y modulo Z.

This pairing is well-defined: replacing x by x+Mu adds -u^T y, and replacing
y by y+Mv adds x^T v. It is alternating because M^-1 is skew-symmetric over
Q. It is perfect: if x pairs integrally with every standard basis vector,
then M^-1 x is integral, so x belongs to MZ^n. Its order is |det M|.

A finite abelian group with a perfect alternating Q/Z pairing is an
orthogonal sum of hyperbolic pairs C_d × C_d. Indeed choose an element x of
maximal order d. Perfectness gives y with b(x,y)=1/d. Every element has
order dividing d, so <x,y> is C_d² and its pairing is perfect. It splits off
by subtracting the uniquely determined pairwise coordinates. Iterate on its
orthogonal complement. In particular the group order is a square.

For a tournament dominance matrix, M_ij is±1 off the diagonal and zero on
it. At even n, M modulo2 equals J-I and

    (J-I)^2=I over F2,

because n is even. Thus M is invertible modulo2; det M is odd. Combined
with the Pfaffian identity det M=Pf(M)², its discriminant order is an
**odd square**. This is stronger than a square-or-twice-square restriction:
this carrier has no2-primary group at all.

Consequently no finite linking group made this way can preserve the
nonalternating2-primary obstruction responsible for twice-square Sha orders.
An extra factor C2 with pairing b(1,1)=1/2 gives an abstract antisymmetric
nonalternating model of twice an odd square, but it is added data, not a
feature produced by tournament growth. No assertion is made that such a
finite model is realized by an abelian variety over a global field.

## 3. The n+2 endpoint padding preserves the entire pairing

Reorder the padded tournament so that the old vertices come first, then its
new source and sink. With e the all-ones column, its matrix is

    P(M) = [ M    -e    e ]
           [ e^T   0    1 ]
           [-e^T  -1    0 ].

Write U=[-e,e], J2=[[0,1],[-1,0]]. The integral matrix

    Q = [ I                 0 ]
        [ J2^(-1) U^T       I ]

has determinant1. Direct multiplication gives

    Q^T P(M) Q = (M + U J2^(-1) U^T) directsum J2
                = M directsum J2.

The correction vanishes because the two columns of U are opposite.
Therefore G(P(M)), with its alternating pairing, is isometric to G(M).
The extra hyperbolic unit block has trivial discriminant group. The matrix
congruence is valid even if M is singular; the finite-group interpretation
above is invoked only when M is invertible. Vertex growth n->n+2 is thus
**stabilization**, not a doubling of discriminant order.

For a general two-vertex border with columns u,v and twin arc epsilon=±1,

    Pf([M,u,v; -u^T,0,epsilon; -v^T,-epsilon,0])
      = Pf(M) [epsilon + u^T M^(-1)v].

This rational Schur expression is an integer because the left side is a
Pfaffian. It exposes the joint boundary response. Endpoint padding is the
special zero-contraction case; arbitrary n+2 extension need not preserve
the group or its order. For M=J2,u=(1,1),epsilon=1, using v=(1,1) gives
Pfaffian1 while v=(1,-1) gives3. Equal vector lengths do not determine this
joint contraction.

## 4. Self order join squares the order but not necessarily the group

For an even-order tournament put A=M^-1 and q=Ae. Its self order join has

    K(M)=[ M       ee^T ]
         [-ee^T    M   ].

Since e^T A e=0, block elimination gives

    det K(M)=(det M)^2,
    K(M)^(-1)=[ A       qq^T ]
               [-qq^T   A   ].

If d is the exponent of G(M), then dA is integral and d is its least
positive denominator. The displayed inverse proves

    d divides exponent G(K(M)), and exponent G(K(M)) divides d².

If the named port e is trivial in G(M), equivalently q is integral, the
block elimination is an **integral congruence**, so the full paired group
is the orthogonal direct sum of two copies of G(M). Without this port
condition, cardinality squaring need not mean a direct square of groups.

### Exact four-vertex hostile

Take

    M = [ 0  1  1  1 ]
        [-1  0 -1  1 ]
        [-1  1  0 -1 ]
        [-1 -1  1  0 ].

Then det M=9, q=(-1,1/3,1/3,1/3), and G(M) is C3². Its self order join
has determinant81 and exponent9. Its nullity modulo3 is2, so it has exactly
two cyclic3-primary factors; hence

    G(K(M)) is C9², not C3^4.

The exponent and nullity are directly checked with rational inverse and
finite-field elimination, not inferred from a Smith-form package. In the
same operation H squares, det squares, and Pfaffian magnitude squares, yet
the abelian group acquires an extension. Retaining the marked port class
[e] is the precise additional condition that makes direct-sum elimination
legal in the sufficient case above. We do not assert that its vanishing is
necessary for every possible abstract group isomorphism.

## 5. Skew-Sylvester doubling has an exact integral quotient map

For any integral skew M, including odd size or singular matrices, define

    D(M)=[ M     M+I ]
         [ M-I   -M ].

The inherited identity D(M)²=I2 tensor(2M²-I) proves nonsingularity overQ
because I-2M² is positive definite overR. Put B=I-2M². Add the second block
column to the first, giving

    [ 2M+I  M+I ]
    [ -I     -M ].

The unimodular pivot -I reduces this matrix by integral row/column operations
to diag(I,B). More explicitly the exact quotient map is

    G(D(M)) -> Z^n/BZ^n,
    [(x,y)] -> [x+(2M+I)y],

with inverse [z] -> [(z,0)]. After the displayed block-column addition,
the first block relation eliminates y; the second block relation becomes
exactly B. Equivalently, (-(2M+I)y,y)=D(M)(-y,-y), which proves the
representative replacement directly. Thus the map is a genuine
integral presentation equivalence, not merely equality of determinants.

The transported alternating pairing is

    b([z],[w])=-z^T M B^(-1) w modulo Z.

Indeed D^-1 has top-left block -M B^-1. The formula is well-defined and
alternating. It is perfect: on the quotient2M²=I, so M is invertible with
inverse2M, while B^-1 supplies the ordinary perfect discriminant duality.
The named action of M must be retained; B alone as a symmetric matrix has
a different natural symmetric linking pairing B^-1.

This sharpens the spectral law z->2z+1 for squared skew eigenvalue magnitude:

    |G(D(M))|=det(I-2M²),

but the full pairing is not determined merely by that scalar spectral shift.
For tournament M, B is odd modulo2, so this construction again stays entirely
in the odd-order alternating category. It cannot create the Sha defect2.

## 6. Connection contract and reproduction

Source: a named tournament skew lattice, its boundary port e, and a specified
growth operation. Target: a finite abelian discriminant group with perfect
alternating pairing, plus the exact marked maps above. Preserved: integral
relations and pairing under the stated transports. Lost by determinants:
cyclic exponents, port classes and extension data. Required extra data for a
Sha comparison: an actual global abelian variety, its polarization and
Cassels–Tate pairing, including any2-primary characteristic obstruction.
None is supplied by the tournament construction alone.

Run `python check_pairing_growth.py` and `python -O check_pairing_growth.py`.
The standard-library verifier checks Pfaffian oddness and determinant
squares, exact padding congruence, generic bordering, the skew-doubling
quotient/pairing, and the C3²->C9² self-join hostile. The independent session
audit implements separate checks. The producer has99934 active checks,
including all32834 labelled even-order tournaments of sizes2,4,6 and107
exact transport cases; normal and optimized outputs agree. Universal results rely on the algebraic
proofs, not finite samples.
