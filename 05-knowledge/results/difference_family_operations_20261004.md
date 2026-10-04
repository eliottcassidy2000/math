# Integers as difference families, with retained arithmetic and scaling memory

2026-10-04 (America/Denver).

**Status:** PROVED elementary difference-family, triangular, and stabilization
statements below. FINITE-EXACT for the declared tests. The Syracuse fibre and
two-drop results are inherited canon, independently exercised here. The user's
unparenthesized f expression is not assigned an invented association or law.
The scalar g equation is analyzed literally, then given explicitly different
quotient and family realizations. No Collatz convergence claim follows.

## 1. Inheritance and definitions recovered before deriving

The earlier [operation note](arithmetic_seams_20260921_operations.md), section 1,
explicitly says f and g had not yet received unambiguous definitions. Its
precise proposed model is

    a boxplus b=a+b+1,       a boxtimes b=ab+a+b,
    h(a)=a+1.

These are ordinary addition and multiplication transported through h on
{-1,0,1,...}; their units are −1 and 0. This is a useful existing construction,
not authorization to rename either operation f or g. The
[scaffolding audit](collatz_mod6_20260917_scaffolding_audit.md), sections 4 and 10,
likewise distinguishes its conditional factorial repair from the undefined
four-operation string. The owner's present explicit triangular and scalar
identities are handled in sections 4 and 5 below.

The [paired-prime note](paired_prime_collatz_20261003.md), section 6, already
proves the all-one 12−8 difference as the interchange defect
2(ad+bc) for star(a,b)=a+b+2ab. The
[golden-carrier note](collatz_golden_carriers_20261004.md), section 2, separates
that identity from the accelerated arrows 7->11 and 9->7. The closest
proved mechanism for whole dynamical difference fibres is
[THM-4530, Syracuse drop multiplicities and pair injectivity](../../01-canon/theorems/THM-4530-syracuse-drop-multiplicities-two-copies-over-z-densities-and-pair-injectivity.md).

The canonical hostile is two translated pairs with the same gap but different
componentwise products. The corrected near miss is treating equality of a
difference, a scalar scaling equation, and setwise scaling invariance as the
same predicate. The useful sidecar is the diagonal translation, equivalently
the sum of the pair. The live board is difference classes, raw arithmetic,
triangular decompositions, scaling quotients, and exact dynamical guards.

## 2. A complete arithmetic on whole difference families

For every d in Z let

    D_d={(a,b) in N_0^2 : b−a=d},
    p(d)=(max(−d,0),max(d,0)).

Every member of D_d has the unique form p(d)+(t,t), t>=0. Thus the pairs
(7,11) and (8,12) are representatives of the same integer 4 with translation
memories 7 and 8. Two pairs have the same gap iff suitable common upward
diagonal translations make them equal. No selected representative is needed
to define the class.

Componentwise addition descends to these classes. Componentwise multiplication
does not: multiplying (7,11) or (8,12) by (1,2) gives gaps 15 and 16.
The exact replacement is the crossed product

    (a,b) tensor (c,e) = (ae+bc, ac+be).                 (D1)

Together with componentwise addition, this makes N_0^2 a commutative semiring
with zero (0,0) and multiplicative identity (0,1). Multiplication by (1,0)
exchanges the coordinates. Its difference projection satisfies

    gap(P+Q)=gap(P)+gap(Q),
    gap(P tensor Q)=gap(P)gap(Q).                       (D2)

**Proof and retained coordinate.** Put s=a+b and d=b−a. This is an injective
map to pairs (s,d) satisfying s>=|d| and s=d mod2, with inverse
((s−d)/2,(s+d)/2). Addition is componentwise in these coordinates, and (D1)
becomes

    (s,d) tensor (u,e)=(su,de).

Associativity and distributivity follow immediately; the stated cone and
parity condition are preserved. Therefore equal-gap pairs are an arithmetic
congruence, and the quotient is exactly the ordinary integer ring Z.
The additive inverse exists in the quotient because exchanging coordinates
negates the gap; it is not an additive inverse of a raw nonnegative pair.

This is an operation on equivalence classes. It is not a claim that every
raw product image fills the whole target fibre: D_0 tensor D_0 consists of
pairs (2ac,2ac), so it misses (1,1), although every product has zero gap.

The lost translation has an especially simple exact memory law. Write a raw
pair as (d,t), meaning p(d)+(t,t). Then

    (d,t)+(e,u) = (d+e, t+u+kappa(d,e)),
    kappa(d,e)=(|d|+|e|−|d+e|)/2,

    (d,t) tensor (e,u) = (de, |d|u+|e|t+2tu).            (D3)

The first correction is cancellation between opposite signs. The second
follows from (|d|+2t)(|e|+2u)=|de|+2t_new. These formulas recover the complete
raw sum or product, not a chronological history of how it was produced.
They are the difference-family counterpart of the earlier
[decorated Laurent carry law](golden_primes_carries_route_compiler_20261003.md),
section 2, with a different value map and different kernel.

## 3. Precisely which endpoint operations forget the origin?

**PROVED classification.** A function F:Z^2->Z has the property that

    F(a+d,b+e)−F(a,b)

depends only on (d,e), for every a,b,d,e in Z, iff

    F(a,b)=alpha*a+beta*b+gamma

for fixed integers alpha,beta,gamma. No polynomial assumption is needed.

For the proof, subtract F(0,0) and set the base point to zero. The resulting
function H satisfies H(x+y)=H(x)+H(y) on Z^2. Decomposing an integer vector
along the two standard basis vectors gives H(a,b)=alpha*a+beta*b. The converse
is substitution. Thus nonlinear endpoint operations generally need an origin
coordinate or a changed lift such as (D1).

For ordinary endpoint multiplication, the exact missing information is visible
in the sum/gap chart:

    gap((a,b) component-product (c,e))
       = ((a+b)(e−c)+(c+e)(b−a))/2.                     (D4)

The two sums suffice to restore the product difference, subject to their
parity constraints. This is an explicit sufficient coordinate repair, not a
claim that two integers are minimal under arbitrary encodings.

The inherited boxplus law is affine and its constant cancels in endpoint
differences. The inherited boxtimes law is nonlinear and retains the mixed
terms from (D4). The four operation names alone therefore do not imply a
common quotient arithmetic; the evaluation map decides what descends.

## 4. The owner's triangular identity is a lossless family chart

Let T(n)=n(n+1)/2. The explicit proposed equality is correct:

    A+B=N−2,
    T(N−1)=T(A)+T(B)+(A+1)(B+1).                       (T1)

Take 0<=A<B and N>=3, and set d=B−A. The exact domain and inverse are

    0<d<=N−2,       d=N mod2,
    A=(N−2−d)/2,    B=(N−2+d)/2.

There are floor((N−1)/2) such nonnegative splits. If both inputs must instead
be positive, omit the A=0 boundary member. The two parts of (T1) become

    rectangle=(A+1)(B+1)=(N^2−d^2)/4,
    triangles=T(A)+T(B)=(N^2−2N+d^2)/4.                 (T2)

Their sum is T(N−1), proving (T1). The quantities N and d retain the whole
split, while N alone labels its finite family. If A is increased by one and
B decreased by one, the gap becomes d−2; the rectangle gains exactly d−1
and the triangular part loses exactly d−1. This remains inside the strict
split family when d>2.

The unshifted companion is T(A+B)=T(A)+T(B)+AB. The shifted rectangle in
(T1) is required when the triangular argument is A+B+1. These established
identities can support a definition of f, but do not determine the grouping
or meaning of the owner's string `(A+B)f1fAfB` by themselves.

## 5. Exact interpretations of the g stabilization equation

With c=AB, the literal scalar equation is x=cx. Over Z or Q it means
(c−1)x=0, hence x=0 unless c=1. In particular, for positive A,B it permits
an arbitrary scalar only at A=B=1. A nonzero realization must change the
object or equality relation explicitly.

**Modular repair.** Modulo a specified positive integer m, the fixed residues
are exactly

    x multiple of m/gcd(m,c−1),

and there are gcd(m,c−1) of them. Indeed divide (c−1)x=0 modm by that gcd
and cancel the remaining coprime coefficient. Taking m=c leaves only zero.
Taking m=c−1, for c>1, makes every residue fixed; c=2 gives the zero ring.

More intrinsically, Z/(c−1) is the universal ring quotient forcing c=1:
every unital ring map out of Z with image(c)=1 factors uniquely through it.
The proof is exactly that its kernel contains c−1. A faithful map retaining
all integer values cannot satisfy this relation when c!=1.

For any commutative ring R and element u, the same statement is R/(u−1).
The golden counterpart is O/(phi^L−1), where O=Z[phi]: it identifies shifts
by L digit positions. This is the precise common operation-minus-identity
construction behind scalar stabilization and the inherited golden monodromy
quotients. It does not identify their moduli, their exact values, or their
orbit predicates. When c=AB varies with the pair, its quotient ring also
varies and must be part of the object's type.

**Family repair and boundary.** For an integer c>1 there is no nonempty set
S of positive integers with cS=S. A least member r would require r/c in S,
contradicting minimality. For positive rationals there are invariant sets:
they are precisely unions of two-sided orbits r*c^Z. Equality implies closure
under both multiplication and division by c, which proves the classification.

The rooted forward ray R={r*c^k:k>=0} instead satisfies

    cR=R without its root r.

Its shift is a bijection onto a different tail, not equality of sets of
integer values. Keeping the root and the nonnegative height k permits a
route back to r; quotienting by the two-sided scale relation discards that
height and boundary. A finite truncation has two visible boundaries: it
loses r and gains r*c^(L+1). The script checks both, avoiding a spurious
finite-set invariance claim.

## 6. Dynamical difference fibres are guarded subsets, not full diagonals

For positive odd A, let S(A)=oddpart(3A+1) and K(A)=(A−S(A))/2. This is the
fully accelerated map. A positive rise of gap 2r, r>=1, is necessarily the
unique edge

    (A,S(A))=(4r−1,6r−1).                               (C1)

Indeed any valuation v>=2 gives S(A)<A for A>1; valuation one gives the
displayed formula. Thus D_4 contains the unique positive odd Syracuse rise
(7,11). Its translate (8,12) remains in D_4 but is not an edge of the odd map.
The quotient preserved the difference while forgetting the legal-edge test.

The inherited THM-4530 supplies the complete fibre over a drop d in Z, on
all signed odd sources:

    {8d+1, −4d−1}
      union {2d+(6d+1)/(2^v−3) : v>=3, (2^v−3) divides 6d+1}.       (C2)

Its mechanism is 6d+1=(2^v−3)S(A). Each displayed quotient is odd, so the
valuation is exact and source reconstruction A=S(A)+2d is reversible.
Only finitely many v occur because 2^v−3<=|6d+1|. For |d|<=D all displayed
sources lie in |A|<=8D+1, which makes the script's independent brute census
complete for its stated fibre universe.

A drop fibre does not itself support ordinary source multiplication:

    K(13)=K(33)=4,
    K(3*13)=−10,       K(3*33)=−25.

THM-4530 proves that retaining the next drop, (K(A),K(S(A))), uniquely
determines A on positive odd sources, and separately on the positive 3A−1
sheet. Here the two examples have signatures (4,2) and (4,3). One may decode
by listing (C2), keeping the positive candidates, and checking the second
drop. The signed domains must not be merged without their sheet/sign:
the signed + map has both A=1 and A=−1 with signature (0,0).

Two drops are therefore a proved source decoder on the stated domain, not a
descent proof or a ring homomorphism. Once the source is recovered one may
replay its unique next step. The whole diagonal family, its arithmetic
quotient class, and its finite legal-edge intersection are different objects.

## 7. Connection map, test universe, and next operation

| Source | Map and preserved predicate | Information lost or retained separately |
| --- | --- | --- |
| Nonnegative pair | Its gap d; crossed arithmetic gives ordinary Z | Translation t, needed to recover the original pair |
| Fixed triangular total N | Admissible gap d and formulas (T2) | N alone forgets which split was chosen |
| Scalar multiplication by c | Quotient by c−1 makes the action identity | Exact height, quotient and the pair-dependent modulus |
| Positive odd Syracuse source | One drop gives its finite fibre; two decode it | One drop loses the source; signed domains need their tag |

The script uses exact integer arithmetic and explicit failures, including
under Python -O:

    python 04-computation/experiments/difference_family_operations_20261004.py

Saved [output](difference_family_operations_20261004.out): 256 pair and 4096
triple semiring controls; 14450 decorated sum/product checks for gaps −8..8
and translations 0..4; 2401 affine controls; all 9900 strict nonnegative
triangular splits for N=3..200; 3321 scalar controls, 1088 modular fixed-set
checks, and 432 finite ray-boundary controls; all 2001 signed drop fibres
with |d|<=1000; and two-drop decoding of every positive odd source <=20001.
No primality filter or random sampling is applied.

The useful next step is to give each proposed operation its explicit value
map and retain a minimal *specified* missing coordinate: translation for
raw difference pairs, the modulus for stabilization, and the legal source
or a complete decoder for dynamical fibres. This lets a whole family remain
an integer object without pretending that every predicate descends to it.
