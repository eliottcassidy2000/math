# A complete quadratic-field atlas for the exceptional norm conic

2026-10-26. **PROVED deduction from the cited quaternion splitting theorem
and inherited resolved-fibre equations; independent audit PASS.**
The [independent audit](../certificates/AUDIT_QUADRATIC_FIBRE_ATLAS.md)
reviews the pre-promotion proof hash and includes a separate 117-gate checker;
this status-and-link amendment changes no mathematical statement.
The local-global quaternion theorem is classical. The contribution is its
exact application to every rational pair (D,delta) in the existing norm-seven
fibre family, with an explicit decision rule and minimum lifting degree.

## 1. Inheritance, carrier, and the new question

The closest proved mechanism is
[RESOLVED_FIBRES, R1–R4](../../normseven_descent_moments_20261026/RESOLVED_FIBRES.md).
For nonzero rational D,delta it uses the branch order

    (E,E+x,E-x,E+z,E-z,E+x+z,E-x-z)

and twists

    (2/delta,1,1,1,-D,1,-1).

For either marked complementary norm-seven pair, the quotient product has
a rational point y over [E:x:z]=[0:0:1]. Normalization has eight geometric
nodes over y. If D is nonsquare, these are four degree-two points with
residue field K=Q(sqrt(D)); if D is square, all eight are rational.
Resolving each node inserts

    C_delta: A^2+B^2=delta C^2.                          (1)

For nonsquare D, the complete reduced resolved fibre is four copies of this
K-conic regarded as a Q-scheme by composition with Spec(K)->Spec(Q).
For square D it is eight Q-conics. This is not Weil restriction.

That source proves examples with minimum point-field degrees two and four,
and identical good-prime fibre zeta functions. It does not classify every
quadratic K. We retain its geometry and close precisely that arithmetic
decision problem here.

Canonical hostile: delta=3,D=17 passes the earlier prime-three obstruction
test but fails at two. Corrected near miss: listing only odd ramified primes
misses a genuine obstruction. MISTAKE-506 separately warns that denominator
support is not a quaternion ramification set; our set is computed from the
actual Hilbert symbol. Least-used sidecar: the Brauer class of the tangent
conic over the node's residual field.

Board: residual quadratic field; transverse squareclass; local ramification;
minimum point-field degree; good-prime counts; relative split torification.
Anchor: decide the resolution lift. Niche: local-global conic splitting.
Wildcard: the prime-two and real-place gates invisible to a first odd-prime
probe. None of these statements classifies interior hourglass points.

## 2. An explicit finite list of local gates

Replace delta by its unique signed squarefree integer representative q in
Q*/Q*^2. Write q=2^e u with e in {0,1} and u odd, possibly negative. Define

    R(delta) = {odd primes p dividing q with p=3 mod4}
               union {2 if u=3 mod4}
               union {infinity if q<0}.                 (2)

Replace D similarly by signed squarefree d. If d!=1, define

    Split(d) contains an odd prime p iff p does not divide d
                     and d is a nonzero square modulo p;
    Split(d) contains 2 iff d=1 mod8;
    Split(d) contains infinity iff d>0.                  (3)

These are exactly the rational places that split in Q(sqrt(d)). In
particular a prime ramified in that quadratic field is not split.
Only the finitely many places in R(delta) need testing.

**Theorem QA (resolved-fibre atlas).** Let D,delta be nonzero rationals.

1. If d!=1 and K=Q(sqrt(D)), then

       C_delta(K) is nonempty iff R(delta) intersect Split(d) is empty.

2. Let m(D,delta) be the smallest degree [F:Q] of a number field giving
   a point of the specified complete resolved fibre above y. Then

       d=1:   m=1 if R(delta) is empty, and m=2 otherwise;
       d!=1:  m=2 if the intersection in part 1 is empty,
                        and m=4 otherwise.              (4)

3. Over the actual node field K (K=Q when d=1), each exceptional conic
   has a finite decomposition into split K-tori iff it has a K-point,
   equivalently iff the corresponding test above passes.

This is a minimum-degree statement for a fixed rational point y and fixed
twists. It is not a statement about the least degree of any point anywhere
on the whole surface, or about integral lifts.

## 3. Proof and exact classical inputs

Put B_delta=(-1,delta)_Q. A direct coordinate change relates (1) to the
conic of pure quaternions: its pure reduced norm is

    X^2-delta Y^2-delta Z^2,

and [A:B:C] -> [delta C:A:B] takes (1) to its zero locus.
Thus over any characteristic-zero field extension F, C_delta(F) is
nonempty iff B_delta tensor F is split. The equivalence is the classical
conic/split-quaternion theorem, [Voight, Theorem 5.5.3](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_5).

For odd p the Hilbert symbol of (-1,delta) is

    (-1,delta)_p=(-1)^((p-1)v_p(delta)/2).

At two it is (-1)^((u-1)/2), because (-1,2)_2=1 and the symbol on odd
units has that value. At the real place it is negative exactly when
delta<0. These give precisely (2). They follow from the explicitly checked
local formulas in [Voight, §12.4.12–16](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_12).
The symbol is unchanged by rational squares; factors in denominators
contribute through valuation parity, not through their mere presence.

The classical quadratic splitting criterion says that a quadratic number
field K splits a quaternion algebra B over Q iff no ramified place of B
splits in K. This is [Voight, Proposition 14.6.7](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_14).
Its hypotheses apply because K/Q is separable quadratic. At a ramified
place, a degree-two local extension splits the local division quaternion;
a split place leaves two copies of the original obstruction. The global
local-global theorem then gives equivalence, not just necessity. The
elementary quadratic-field splitting tests are (3): odd-prime Hensel
lifting, the odd-unit square test modulo eight at two, and the sign at
infinity. This proves part 1. For d=1 the same global splitting theorem
says B_delta splits over Q iff R(delta) is empty.

An F-point of the complete fibre must choose a node-field embedding K->F
and then a point of its conic, by the inherited fibre description. If the
conic splits over K, the least possible degree is [K:Q]. If it does not,
no F of that degree works. In all cases adjoining i suffices:
[A:B:C]=[1:i:0] lies on (1). In the nonsplit case i is not already in K,
so [K(i):K]=2. Degrees of fields containing K are multiples of [K:Q].
This proves (4), including exactness of each lower and upper bound.

The infinity point above is legitimate: this question asks for points of
the complete exceptional conic. It is not an affine norm equation with
C required nonzero. A soluble smooth conic is P1 over its field, and has
the split-torus decomposition G_m plus two points. Conversely every
nonempty split torus has a rational identity point, so an anisotropic
conic has no such finite decomposition. This proves part 3.

The local-global input is CITED, not re-proved by the finite checker.
The application to the inherited fibre and the minimum-degree consequence
are the deductions here. No Brauer-Manin obstruction for the whole root
surface, and no description of its full Brauer group, is asserted.

## 4. Delta three: a closed congruence rule and two different hostile places

For delta=3, R(3)={2,3}. Hence, for signed squarefree d!=1,

    C_3(Q(sqrt(d))) has a point
        iff d != 1 mod8 and d != 1 mod3.                 (5)

The second condition includes multiples of three, which are ramified and
therefore permitted. No positivity condition on d is needed here because
the real place is not ramified in (-1,3).

| d | Gate that survives over K | Minimum fibre point degree |
|---:|---|---:|
| 2 | neither | 2 |
| 3 | neither | 2 |
| 5 | neither | 2 |
| 13 | prime 3 | 4 |
| 17 | prime 2 | 4 |
| 73 | primes 2 and 3 | 4 |
| -1 | neither | 2 |
| -7 | prime 2 | 4 |

For d=2 use [1:sqrt(2):1]; for d=3 use [sqrt(3):0:1]; for d=5 use
[(1+sqrt(5))/2:(1-sqrt(5))/2:1]. For d=-1 the point [1:i:0] suffices.
These are exact positive controls.

For d=13 the old obstruction occurs at three. For d=17, the field embeds
in Q_2 because 17 is a 2-adic square, but (1) with delta=3 has no nonzero
Q_2 solution. Scale to integral coordinates with at least one unit. Modulo
four, odd C would demand a sum of two squares congruent to three; even C
forces both A and B even. Both cases are impossible. Thus a test looking
only at primes congruent to three modulo four is insufficient over K.

The real place is independently necessary: delta=-1 over any real
quadratic field fails, whereas K=Q(i) supplies [1:i:0]. For example d=2
has no finite obstruction left but still fails at infinity. Conversely
delta=2 has no ramified place, so every node field passes. The ramification
set of (-1,-1) is {2,infinity}, matching the corrected MISTAKE-506 boundary.

## 5. What the counts and F1 language forget

For fixed nonsquare D, the inherited good-prime fibre zeta formulas are
independent of nonzero delta once the bad primes are excluded. The atlas
shows exactly what these counts omit: the intersection of the conic's
ramification places with the splitting places of its residual field.
Equal smooth finite-field conic counts are compatible with different
minimum rational lifting degrees and relative torifiability.

Over Q the complete fibre for nonsquare D has no rational point for any
delta, so none is a finite union of nonempty split Q-tori. A passing
relative test over K must not be promoted to a split Q-torification.
These are local closed-fibre assertions, not whole-surface F1 models.

The incoming [Gale result](../../normseven_arithmetic_f1_20261026/gale/RESULT.md)
uses a quadratic Galois character on NS to distinguish constant K3 twists.
The present atlas uses a different invariant, a quaternion class on an
exceptional conic. They are compatible examples of retaining arithmetic
data after geometric splitting, not an identification of their carriers.

Source -> target: fixed quotient point and root twists -> node field plus
tangent conic -> finite list of local splitting gates. Preserved: the
specified fibre, field of definition, and exact minimum lifting degree.
Lost by good-prime counts: bad-place and real-place norm obstructions.
Required sidecars: D and delta squareclasses, and the fact that the conic
is projective. Cheapest decisive test: delta=3,D=17 modulo four.

## 6. Reproduction and remaining scope

Run `python3 -B check_quadratic_atlas.py` and
`python3 -B -O check_quadratic_atlas.py`; the saved JSON receipts must be
byte-identical. The standard-library checker compares Hilbert predictions
with direct primitive congruence searches, computes quadratic splitting by
an independent root search, verifies exact quadratic-field point witnesses,
and records all minimum-degree decisions in its declared finite universe.
Negative bounded searches are used only as explicit modular obstructions,
never as evidence of global absence without the proved criterion.

The local modular controls themselves are decisive in their declared
squarefree-integral scope. At odd p, a primitive solution modulo p^2 has
a unit partial derivative: if p divides delta and both A,B were divisible
by p, then C would be a unit and the equation would fail modulo p^2.
Otherwise the unit derivative is immediate. Hensel lifting applies.
At two, a primitive solution modulo eight has a partial derivative of
valuation one; for even squarefree delta, the equation modulo four forces
at least one of A,B odd. The one-variable Hensel inequality 3>2*1 then
lifts the solution. Conversely any p-adic point scales to a primitive
integral one. These controls check the cited Hilbert calculation through
actual equations, not through a second implementation of its formula.

The final receipts have **16,773 explicit gates**. Their SHA256 is
`296406290b1003a43bc4feae0a6e927f63121dc967a2e3620338ddd526e9c5ec`.

Remaining questions include nonconstant K3 twists, integral saturation and
arithmetic of the entire resolved quotient map. This packet neither solves
those questions nor claims an implication for RH or integer hourglasses.
