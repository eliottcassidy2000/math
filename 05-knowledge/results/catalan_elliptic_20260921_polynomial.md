# The cubic polynomial: exact identities, an elliptic curve, and distinct symmetry actions

**2026-09-21. Status: PROVED under the stated reading of the pasted
formulas; CITED basic elliptic-curve inputs; FINITE-EXACT controls.**
The user's asterisks are ambiguous. A clarification was requested; until
answered, the following explicitly uses multiplication, not an inferred
exponentiation. No global Collatz result or complete integral-point
classification is asserted.

## Inheritance and concept board

The closest inherited mechanism is the
[shifted arithmetic model](arithmetic_seams_20260921_operations.md): first
specify an operation and its actual coordinate map. The
[rational-cycle lift](arithmetic_seams_20260921_dynamics.md) provides the
canonical hostile: the same abstract S_3 can act on different sets with
different preserved predicates. The corrected near miss is equating a
three-point symmetry with a global conjugacy. The least-used coordinate
here is the Galois action on the three nonzero two-torsion points.

| Live object | Exact predicate | Cheap decisive probe |
|---|---|---|
| Pasted b,c,s,t | Polynomial identity | Expand before imposing any curve equation |
| Cubic P(x) | Square values | Parity modulo8 and a nontrivial integer point |
| Two elliptic models | Rational or geometric equivalence | j-invariant, then good-prime point counts |
| Three roots | Galois symmetry | Irreducibility and discriminant squareclass |
| Rational group orbit | Integrality versus rationality | A small multiple acquires denominators |

Anchor: a lawful map to the user's elliptic curve. Niche: square-value
obstructions. Wildcard: distinguish the relevant S_3 actions explicitly.

## 1. The polynomial identities do simplify, with one retained difference

Over Q, interpret the input as

```text
b=x+1,
c=xb/2+b(b+1)/2,
s=2x^3+3x^2-x-1+c,
t=-(xb^2-x)-(bx^2-b).
```

Direct expansion proves

```text
c=(x+1)^2,
s=x(2x^2+4x+1),
t=-(2x+1)(x^2+x-1),
s+t=c.                                                (P1)
```

The last identity reflects `t=(x+b)(1-xb)` and `s=c-t`.
For the separately proposed polynomial `P(x)=2x^3+4x^2+1`,

```text
P(x)-s(x)=1-x.                                        (P2)
```

Thus P and s agree only at x=1, not as polynomials. The identity (P1)
is a valid parametrized additive relation; it does not itself put the
three values on the fruit cubic or impose any Collatz successor rule.

## 2. A natural square-value interpretation gives a different elliptic curve

If the intended connection asks when P(x) is a rational square, introduce
`y^2=P(x)`. The rational change `X=2x,Y=2y` gives

```text
E_0: Y^2=X^3+4X^2+4.                                  (P3)
```

This interpretation is explicit, rather than an assertion that the
user supplied the extra equation. Compare it with

```text
E_4: Y^2=X^3+109X^2+224X.
```

For `Y^2=X^3+aX^2+bX+c`, use
`c4=16(a^2-3b)` and elliptic discriminant equal to16 times the cubic
discriminant. Exact arithmetic yields

| Invariant | E_0 | E_4 |
|---|---|---|
| Elliptic discriminant | `-2^8*7*13` | `2^14*5*7^2*13^3` |
| j | `-65536/91` | `1408317602329/2153060` |
| Number of F_11 points, including infinity |16|12|

Both curves are nonsingular. Distinct j-invariants show they are not
isomorphic even over an algebraic closure of Q. More strongly, the
different point counts at the common good-reduction prime11 show they
are not Q-isogenous. This excludes even a nonconstant Q-morphism between
their smooth projective curves: translating the image of the origin
would turn such a morphism into an isogeny.

**CITED structural inputs:** j classifies geometric elliptic isomorphism;
isogenies preserve point counts at common good primes; a nonconstant
morphism of genus-one curves with origins is a translated isogeny. See
[Elkies, Math223 course notes, Chapters III, V and VII](https://people.math.harvard.edu/~elkies/M223.24/index.html).
For the finite-field count mechanism,
[Sutherland, Lecture7, Theorem7.3](https://math.mit.edu/classes/18.783/2025/LectureNotes7.pdf)
identifies `#E(F_q)` with `deg(Frobenius-1)`. An isogeny intertwines
Frobenius, so multiplicativity of degree makes these degrees equal.
The two explicit invariant calculations and point censuses are our exact
controls, not imported numerical database claims.

The shared primes7 and13 in discriminants are insufficient to construct
a map. This does not exclude useful correspondences through a higher
genus curve or a relationship between selected arithmetic predicates;
those would require their own map and hypotheses.

## 3. There is a real S_3, acting on the two-torsion rather than the fruit labels

The monic cubic `X^3+4X^2+4` has no root modulo5 and is therefore
irreducible over Q. Its discriminant is `-1456=-16*91`, not a rational
square. Its splitting-field Galois group is consequently S_3: an
irreducible cubic has a transitive subgroup of S_3, and the squareclass
of its discriminant detects containment in A_3. The quadratic resolvent
is `Q(sqrt(-91))`.

For a Weierstrass model the nonzero two-torsion points are `(r,0)`, where
r runs through those three cubic roots. Thus this is a concrete
Galois S_3 action on `E_0[2]` minus its origin. In particular E_0 has no
nonzero rational two-torsion point.

For E_4 the cubic is `X(X^2+109X+224)` and the quadratic discriminant is
`10985=13^2*65`. Its two-torsion splitting field is `Q(sqrt(65))`;
Galois fixes `(0,0)` and exchanges the other two points. That action is
C_2, not S_3. The separate S_3 that permutes the fruit coordinates is
an automorphism action on the whole fruit cubic, explained in the
[elliptic lane](catalan_elliptic_20260921_elliptic.md). This distinction
identifies exactly which three objects are being braided in each case.

## 4. A small rational point generates the nontrivial square49^2

On E_0 start with Q=(0,2). Exact chord-and-tangent addition gives

```text
 Q=(0,2),              2Q=(-4,-2),
3Q=(1,-3),             4Q=(20,98),
5Q=(-24/25,326/125).                                  (P4)
```

Returning through (P3), `4Q` gives

```text
P(10)=2*10^3+4*10^2+1=2401=49^2.
```

There are infinitely many rational points: E_0 has trivial rational
torsion. To see this, the common odd-good-prime reduction bound injects
rational torsion into the groups over F_3 and F_11, whose orders are6
and16. Its order divides2; the irreducible cubic excludes order2.
The reduction bound is the **CITED** input discussed in Chapter VII,
Application3.2 of the Elkies notes above. Therefore Q has infinite order.
No rank-one or generator-of-the-full-group assertion is needed here.

The denominators at5Q show why an elliptic group orbit cannot simply be
read as a sequence of integer values. Clearing denominators changes the
affine coordinates, and does not preserve a proposed integer successor
rule. This is a concrete height/integrality boundary for the analogy.

## 5. An adjacent-square obstruction survives without a global Catalan claim

**PROVED.** Every integer point on `y^2=P(x)` has even x. For odd x,
`P(x)` is3 or7 modulo8, neither a square. Put x=2z and y=2k+1. Then

```text
k(k+1)=4z^2(z+1).                                     (P5)
```

This is a genuine factorization into adjacent integers, with
`gcd(k,k+1)=1`. If z>0 and z+1=q^2 for an integer q, then
`P(2z)=(4zq)^2+1`, strictly between two consecutive squares. Thus this
entire infinite family of x is excluded. It is an elementary
adjacent-square obstruction, and does not require Mihailescu's much
more general perfect-power theorem.

For x<=-3, `P(x)=2x^2(x+2)+1<0`, so there are no real y. The exact
census for `-2<=x<=200000`, y>=0, finds only `(-2,1),(0,1),(10,49)`.
**FINITE-EXACT:** this last sentence is a bounded census. Completeness
of integral points at all heights is not proved here. The full
Diophantine problem can be pursued as an integral-point computation
with a certified height bound or a suitable descent; a finite scan
cannot replace that obligation.

## Reproduction

```text
python 04-computation/experiments/catalan_elliptic_20260921_polynomial.py
python -O 04-computation/experiments/catalan_elliptic_20260921_polynomial.py
```

The [script](../../04-computation/experiments/catalan_elliptic_20260921_polynomial.py)
and [JSON](../../04-computation/experiments/catalan_elliptic_20260921_polynomial.json)
give explicit universes, exact fractional group operations, independent
finite-field counts, and optimized-safe checks. Polynomial expansion
and factorization prove the all-parameter identities; rational grids
only corroborate them. No numerical approximation is used to decide a
square, equality, discriminant, or point count.
