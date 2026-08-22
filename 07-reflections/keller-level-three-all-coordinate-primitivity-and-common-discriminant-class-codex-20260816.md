# The three coordinate views still share one sign class at degree 27

**Status: VERIFIED-EXACT theorem candidate for the fixed sporadic map at
level three, pending independent proof audit.**  The exact witness is
`04-computation/keller_level_three_three_coordinate_primitive_finite_field_probe_20260816.py`.
This note does not classify the outside counterexample family and has no
`JC(2)` or LRC consequence.

## 1. Inheritance and the gate that was actually missing

At level one, proved THM-2546 gives the three cubic coordinate discriminants

```text
Delta_x=-4 A_x^2 L,
Delta_y=-4 A_y^2 L,
Delta_z=-4 A_z^2 L.
```

Their common square class is `[-L]`.  At level two, proved THM-3508 first
establishes that all three actual source coordinates are primitive in the
generic degree-nine extension and then uses trace-form congruence to put all
three discriminants in `[H]`.  The canonical hostile is an intermediate
coordinate: its minimal degree is only three, and its apparent degree-nine
eliminant is a cube with zero discriminant.

Proved THM-3495 supplies the corresponding level-three facts only for the
`x` view:

```text
[K_3:K]=27,                 [Disc(m_x)]=[-2J].         (1)
```

The missing gate was not another discriminant expansion.  It was the
primitivity of the actual source coordinates `y` and `z` in the same
degree-27 algebra.

## 2. The first lawful target is a useful hostile

THM-3495's independent finite-field audit uses the target

```text
(a,b,c)=(93,28,83) in F_101^3.                         (2)
```

It splits through two levels into nine second points and proves the product
of their nine `x` cores squarefree of degree 27.  Reconstructing `y` and `z`
inside those same cubic quotients gives derivative-gcd degrees

```text
(x,y,z)=(0,1,0).                                       (3)
```

Thus the target in (2) proves `x` and `z` separation but not `y` separation:
one pair of geometric sheets collides under the `y` projection.  This does
not refute generic `y` primitivity.  It refutes the shortcut “one lawful
degree-27 fibre that works for `x` automatically works for every coordinate.”

## 3. A simultaneous all-coordinate witness

Continue the parent audit's deterministic permutation of `F_101^3`.  After
5,937 targets, the exact search reaches

```text
(a,b,c)=(77,62,4).                                      (4)
```

Its first fibre is

```text
(13,36,2), (39,84,75), (49,74,79),                    (5)
```

and the three blocks of second points are

```text
((35,40,46),(83,54,93),(84,61,87)),
(( 2,91,44),(35, 7,50),(64,28,37)),
((24, 2,74),(84,50,58),(94,59, 2)).                   (6)
```

Every point in (4)--(6) lies off `V(L)`, every inverse denominator used in
the two split stages is nonzero, and every final cubic core is separable.
Thus the specialized fibre algebra is the honest rank-27 etale algebra

```text
A=A_1 x ... x A_9,             dim_F101 A_i=3.         (7)
```

No assumption that the last cubics split over `F_101` is made.

## 4. Exact quotient-algebra computation

For a second point `(a,b,c)`, let

```text
g(X)=L X^3+(4-3bc)X-2c,
A_i=F_101[X]/(g).                                      (8)
```

Inside (8), the companion reconstructs the remaining source coordinates by
the exact inverse formulas

```text
D=((12a-b^2)X^2+bX+2),
y=b-3aX((9ac-b)X+2)/D,
z=(2X-c-3X^2y)/X^3.                                   (9)
```

It proves `D` and `X` are units by extended gcd, substitutes `(X,y,z)` into
all three displayed coordinates of the original map, and recovers the
second point exactly in the quotient.  This is stronger than evaluating a
chosen rational root.

For each of `X,y,z`, form the `3 x 3` multiplication matrix in the basis
`(1,X,X^2)`.  Its characteristic polynomial is computed from trace,
`trace(M^2)`, and determinant, and Cayley--Hamilton is checked back inside
the quotient.  All 27 local coordinate polynomials have degree three and are
squarefree.

Multiplying the nine local polynomials coordinate by coordinate gives three
monic degree-27 polynomials.  Their derivative gcds are all one, with ordered
coefficient hashes

```text
x: 306dc5831a989a0fe953ae5eb6127d26bd6ff07f8cc4e25d54dc0f65e53e5e19
y: aba63853d168e670461e67b67e42d26380e18c62b363c8931aa187f1a97d4e51
z: b8a09c5bfed7affdabbbb3fcf222002a3a7c11c8d55cb19bcc949494725bb51d. (10)
```

Hence every coordinate has 27 distinct eigenvalues on the geometric fibre
and generates the full algebra (7).

## 5. Generic primitivity and the common class

On the generic open set where the three cubic leading coefficients and the
inverse denominators are units, use any fixed rank-27 tower basis.  For
`xi in {x,y,z}`, the determinant of the power matrix

```text
(1,xi,xi^2,...,xi^26)                                  (11)
```

is a rational function of the final target.  The good-reduction witness
(4)--(10) is a specialization where this determinant is nonzero for all
three coordinates.  Therefore none of the three generic determinants is
identically zero, and

```text
K(x)=K(y)=K(z)=K_3.                                    (12)
```

For a primitive coordinate, its power-basis discriminant differs from the
trace discriminant of the fixed tower basis by the square of the basis-change
determinant.  Combining (12) with THM-3495's exact full class (1) gives the
candidate theorem

```text
[Disc(m_x)]=[Disc(m_y)]=[Disc(m_z)]=[-2J]              (13)
```

in the generic target function field modulo squares.  As in THM-3508, any
actual saturated degree-27 eliminant differs from its monic minimal
polynomial by a scalar whose discriminant exponent is `2*27-2=52`, an even
power.  Thus (13) is unchanged by denominator clearing or nonzero target
scaling.

## 6. What this says about the three `-4` cubics

The constant `-4` at level one is not three independent classifications.
It is a square times the sign unit, so after quotienting by squares all three
cubics expose the same quadratic permutation character.  The fixed tower now
has the following candidate atlas:

```text
level 1, degree  3:  all three classes [-L]       (PROVED THM-2546),
level 2, degree  9:  all three classes [ H]       (PROVED THM-3508),
level 3, degree 27:  all three classes [-2J]      (this exact candidate).
```

The carrier changes under composition: old components cancel in parity and
the sign quotient selects the newest image prime, with a constant unit that
must not be discarded.  This is a dynamical discriminant pattern for one
fixed map, not a classification of all counterexamples.

In particular, equality of (13) retains only the common `C2` sign quotient
of three degree-27 root actions.  It forgets labelled roots, the full wreath
monodromy, block systems, exact discriminant multiplicities, the other
Jelonek components `L,H`, boundary valuations, and effective sections.  The
outside infinite family has not been shown to share (12), (13), or the same
Jelonek tower.  Its odd fibre count `1+2k` is an involution-orbit statement,
not a discriminant-class theorem.

## 7. Reproduction and boundary

Run

```text
python -B 04-computation/keller_level_three_three_coordinate_primitive_finite_field_probe_20260816.py
python -B -O 04-computation/keller_level_three_three_coordinate_primitive_finite_field_probe_20260816.py
```

The two transcripts agree byte-for-byte.  The semantic digest is

```text
61ddcb9dbd3c7582a514bf1a42bbff3cd9ca7d052ea9620721979ed6c8138245.
```

Independent audit must still verify the quotient inverse, map substitution,
characteristic-polynomial calculation, good-reduction specialization, and
trace-form inference before canon promotion.  Even after that audit, no
exact positive discriminant multiplicities, level-four three-coordinate
primitivity, all-level induction, outside-family classification, arbitrary
Keller theorem, `JC(2)`, or LRC result follows.
