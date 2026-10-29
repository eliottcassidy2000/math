# Fruit symmetry quotient: a cubic reconstruction guard, not a free bit

Session97, 2026-10-29. **PROVED** algebra below; **FINITE-EXACT** checks in
`check.py`; **CITED** historical rank/generator statements. No global
least-height certificate, rank computation, Sha computation, or Collatz theorem.

## Recovery rather than rediscovery

Pinned source commit: `5d6eb91dd76eee86b91e04d856daa32c32b3d61d`.
The relevant prior paths are
`05-knowledge/results/catalan_elliptic_20260921_elliptic.md`,
`05-knowledge/results/prime_shells_20260921_elliptic.md`, and
`05-knowledge/results/collatz_mod6_20260921_fruit_rank_positive_multiples.md`.
They already establish the six torsion guards, pair-sum reciprocation,
ternary index tree and failure of positivity heredity. In particular the
last note distinguishes finite saturation from a full generator proof,
and repairs a reversed BSD index factor. Those cautions remain in force.

The primary [Bremner–MacLeod paper, §2 and Remark2.2](https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_43_from29to41.pdf)
was opened directly this session: it supplies the curve/map, generator
`G=(-4,28)`, signed solution `(11,4,-1)`, and ninth-multiple positive
example. Its Lemma2.6 supplies the degree-three isogeny identified below.
The [2025 corrigendum record](https://publikacio.uni-eszterhazy.hu/8862/)
is inherited from the source note, not newly checked; no table-wide
leastness assertion is imported. Rank and first positive height are
separate quantities: this rank-one example already has 79–81-digit
positive coordinates. It supplies no bound21 for other curves.

## 1. Original rational expression and its exact exceptional set

Write `s=a+b+c`, `q0=ab+ac+bc`, `r0=abc`. Clearing denominators gives

    s³ − 6s q0 + 7r0 = 0,

**equivalent** to the original sum being4 only when every pair sum is
nonzero. The excluded projective points are precisely the six permutations
of `[1:-1:0]` and `[1:1:-1]` (three of each). Indeed imposing `a+b=0`
in the cubic gives `c(c²−a²)=0`; symmetry gives
these six points. A valid point cannot have `s=0`: then `abc=0`, giving
one of the excluded points. Rational solutions scale to primitive integer
triples by a common denominator and gcd; choose `s>0` for a unique sign.
A signed solution is not a positive one.

For `E: y²=x³+109x²+224x`, use

    [a:b:c]=[56−x+y:56−x−y:−56−12x],
    x=−28(a+b+2c)/(6a+6b−c),
    y=364(a−b)/(6a+6b−c).

The inverse denominator vanishes only at the invalid infinity point on
the cubic: substitute `c=6(a+b)` into the equation to obtain a nonzero
constant times `(a+b)³`, then `a+b=c=0`. The curve's six rational torsion
points `O,(0,0),(4,±52),(56,±728)` map to the six excluded points.
Everything else is an exact rational/integer solution. The point at
infinity maps to `[1:-1:0]` by extending the map, not by plugging infinite
numbers into affine formulas.

## 2. The symmetry quotient and the missing cubic guard

Normalize `s=1` and write `q=q0/s²`, `r=r0/s³`. Then `r=(6q−1)/7`.
The unordered triple is the root multiset of

    f_q(t)=t³−t²+qt−(6q−1)/7.                    (1)

Thus **q is a complete invariant of an unordered triple on this curve**.
It does not follow that every rational q reconstructs a rational triple.
The discriminant calculation is

    49 Disc(f_q)=(q+1)(1+29q−196q²).

For the oriented alternating invariant
`W=7(a−b)(b−c)(c−a)/s³`, this gives the intermediate cover

    W²=(q+1)(1+29q−196q²).                       (2)

The six-sheeted generic ordered-to-unordered map factors as a
three-sheeted cyclic-permutation quotient followed by a two-sheeted
orientation quotient. Equations(1),(2) distinguish **geometric degree**
from **rational reconstruction**: over Q, a point on (2) need not lift.

**Explicit hostile.** `(q,W)=(0,1)` satisfies (2), but (1) becomes
`7t³−7t²+1=0`. Its only possible rational roots are `±1,±1/7`, none of
which is a root. It is irreducible. A discriminant square has certified
cyclic Galois symmetry here, not three rational roots. Adding an orientation
bit cannot repair this arithmetic failure.

**Exact restoration guard.** For rational `(q,W)` away from the
repeated-root divisor, it is enough to retain one rational root of(1).
The remaining quadratic has square discriminant: the cubic discriminant
is the quadratic discriminant times the square of the product of the
root's differences from the other two. Therefore it splits over Q.
After this root witness, orientation chooses the other two labels.
For merely unordered reconstruction, one rational root plus the rational
square discriminant is sufficient; named-port reconstruction still owes
the chosen root and orientation. This is an actual field-theoretic
instance of the repository's faithful-carrier warning.

## 3. The cover is the inherited three-isogeny, explicitly

Set `Z=−196(q+1)`, `Y=−196W`. Equation(2) becomes

    Y²=Z³+421Z²+43904Z.                          (3)

From the original Weierstrass point the map is

    Z=x(x−56)²/(x−4)²,
    Y=(x−56)(x²+44x+224)y/(x−4)³.              (4)

These identities follow by substituting the fruit coordinates and
expanding. They also match Bremner–MacLeod Lemma2.6's degree-three
isogeny at N=4. Its kernel is `O,(4,52),(4,−52)`, exactly the cyclic
coordinate-permutation subgroup. Thus the cubic guard is precisely a
rational isogeny-lifting obligation, not an invented analogy. The hostile
`(-196,-196)` on(3) does not lie in this isogeny's rational image, because
such a lift would yield rational roots of its irreducible cubic.
Nonlifting here is not a proof of a nontrivial Sha element: local
solubility of the associated torsor has not been established.

## 4. Positivity is an exact quotient chamber—but rationality is extra

For a split real cubic with sum of roots1, all roots are positive iff
`q>0` and `r>0`: if two roots are `−u,−v` with u,v>0, the third is
`1+u+v` and q is negative. A negative product already excludes three
positive roots. Consequently the positive chamber is

    1/6 < q < (29+5 sqrt(65))/392,               (5)

together with rational splitting of(1). To prove the upper bound use
(2), `q>1/6`, and the roots of `1+29q−196q²`. Equality at the upper
endpoint is irrational and never occurs for a rational triple. This
makes positivity cheap to decide *after* reconstruction, not a shortcut
past rational splitting. The interval endpoints and the number7 come
from this cubic; no tournament `{7,21}` identification follows.

## 5. Exact point recurrence, scope and further direction

On E, for P=(x,y), Q=(u,v), the secant slope is `(v−y)/(u−x)`, or
`(3x²+218x+224)/(2y)` for doubling. The sum has coordinates
`z=m²−109−x−u`, `m(x−z)−y`; inverse pairs give O. This generates all
multiples and six torsion cosets exactly without floating point.
`check.py` verifies 120 such points, both maps, the cubic, the quotient,
the chamber, and the hostile irreducibility; it reproduces 9G's primitive
79/80/81-digit triple. First positive among the tested multiples is9,
not an independently certified globally least integer triple.

Run `python3 check.py` and `python3 -O check.py`. Checks use explicit
exceptions and standard-library Fractions. The next useful extension is
an exact cubic-splitting/isogeny-image descent certificate: distinguish
local obstruction from genuinely everywhere-locally-soluble failure.
That could connect to the Sha session through actual torsors rather than
shared integers, but it is OPEN here.

### Second checkpoint: a positive-real, locally impossible lift

The initial hostile q=0 lies outside the positive chamber. A stronger
FINITE-EXACT hostile is the point `35H` on(3), where `H=(-196,196)`.
Its q coordinate lies strictly in(5), and its W is rational, so all three
roots of(1) are distinct positive real numbers. Nevertheless q has even
numerator and odd denominator. Reducing the monic cubic over Z_2 gives
`t³+t²+1`, which has no F_2 root; thus there is no Q_2 root and no
rational fruit lift. A root of a monic polynomial over Z_2 must be
integral (otherwise its leading term has uniquely lowest valuation),
justifying the reduction obstruction. The checker reconstructs35H
from its tiny generator rather than embedding hundreds of digits.
This separates all three obligations: real positivity, rational
orientation, and local arithmetic lifting. It is explicitly **not** a
Sha witness, since it fails already over Q_2. The q=0 hostile has the
same mod2 obstruction.


## 6. An all-order ternary arithmetic return controller

**PROVED**, not an inference from the tested multiples:

    nH belongs to phi(E(Q))  iff  3 divides n,           (6)

where phi is(4), restricted to the cyclic subgroup generated by H.
Exact addition gives `3H=(-225,210)=phi(G)`. Thus the class of H in the
abelian quotient `E'(Q)/phi(E(Q))` has order dividing3. It is not zero:
the q=0 rootless-mod2 certificate excludes H from the rational image.
Its order is therefore exactly3, proving(6) for every integer n.
This proves neither that the full cokernel has order3 nor the rank or
generator status of E'. No such extra claim is needed.

Consequently the named ternary transitions `n -> 3n+r`, r in
`{-1,0,1}`, have a genuine local lifting controller: precisely the
r=0 child has a rational phi-lift. For that child the explicit source is
`nG`, since `phi(nG)=3nH`; retain the integer n, not merely its residue,
to reconstruct that source and its height. The three lifts differ by the
three kernel points, corresponding to the three cyclic coordinate
rotations. An ordered root owner chooses one of them. The endpoint
residue does not provide the omitted source/height or ensure positivity.
In particular `27H=phi(9G)` lifts to the positive example, whereas35H is
positive-real but not rationally liftable. This is an arithmetic
return-controller connection with the repository's ternary grammar,
not an identification with Collatz or a proof of descent.

### Positive-real nonlifting is infinite, not an isolated hostile

Using the CITED infinite order of G, the identity `3H=phi(G)` implies
that H has infinite order. On a real elliptic curve the identity component
is a circle, and every infinite cyclic subgroup of that circle is dense:
a closed proper subgroup of a circle is finite. The point6H lies on that
component and has infinite order. Therefore `{35H+6kH:k in Z}` is dense
in the component containing35H. The chamber(5) is open and contains35H,
so infinitely many members of this progression lie in that chamber.
Every index is `2 mod3`, hence none of these points lifts rationally by(6).
Thus even real positivity together with rational orientation leaves
infinitely many arithmetic false positives. This uses standard real
elliptic-curve topology (the background source is the inherited
[Elkies course](https://people.math.harvard.edu/~elkies/M223.24/index.html));
no effective bound on the next positive index or distribution error is
claimed. The finite35H computation supplies a nonempty open target, not
the infinitude proof.
