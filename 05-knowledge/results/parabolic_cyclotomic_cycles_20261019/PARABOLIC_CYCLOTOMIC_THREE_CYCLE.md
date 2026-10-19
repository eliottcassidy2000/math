# The parabolic three-cycle, the real seventh cyclotomic field, and an adelic boundary

2026-10-19. Status: PROVED elementary identities and consequences of the cited
rational-cycle denominator theorem. These are scoped connections, not a new
resolution of the Flynn–Poonen–Schaefer conjecture or a priority claim.

## 1. Inheritance and source

The rational marked three-cycle parameterization is inherited from
THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md and
appears in Morton–Raianu, *Arithmetic properties of 3-cycles of quadratic
maps over Q*, equations (3)–(4): https://arxiv.org/pdf/2105.07435 . That paper
also identifies a seventh-cyclotomic norm in the parameter numerator and
uses other cubic norm equations. The following particular identities are
proved directly; no conjecture about its numerator graph is imported.

Put

    D=2t(t+1),
    x0=(t^3+2t^2+t+1)/D,
    x1=(t^3-t-1)/D,
    x2=-(t^3+2t^2+3t+1)/D,
    c=-(t^6+2t^5+4t^4+8t^3+9t^2+4t+1)/(4t^2(t+1)^2).

For rational t!=0,-1 these are distinct rational points and x0->x1->x2->x0
under x^2+c. Distinctness follows, for example, from the nonzero rational
factor t^2+t+1 in their differences. Let lambda=8x0x1x2 be the multiplier.

## 2. Exact square-root unfolding at c=-7/4

Define g(t)=t^3+t^2-2t-1. Literal expansion gives

    c+7/4 = -g(t)^2/[4t^2(t+1)^2],
    lambda-1 = -(t^2+t+1)^3 g(t)/[t^3(t+1)^3].          (1)

The cubic g has discriminant49 and roots 2cos(2pi j/7), j=1,2,3. One
self-contained verification sets theta=zeta+zeta^(-1), with zeta a primitive
seventh root. Expanding g(theta) and multiplying by zeta^3 gives the sum of
all seven powers of zeta. The three real conjugates are distinct, and the
rational-root test makes g irreducible over Q.

In particular g has no rational zero, so every rational marked three-cycle
has c<-7/4. The parameter can approach -7/4 arbitrarily closely through
rational three-cycles, by taking rational t approaching any root of g.

A convenient cycle-level derivation retains the trace T=x0+x1+x2. Direct
symmetric elimination gives

    T^2+T+c+2=0,
    lambda=8[c(T+1)+1].                                (2)

The two three-cycle multipliers therefore satisfy

    lambda^2-(8c+16)lambda+64(c^3+2c^2+c+1)=0.

For c=-7/4-epsilon with epsilon>=0, this becomes the exact formula

    lambda_±=1-4epsilon ± (14+8epsilon)sqrt(epsilon).    (3)

The sign only labels the two cycle branches. Thus the multiplier splits
linearly in sqrt(epsilon), while the rational marked-parameter map in (1)
is quadratically ramified at each real seventh-cyclotomic root.

For a fully rational attracting control, take t=5/4. Then

    c=-226801/129600,
    lambda=5605019/5832000,     0<lambda<1.

Infinitely many rational t on the same side sufficiently near a root of g
give distinct rational parameters with attracting three-cycles and
lambda approaching1. Therefore there is no uniform positive real gap
between rational three-cycle multipliers and the parabolic value1.

## 3. Rational multiplier is not rational pointwise periodicity

For a chosen trace T with c=-T^2-T-2, the three cycle points are roots of

    G_T(X)=X^3-T X^2-(T^2+2T+3)X+T^3+2T^2+3T+1.        (4)

Its discriminant is

    disc G_T=(4T^2+6T+9)^2.                             (5)

The roots are permuted by X->X^2+c. A rational trace, multiplier, or square
discriminant therefore need not give three rational points: an irreducible
cyclic cubic can remain.

The exact hostile is c=-2. Its trace values are T=-1 and0; the two cubics
are X^3+X^2-2X-1 and X^3-3X+1, both irreducible over Q by the rational-root
test. Their multipliers are8 and-8, both rational. These are algebraic
three-cycles, not rational three-cycles. In particular the necessary square
condition -4c-7 being a rational square is not sufficient for rational
pointwise periodicity.

## 4. An actual inverse-orientation map at the parabolic parameter

Let theta range over the three roots of g. Squaring the underlying seventh
root gives the Galois permutation

    sigma(theta)=theta^2-2.

It is the restriction of the Chebyshev map f_-2 to these three points.
In Q[theta]/(g),

    sigma^2(theta)=-theta^2-theta+1.

Now set L(theta)=-theta-1/2. A literal identity modulo g gives

    f_(-7/4)(L(theta)) = L(sigma^2(theta)).              (6)

Thus the irrational parabolic three-cycle at c=-7/4 is an affine image of
the real seventh-cyclotomic three-cycle, with the Galois direction reversed.
Its point polynomial is

    8X^3+4X^2-18X-1,

and its multiplier is1. This is a bijection of the specified three-point
sets, not a conjugacy of the two quadratic maps on the entire real line.
It supplies a genuine order-three/cyclotomic connection rather than an
identification based only on the number7.

## 5. Rational neutral cycles are completely excluded beyond fixed points

Use the denominator result of Eliahou–Fares, *Some results on the
Flynn–Poonen–Schaefer conjecture*, Proposition2.1 and Theorem2.12 in arXiv v1:
https://arxiv.org/pdf/1903.08865 . Any rational cycle of exact length n>=3
has x_i=a_i/b, all gcd(a_i,b)=1, and 4|b. Consequently its multiplier obeys

    v_2(product_i 2x_i)=n(1-v_2(b))<0.                 (7)

At every odd prime p|b its valuation is -n v_p(b). Thus such a cycle is
repelling at every denominator prime, even if it is attracting over R.
In particular its rational multiplier cannot equal either1 or-1.

For an exact rational two-cycle u,v, subtracting the two equations gives
u+v=-1, and the multiplier is4uv. Multiplier1 forces u=v=-1/2, contrary
to exact period2. Multiplier-1 would give a quadratic with discriminant2,
so u,v would not be rational. For a fixed point the multiplier is2x.

Therefore the only rational periodic points, for rational c, with real
multiplier of absolute value1 are the fixed points

    x=1/2, c=1/4, lambda=1;
    x=-1/2, c=-3/4, lambda=-1.

This excludes every higher rational neutral cycle without assuming the
full period bound conjecture. It does not exclude hyperbolic rational
cycles of higher period. The attracting rational three-cycle above shows
why replacing this exact neutral obstruction by a uniform real gap would
be false.

## 6. Connection and remaining obstacle

Source: the marked rational three-cycle curve. Target: trace/multiplier
coordinates and a real cyclotomic boundary. Maps: symmetric elimination,
(1), and the three-point affine map (6). The symmetric formulas retain the multiplier and determine the
ramification order. The affine three-point map intertwines the stated
permutations but does not preserve derivatives or multipliers: the
cyclotomic c=-2 multiplier is 8 and the parabolic c=-7/4 multiplier is 1.
Lost by trace or multiplier: pointwise rationality of the cubic roots. Restoring sidecar:
the cubic splitting field, together with denominator valuations.

This is compatible with the earlier parabolic escape-clock work, but does
not turn a real bottleneck asymptotic into an all-denominator arithmetic
period bound. The scalar parabolic clock, its integral rounding error,
and the rational point-splitting condition are different objects.
