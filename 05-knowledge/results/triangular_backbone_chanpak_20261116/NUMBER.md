# From reflected sums to a square-index Pell subfamily

The two arithmetic identities (proved in [BALANCE.md](BALANCE.md))
come from one reflection with an unmatched center. For distinct a<b
write d=b−a and s=a+b:
$$
s=2a+d=2b-d,qquad ab=a^2+ad=b^2-bd=(s^2-d^2)/4.
$$
The last form is a useful inverse check for the summand and
multiplicand graphs: if one prescribes both sum and product, the
discriminant s²−4ab must be a nonnegative square with the right
parity (strict positivity for distinct vertices). Matching one
coordinate is not an inverse for the two-coordinate interface.

## Two Pythagorean rays meet once

The linear consecutive balances with length difference one produce
exactly the primitive consecutive-parameter triples
`(2n+1,2n(n+1),2n(n+1)+1)`, as proved in BALANCE.md. This
reuses the consecutive parameter ray in
[THM-3335](../../../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md).
The other natural ray is consecutive *legs* a,a+1 with hypotenuse c.
Put z=2a+1. The latter condition is z²−2c²=−1. These rays
meet exactly at (3,4,5), since
`2n(n+1)-(2n+1)=2n²-1=1` only at n=1 for positive n.

There is a clean two-way way to find the second ray inside the first
Pell system: **the positive square-triangular numbers whose triangular
index is itself a square are exactly the triangles T_(z²), where
z=2a+1 comes from a consecutive-leg Pythagorean triple**, also
allowing the degenerate first case a=0,c=1. Indeed,
$$
T_{z^2}=z^2(z^2+1)/2=(zc)^2
$$
if and only if (z²+1)/2 is an integer square. In the backward
direction a rational whose square is integral is integral; the
triangular-square index cannot be an even square without making
(z²+1)/2 a half-integer, which is not a rational square. Thus z
is odd. In Pell units this is exactly squaring
`z+c√2`, yielding `(2z²+1)+2zc√2`.
The (3,4,5) seed gives T_49=35²; (20,21,29) gives T_1681=1189².
By contrast T_8=36 is not on this square-index ray. The Pell
iteration and THM-3335 cover the generic square-triangular family,
not just this intersection.

## What the 8–9 chain is measuring

For leg-square pair s=a² and m=b² with a<b, the three
sums in THM-4510 are s,m,s+m=c². The composition of the
sum reflections is a translation by s or −s modulo m.
For (3,4,5) take s=9 and m=16, with exposed endpoints 8 and 9.
On the consecutive-leg Pell family a=b−1, there is an exact
two-scale arithmetic rule:
$$
s=(b-1)^2\equiv1-2b\pmod{b^2},\qquad
(s-1)^2\equiv0\pmod{b^2},\qquad bs\equiv b\pmod{b^2}.
$$
In coordinates x=u+bv, with 0≤u,v<b, one step of the
translation is u→u+1 mod b and v→v−2+[u=b−1] mod b.
The high-digit change is exactly one after a complete low-digit
turn: the −2 drifts cancel modulo b, the carry does not.
The exact moduli, the side of the
reflection, and the exposed endpoints matter. This is a model for
the tournament task: a phase that is invertible modulo a period
need not return at the prescribed absolute count.

## Small cutoffs versus transferred proofs

For Schur's problem two colors can avoid a monochromatic
x+y=z through 4 (put {1,4} together and {2,3} together).
With n=5, fix color(1)=A. Then 2=B,4=A. If 3=A,
1+3=4 fails; hence 3=B. Then 2+3=5 forces 5=A,
while 1+4=5 forces 5=B. This is the microcosm: the maximum
convention gives S(2)=4 and the first impossible prefix is 5.
In the [Boolean Pythagorean proof](https://arxiv.org/html/1605.00723)
a much larger exact forced-color clash plays the analogous
role; it is not a theorem that the two cutoffs scale together.
For tournament counts the productive analogy is the conditional
clause with two incompatible graph completions and a valid
local-to-global rule, developed in [CORE.md](CORE.md).
