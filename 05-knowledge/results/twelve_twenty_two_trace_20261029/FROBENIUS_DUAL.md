# Seventh-root traces: exterior dual, Frobenius, and Paley orientation

**PROVED scoped.** The exterior conjugacy and root identification are
inherited from Session104. The clock and orientation maps below use that
specific algebra and fixed chart; they are not maps between arbitrary
tournament atoms.

## F1. Three distinct uses of squaring

For the counter C in [the README](README.md), let D=C-I. Its characteristic
polynomial is q(d)=d^3+d^2-2d-1. In the ordered wedge basis(01,02,12),

\[
\Lambda^2 C=\begin{pmatrix}1&0&-1\\1&-1&-1\\0&-1&-1\end{pmatrix},
\quad S=\begin{pmatrix}0&-1&2\\-1&-2&4\\-1&-3&5\end{pmatrix},
\quad (\Lambda^2 C)S=SD,\quad\det S=1.
\]

Thus the synthetic exterior square of1 direct-sum C has the two integer
blocks C and D. This is the existing
[Session104 conjugacy](../session104_descent_exterior_20261029/exterior/K4_EXTERIOR_BOUNDARY.md),
not a new six-letter legality representation.

If zeta is a primitive seventh root of unity, the three roots of q are
zeta^k+zeta^(-k)=2cos(2pi k/7), for k=1,2,3. The order-three Galois
automorphism is

\[
\sigma(d)=d^2-2,\qquad\sigma(u)=u^2-2u\quad(u=1+d).
\]

The subtraction of2 matters: **angle doubling is d^2-2, not d^2**.
This operation permutes the conjugates. In contrast, the clock replacement
C to C^2 acts on the same state module; and state squaring x to x^2 is
nonlinear. The three operations must not be identified. P1's integer
formula recovers C from C^2 without making them Galois conjugates.

An older exact map uses this same cubic: the affine map L(d)=-d-1/2
intertwines sigma^2 with x to x^2-7/4 **on the three roots of q**.
See [the October19 parabolic three-cycle](../parabolic_cyclotomic_cycles_20261019/PARABOLIC_CYCLOTOMIC_THREE_CYCLE.md)
and its [independent audit](../parabolic_cyclotomic_cycles_20261019/INDEPENDENT_AUDIT.md).
It is not a conjugacy on neighborhoods: the polynomial-cycle multiplier
is1, whereas the angle-map cycle has multiplier8. The finite root set and
the derivative behavior carry different information.

## F2. The11 chart explains clocks12 and22

The discriminant is49, so reduction at a prime p other than7 is etale.
On A/pA the arithmetic Frobenius equals sigma^r, where

|p modulo7|r|
|---|---:|
|1 or6|0|
|2 or5|1|
|3 or4|2|

This follows by raising zeta to p and identifying exponents up to sign.
It is valid also for split primes, when A/pA is a product of fields.
Exact algebra gives u sigma(u)=d, hence u sigma^2(u)=sigma^2(d).
Consequently

\[
u^{p+1}=\begin{cases}
u^2&r=0,\\d&r=1,\\\sigma^2(d)&r=2,
\end{cases}\qquad
u^{2p}=(\sigma^r u)^2.
\]

At p=11, r=2 and the cubic is irreducible. Therefore in F1331,

\[
\boxed{u^{12}=\sigma^2(d),\qquad u^{22}=\sigma^2(u^2).}
\]

These are literal source-preserving equalities connecting12=11+1 and
22=2*11 to the exterior dual and square. The element u has order266;
each of u^12 and u^22 has order133=11^2+11+1 and generates the norm-one
subgroup. The checker verifies order266 by all its prime divisors2,7,19.

The degree is decisive. For the Pell unit v=1+sqrt(2), the analogous
11+1 power is the quadratic norm v^12=-1, collapsing its phase. For
the cubic u the full field norm uses1+11+11^2, not11+1, so u^12 instead
retains a nontrivial conjugate of d. Neither formula transfers merely
because the integer12 occurs in both constructions.

At p=13 the cubic splits, and u^12=1. This produces P4's latency
obstruction and P5's growing precision layers. At7 the etale/Frobenius
table is excluded; the nilpotent trace kernel must be retained explicitly.

## F3. A literal F8 map and its orientation character

Modulo2, the counter polynomial is x^3+x+1 and u generates F8^* of
order7. Fix the logarithm chart

\[
u^i\longleftrightarrow i\in\mathbb Z/7\mathbb Z.
\]

Define a tournament on these seven points by i to j iff j-i belongs to
{1,2,4}. Vertices are nonzero field elements; the observable is the
quadratic character of their logarithmic difference. The chosen generator
fixes the orientation gauge; there are no ties for distinct vertices.

Multiplication by u induces i to i+1. Field Frobenius induces i to2i
and preserves the arcs. Together these give the affine group of order21.
The exact checker compares all7! permutations with this group and
independently counts189 Hamilton paths, recovering the inherited Paley7
controls rather than claiming a new census.

But d=u-1=u+1=u^3 in characteristic2. The chart change from generator u
to generator d is multiplication by3 on exponents (or by its inverse5,
depending on the direction of chart identification). Both3 and5 are
quadratic nonresidues modulo7, so **every arc is reversed**. The synthetic
exterior dual therefore supplies an actual converse map in this fixed
F8 chart. Frobenius supplies an orientation-preserving transport.

The retained sidecar is the generator's square/nonsquare class, a single
orientation bit once this source and chart have been fixed. This does
not assert that one bit reconstructs arbitrary attachments, arbitrary
F8 state data, or the twelve-vertex81/27 pair. The prior Fano/octonion
construction needs an additional chosen Fano plane;14 cyclic triangles
are not automatically the seven lines of that plane. See
[the September22 Paley/Fano packet](../collatz_mod6_20260922_paley_fano_octonion_design.md).

## Transfer contract

The accepted chain is the specified integer counter algebra, its exterior
dual, reduction modulo2, a chosen logarithm chart, and the Paley arc
predicate. It preserves the stated operators and records when orientation
reverses. It forgets integer precision and the source grammar's individual
letters. The `ee` rank2 obstruction still forbids identifying the entire
legality machine with a four-dimensional exterior-square representation.
No Collatz orbit, graph homology, or all-size strong-atom classification
is transported by this chain without a further map.
