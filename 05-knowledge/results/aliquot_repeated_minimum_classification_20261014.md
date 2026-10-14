# Complete classification of repeated-minimum dyadic three-chamber closures

Research successor, 14 October 2026. This extends the proved subfamily in [the prime-chamber note](aliquot_prime_chamber_cycles_20261014.md). That earlier packet remains unchanged. The classification below is proved for the repeated-minimum family, with no exponent cutoff. Its global novelty has not been established. The unrestricted existence of aliquot 3-cycles remains open.

## The theorem

Let `a>=1`, `c>=a`, `x=2^a`, and `z=2^c`. Consider the affine prime-cofactor template

`x q=(x-1)p+2x-1`,
`z r=(x-1)q+2x-1`,
`x p=(z-1)r+2z-1`.

Its unique rational solution has all three cofactors integral if and only if

1. `c=a`, with `p=q=r=2^(a+1)-1`; or
2. `(a,c)=(1,4)`, with `(p,q,r)=(23,13,1)`.

The first case is a constant itinerary. In the second case the integer itinerary is the virtual loop `46 -> 26 -> 16 -> 46`; cofactor `r=1` is not prime, and the actual aliquot successor of `16` is `15`.

**Aliquot consequence.** If a three-cycle of distinct aliquot numbers has each member of the form `2^a p` with `p` an odd prime, then its minimum 2-adic exponent occurs exactly once. All dyadic prime-cofactor itineraries with a repeated minimum are excluded. This does not exclude the remaining exponent patterns with a unique minimum, nor numbers outside these prime-cofactor chambers.

## 1. The exact divisor reduction

Put

`D=(2x-1)z+(x-1)^2`,
`Q=(2x^2-1)(3x^2-3x+1)=6x^4-6x^3-x^2+3x-1`.

The prior note proves that integral closure is equivalent to `D|Q`. For completeness, its closure numerator is

`N=4x^2 z-4x^2+2xz^2-5xz+4x+z-1`,

and

`(2x-1)^2 N+xQ`
`=D(6x^3+4x^2 z-10x^2-2xz+5x-1)`.

Because `z` is divisible by `x`, `D` is coprime to both `x` and `2x-1`. Thus `p=N/D` is integral exactly when `D|Q`. Integrality of the other cofactors then follows from the three affine equations, as detailed in the prior proof.

When `c=a`, direct substitution gives the constant solution. Henceforth assume `c>a`, and write

`z=x t`, where `t=2^(c-a)>=2`.

Then

`D=(2x^2-x)t+(x-1)^2`.

Assume `D|Q`, and introduce the positive integer quotient `K=Q/D`.

## 2. A second-order dyadic congruence

Modulo `x^2`,

`D=1-(t+2)x`, and `Q=3x-1`.

Since the inverse of `1-(t+2)x` modulo `x^2` is `1+(t+2)x`,

`K=-1+(1-t)x mod x^2`.

Consequently

`H=(K+1)/x`

is a positive integer satisfying

`H=1-t mod x`.                                             (1)

As a positive rational expression before imposing integrality,

`H=(Q+D)/(xD)`.                                            (2)

The congruence and a small upper bound will now exhaust all dyadic cases.

## 3. The large-ratio sector: t>=x

Here `t` is a multiple of `x`, because both are powers of two. Equation (1) gives `H=1 mod x`.

The denominator `D` increases with `t`, and

`((3x-1)D-Q)|_(t=x)=x(4x^2-5x+2)>0`.

Therefore `(2)` gives `0<H<3` throughout `t>=x`. With `x>=2` and `H=1 mod x`, necessarily `H=1`, hence `K=x-1`.

But `Q=1 mod (x-1)`. The equality `Q=DK=D(x-1)` forces `x-1` to divide one, so `x=2`. Now `Q=49`, `D=6t+1`, and `K=1`, hence `t=8`. This is precisely `(a,c)=(1,4)`.

The stronger estimate `H<3x/t` is not used in this sector and is not being asserted globally.

## 4. The small-ratio sector: 2<=t<x

Because `x,t` are powers of two, `x>=2t`. The positive solutions of (1) begin with

`H>=x+1-t`.                                                (3)

In this sector one has

`H<3x/t`.                                                 (4)

To prove (4), the numerator of the required inequality is exactly

`3x^2 D-t(Q+D)`
`=x[-(2x-1)t^2+(3x^2-1)t+3x(x-1)^2]`.

For `0<t<=x`, its bracket is at least

`(x^2+x-1)t+3x(x-1)^2>0`.

Thus (4) is valid throughout the needed sector.

### Case t>=8

Using `x>=2t`,

`x+1-t-3x/t = x(1-3/t)-t+1 >= t-5>0`.

This contradicts (3) and (4).

### Case t=4

Equations (3) and (4) require `x-3<3x/4`, hence `x<12`. The only power of two satisfying `x>=8` is `x=8`. At this point,

`D=529`, `Q=21463=40*529+303`,

so `D` does not divide `Q`.

### Case t=2

Here `x>=4`. The positive congruence solutions in (1) are `x-1,2x-1,...`; the bound `H<3x/2` leaves only `H=x-1`. Thus

`K=xH-1=x^2-x-1`.

But direct subtraction gives

`Q-D(x^2-x-1)=x^2(x^2+3x-1)>0`,

contradicting `Q=DK`.

The three cases exhaust the positive dyadic ratios `2<=t<x`. Together with Section 3, this proves the classification.

## 5. Why this is more than a larger search

The finite census in the earlier note suggested a classification of all dyadic three-chamber itineraries. The present proof settles an infinite two-parameter subfamily by combining:

- exact affine closure;
- an integer divisor gate;
- the quotient modulo `x^2`;
- positivity and the ordering of powers of two.

The only surviving nonconstant algebraic closure fails the original prime-factorization hypothesis. The result thus preserves the distinction between a closed affine itinerary and a genuine aliquot cycle.

The remaining part of the proposed full dyadic classification has a unique minimum exponent. No assertion about that part follows merely from this proof.

## 6. Reproducibility

The companion dependency-free script checks the displayed polynomial identities in an exact two-variable polynomial ring, verifies the quotient congruence algebraically, and checks the case inequalities and divisor classification over the same repeated-minimum window used by the earlier packet (`1<=a<=64`, `a<=c<=3a+2`). The finite checks are corroboration; Sections 2–4 are the unbounded proof.

Run:

`python -B 04-computation/experiments/aliquot_repeated_minimum_classification_20261014.py`

Normal and optimized runs produce identical JSON, and all verification checks remain active under `-O`.
