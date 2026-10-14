# Complete dyadic prime-cofactor three-chamber classification

Research successor, 14 October 2026. This completes the affine-template conjecture proposed in [the prime-chamber note](aliquot_prime_chamber_cycles_20261014.md), using the [unbounded repeated-minimum theorem](aliquot_repeated_minimum_classification_20261014.md) and the two additional arguments below. The earlier packets remain unchanged. Global novelty has not been established.

## The theorem and its scope

For positive integers `a_0,a_1,a_2`, put `t_i=2^a_i`. Consider the cyclic affine equations

`t_(i+1) p_(i+1)=(t_i-1)p_i+2t_i-1`, for `i mod 3`.

Their unique rational solution is integral if and only if, up to cyclic rotation,

1. `a_0=a_1=a_2=a`, with all cofactors `p_i=2^(a+1)-1`; or
2. `(a_0,a_1,a_2)=(1,1,4)`, with cofactors `(p_0,p_1,p_2)=(23,13,1)`.

This is an unbounded classification: no limit on the exponents is assumed.

**Aliquot consequence.** There is no aliquot cycle of three distinct integers all of whose members are even and have prime odd part, that is, each member is `2^a p` with `a>=1` and `p` an odd prime.

Indeed, those factorization hypotheses make the affine equations exactly the true aliquot equations. The first case gives only a constant itinerary, and is an actual perfect number only when its displayed Mersenne cofactor is prime. The second gives the virtual loop

`46 -> 26 -> 16 -> 46`,

whose last cofactor is one rather than a prime. The actual last step is `s(16)=15`, because `sigma(1)=1`, not `1+1`.

This does **not** rule out general sociable 3-cycles. Odd members, composite odd parts, and odd-prime powers with exponent greater than one lie outside the exact prime-cofactor model. The integral classification is a statement about the affine equations; it does not turn those equations into the true divisor-sum function on composite cofactors.

## 1. Affine closure and the dyadic lattice

For a dyadic core `t=2^a`, the formal branch on `n=t p` is

`n -> ((t-1)/t)n+2t-1`.

Every slope lies strictly between zero and one, and every intercept is positive. Composition of three branches therefore has a unique positive rational fixed point. Thus the cyclic system has a unique positive rational solution.

For cores `x,y,z`, its common odd denominator is

`D=xy+yz+zx-x-y-z+1=xyz-(x-1)(y-1)(z-1)>0`.

The first cofactor numerator is

`N(x,y,z)=(z-1)(y-1)(2x-1)+(z-1)y(2y-1)+yz(2z-1)`.

The cofactors are respectively `N(x,y,z)/D`, `N(y,z,x)/D`, and `N(z,x,y)/D`.

For every integral solution, each cofactor is odd: reduce its outgoing affine equation modulo two. A stronger restriction comes from the smallest core `x`. Write all cores as `x u_i`. Reducing the equations modulo `x` gives

`p_i=-1 mod x`.

Hence `H_i=(p_i+1)/x` is a positive integer. Substitution and division by `x` give

`u_(i+1)(xH_(i+1)-1)=(x u_i-1)H_i+u_i`.

Reducing this equation modulo `x` yields

`H_i=u_i+u_(i+1) mod x`.                                  (1)

Equivalently, `p_i=x(u_i+u_(i+1))-1 mod x^2`. The successor position matters: reversing a three-chamber itinerary changes this congruence at the largest core.

## 2. Repeated minimum: the previously proved case

If the smallest exponent occurs twice, rotate the itinerary to `(a,a,c)`, `c>=a`. The [repeated-minimum theorem](aliquot_repeated_minimum_classification_20261014.md) proves that its integral closure is constant unless `(a,c)=(1,4)`, which gives `(23,13,1)`.

Its exact arithmetic gate is

`(2x-1)z+(x-1)^2 | (2x^2-1)(3x^2-3x+1)`,

for `x=2^a,z=2^c`. The unbounded proof uses the quotient modulo `x^2` and dyadic ratio bounds. This section invokes that proved result, not its finite census.

It remains to exclude a repeated maximum and three distinct exponents.

## 3. Repeated maximum: an interval too short for odd cofactors

Rotate the cores to `(x,y,y)`, with `2<=x<y`. Let `q,r` be the cofactors at the two equal cores. From the exact numerator formula,

`q-r=(x-y)(2x+2y-1)/D`,

where `D=y^2+2xy-x-2y+1`.

The difference is negative. Moreover,

`2D-(y-x)(2x+2y-1)=2x^2+4xy-3x-3y+2>0`.

For positivity, the right side is `y(4x-3)+(2x^2-3x+2)`, and both summands are positive for `x>=2`. Consequently

`0<r-q<2`.

But integral cofactors are odd, so their difference must be an even integer. This is impossible. There is no integral closure with a repeated maximum and a distinct minimum.

## 4. Three distinct cores: the two cyclic orders

Write the sorted cores as

`x < xu < xv`,

where `x>=2`, `u>=2`, and `v>=2u` are powers of two. Up to cyclic rotation, there are exactly two orders:

- ascending: `x -> xu -> xv -> x`;
- reverse: `x -> xv -> xu -> x`.

Put

`S=1+u+v`, `T=uv+u+v`, `D=x^2T-xS+1>0`.

At the largest core `xv`, the exact value `H=(p_largest+1)/x` is `A/D` in ascending order and `R/D` in reverse order, where

`A=2x^2uS-x(u+1)(v+2)+v+1`,

`R=2x^2uS-x(uv+v+2u^2+2u)+u+v`.

The decisive interval is

`1 < R/D < A/D < 3`.                                      (2)

It holds throughout the real parameter region `x>=2,u>=2,v>=2u`, so in particular for the dyadic parameters. Here is an explicit proof.

### Lower bound and comparison

Set `B=2u^2+u+(u-1)v`. Then `B>v`, and the identity

`R-D=x(x-1)B-x(v-1)+u+v-1`

is strictly positive, since `x>=2` and `B>v`. Also,

`A-R=(u-1)(2x(u+1)-1)>0`.

These establish `1<R/D<A/D`.

### Upper bound

Define

`K=x^2(u+3)+x(u-2)-1>0`,

`F=14x^2-3x-2+(u-2)[7x^2-2+x(2u-1)]>0`.

For the latter inequality, `14x^2-3x-2=(x-2)(14x+25)+48`, and the remaining bracket is positive. Direct expansion gives

`3D-A=(v-2u)K+F>0`.

This finishes the proof of (2).

## 5. The congruences exclude both orders

### Ascending order

The successor of the largest core has ratio one. Equation (1) forces

`H_largest=v+1 mod x`.

Both `x` and `v` are even, so `H_largest` is an odd integer. The strict interval `1<H_largest<3` contains no odd integer. Ascending order cannot close integrally.

### Reverse order

Here the successor of the largest core has ratio `u`, so

`H_largest=v+u mod x`.

It is an even integer. Equation (2) therefore forces `H_largest=2`, and the largest-core cofactor is `p_largest=2x-1`.

The next affine equation, from the largest core to the middle core, becomes

`u p_middle=v(2x+1)-2`.

Since `u|v`, integrality of `p_middle` forces `u|2`. Thus `u=2`.

But when `u=2`, the exact numerator difference is

`R-2D=(8x^2-6x)-v(2x^2+x-1)`
`      =(4-10x)+(v-4)(1-x-2x^2)<0`,

because `x>=2` and `v>=4`. This contradicts `R/D=2`. Reverse order cannot close integrally either.

## 6. Completion and interpretation

All exponent multisets have now been covered: constant, repeated minimum, repeated maximum, or three distinct exponents. The repeated-minimum theorem supplies the sole nonconstant integral exception, and that exception fails the prime requirement. This proves both the full affine classification and the restricted aliquot no-three-cycle theorem.

The useful mechanism is an arithmetic refinement of ordering: the two cyclic orders place the largest-core quantity in the same open interval `(1,3)`, but force it onto different integer parity classes. One order misses the lattice immediately; the other is forced to its unique admissible point, whose successor divisibility then contradicts its exact value. Merely recording a tournament orientation or an abundance sign would lose those congruence and prime-factorization conditions.

## 7. Exact verification

The companion dependency-free script verifies thirteen displayed identities in `Z[x,u,v]`, including both largest-core formulas, the interval decompositions, the repeated-maximum difference, and the reverse-order boundary contradiction. It also expands `R-D`, `A-R`, and `3D-A` after the substitution `x=X+2,u=U+2,v=2(U+2)+V`: every nonzero coefficient is positive, and the constant terms are respectively `27,11,48`. The JSON records these coefficient lists, providing exact certificates of strict positivity for all `X,U,V>=0`. It independently composes the affine maps with rational arithmetic over the already-tested small exponent window `1,...,12`, covering both orientations. No enlarged census is needed for the theorem.

Run:

`python -B 04-computation/experiments/aliquot_dyadic_three_chamber_classification_20261014.py`

Normal and optimized output agree. The proof is unbounded; the finite checks only corroborate it.
