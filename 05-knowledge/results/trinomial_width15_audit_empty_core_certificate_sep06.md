# Independent audit of the lower-carry trinomial family

**Status: INDEPENDENT ANALYTIC + FINITE-EXACT IDENTITY AUDIT PASS.**
The candidate family is

```text
(-15,2g-15,3g-15),     g>=8 integral, gcd(g,15)=1.
```

For arbitrary nonzero complex coefficients, its first nonzero constant
term is at mass `g` or `2g`. At either first-row cancellation parameter,
the doubled moment divided by the square of the first anchor monomial is
strictly negative. This proves an analytic family beyond the old
smaller-endpoint-at-most-eight strip; it does not close all supports of
smaller endpoint fifteen or the general trinomial problem.

The full producer and proof belong to the
[returns-lane note](trinomial_width15_empty_core_returns_sep06.md).
This audit imports neither that producer nor any symbolic algebra package.
Its [source](../../04-computation/trinomial_width15_audit_empty_core_certificate_sep06.py)
and [output](trinomial_width15_audit_empty_core_certificate_sep06.out) use
standard-library exact rational arithmetic.

## 1. The actual carrier and first-return scope

The support gcd is `gcd(15,g)`. At any support return of mass `m`, the
charge equation implies `g|15m`. Under the stated gcd condition, `g|m`.
At mass `g`, the complete multiplicity channels are

```text
(g-7+j,6-3j,1+2j), j=0,1,2.
```

At mass `2g`, they are

```text
(2g-15+j,15-3j,2j), j=0,...,5.
```

All counts are nonnegative for `g>=8`, with the negative-charge count
strictly positive. These lists follow directly from `2y+3z=15` and
`2y+3z=30`. They certify first support return exactly `g`, first degree two,
doubled degree five and lower carry one. Nonprimitive integral parameters
may still be used to check polynomial identities, but not to infer first
support return `g`.

For `u=g-5`, the normalized first polynomial and the complete canonical
second polynomial are

```text
P=6tau²+20u*tau+u(u-1),
Q=sum_(j=0)^5 (2g)_(15-j)*tau^j/[(15-3j)!(2j)!].
```

The first positive content is `(g)_5/720`. The second canonical anchor is
the square of the first anchor times `tau^-1`. Hence the sign target is
`R=tau^-1 Q mod P`. Since `u>=3`, `P` has positive constant term and
discriminant `8u(47u+3)>0`; both first roots are distinct and negative.

## 2. A separate polynomial-identity certificate

Put

```text
K=(2g)_10/32
 =g(g-1)(g-2)(g-3)(g-4)(2g-1)(2g-3)(2g-5)(2g-7)(2g-9)>0.
```

Reduce powers using the two exact relations

```text
tau²=-(10u/3)tau-u(u-1)/6,
tau^-1=-(6tau+20u)/[u(u-1)].
```

There is no parameter pole after removing `K`: the inverse-term coefficient
satisfies the literal polynomial cancellation

```text
[tau^0]Q/[K*u(u-1)]
 =128(2g-11)(2g-13)(2g-14)/15!.
```

For `j>=1`, `[tau^j]Q/K` is a polynomial of degree `5-j`. The remainder
of `tau^(j-1)` has constant degree at most `j-1` and linear degree at most
`j-2`. Therefore, if `R/K=Cbar+Dbar*tau`, then

```text
deg Cbar<=4, deg Dbar<=3,
deg Tr(R)/K<=4, deg Norm(R)/K²<=8.
```

The independent producer reconstructs these rational values from the two
recurrences at the nine distinct parameters `g=8,...,16`. It matches all
candidate constant, linear, trace and norm formulas exactly. The displayed
degree bounds turn these nine evaluations into identities for every
parameter; they are not an empirical extension from nine samples.

The trace and norm formulas are

```text
Tr(R)=-K(g-5)J(g)/18389170800,
Norm(R)=K²(g-5)H(g)/3757351141239696000000,
```

where `J` has degree three and `H` degree seven. An independent binomial
shift reconstructs every coefficient in `J(s+8)` and `H(s+8)` from their
unshifted arrays; all coefficients are strictly positive. Thus for all
real `g>=8`, trace is negative and norm is positive. The two root values
are real, have negative sum and positive product, and are therefore both
strictly negative. This is the quadratic case of the
[carry-corrected trace criterion](trinomial_trace_sign_empty_core_certificate_sep06.md).

## 3. Consequence, direct controls and reproducibility

At a first-row zero, the nonzero anchor and `R<0` prove that the doubled
constant term is nonzero. All intermediate masses have empty support by
the divisibility argument. Consequently the first nonzero mass is exactly
`g` or `2g`, both attainable: all coefficients one give mass `g`, while
choosing `alpha=tau`, `beta=gamma=1` realizes either nonzero first root and
gives mass `2g`. In particular `2g<3g=15+(3g-15)`.

Independent literal Laurent multiplication checks every earlier empty row
and both complete weighted rows at `g=8,11,13,14`. These are fixed controls
of the carrier interpretation, separate from the nine evaluations used to
certify the polynomial identities. The `g=8` anchored remainder is
`4743488+47087024tau`, recovering the known carry-sign hostile exactly.

```sh
python3 -B 04-computation/trinomial_width15_audit_empty_core_certificate_sep06.py
python3 -B -O 04-computation/trinomial_width15_audit_empty_core_certificate_sep06.py
```

All 69 explicit gates pass; normal and optimized outputs are byte-identical.
Semantic digest:
`a267f91a0e714f37758141ff7432f178e0fae2ae86693f1083856126a48e1bb5`.
Raw SHA-256:

```text
source ca0a5416a4e2c92d590fc3e701a9ac5d3a26ded1535d81d839eeffd2d834d83a
output 45c71b8626b1efc1aa3b97b58b18ba0482591db0c0ace7fd3bb7df1f6579c622
```
