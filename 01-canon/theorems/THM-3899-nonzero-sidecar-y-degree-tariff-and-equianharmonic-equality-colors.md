---
id: THM-3899
title: "Nonzero-sidecar y-degree tariff and equianharmonic equality colors"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every
  nonzero-sidecar square survivor of the
  THM-3881 residual satisfies deg_y(T)>=deg_y(f).  On the positive common
  degree seam, its leading coefficients obey a two-color equianharmonic norm
  law over sqrt(-3), with an odd color parity at x+1.  This is a necessary
  leading filtration, not an existence theorem; all lower coefficients, the
  constant common-degree seam, a Keller atlas, and JC(2) remain OPEN.
source: jc_zero_debt_lift / post-THM-3897 first live f-nonzero gate, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED on 2026-08-23.  The audit independently
  enumerated all fifteen strict deficits, rechecked the two equality-top
  terms, both a-adic obstructions, the primewise v|g step, every boundary,
  and fourteen concrete hostile cells with nontrivial lower coefficients.
  The companion expands the full
  residual, certifies all fifteen strict degree deficits, isolates the two
  equality-top monomials, checks the norm/color identities and a positive
  leading-color payment, and identifies d=2lambda-1 in 28 active gates.
  Normal and optimized runs byte-match the frozen output on every platform.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
related:
  - THM-3897-f-zero-residual-all-degree-global-emptiness
  - THM-3898-equianharmonic-cube-cubic-order-a1-branch-euler-tariff
script: 04-computation/jc2_nonzero_sidecar_y_degree_tariff_thm3899.py
output: 05-knowledge/results/jc2_nonzero_sidecar_y_degree_tariff_thm3899.out
script_sha256: afa3cc0f1a158253d884151dd91bcb2b1fa5509644bdeb22004ad7d34a48982b
output_sha256: 18d6accc6f16f35c3582e23ad67ccbaa883c66cafbe000b05b32a63845068a14
semantic_sha256: 306a6a9a487506858d1461b95d35cec4e69bdd8d152a4d3a9465407fc6cb5af2
hash_basis: raw LF bytes
---

# THM-3899 -- a nonzero sidecar cannot lead in y-degree

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  Retain the THM-3881
coordinates

```text
D=k[x,y],                 a=x+1,
L=9x+4,                   K=y^2-15x^2-15x-4,
P=aL^2,
r=aT+Kf,                  A=KT+aPf.                       (1)
```

The complete residual is

```text
S(T,f)=L^4
 +2(3A+3P+r^2)L^2f
 +(8A+6P+3r^2)(Pf^2-T^2).                                (2)
```

Suppose `T,f,G in D`, `f!=0`, and retain the exact THM-3881 module address:

```text
S(T,f)=G^2,                  T(0,0)=4f(0,0).               (3)
```

Then `T!=0` by the proved `T=0` lane of THM-3881.  Put

```text
m=deg_y(T),                  n=deg_y(f).                   (4)
```

The first exact tariff is

```text
m>=n.                                                        (5)
```

If moreover `m=n>=1`, write `u,v,g in k[x]` for the respective leading
`y`-coefficients of `T,f,G`.  There is an `h in k[x]` such that

```text
g=vh,
h^2=3(aL^2v^2-u^2).                                       (6)
```

For either `d in k` with `d^2=-3`, this is the two-color norm

```text
(h-du)(h+du)=3aL^2v^2.                                    (7)
```

At the prime `a=x+1`, the sum of the two color valuations is odd:

```text
ord_a(h-du)+ord_a(h+du)=1+2ord_a(v).                      (8)
```

Thus one color has odd `a`-valuation and the other has even valuation.  If
`a` does not divide `v`, exactly one color contains the simple `a`-factor.
This is an intrinsic binary leading observable, but it is not sufficient for
the lower coefficients of `(3)`.

## 1. Exact residual support

Expanding `(2)` gives

```text
S=L^4+6aL^4f+12a^2L^4f^2+8a^3L^4f^3
 +6KL^2Tf+12aKL^2Tf^2+6a^2KL^2Tf^3
 -6aL^2T^2-6a^2L^2T^2f+3a^3L^2T^2f^2
 -8KT^3-6aKT^3f-3a^2T^4
 +2K^2L^2f^3+3aK^2L^2f^4-3K^2T^2f^2.                    (9)
```

The high-y filtration sees not the whole quartic at once, but the placement
of the three `K^2` terms relative to the two degrees `(m,n)`.

## 2. Proof of the strict tariff

Assume for contradiction that `m<n`; necessarily `n>=1`.  In `(9)`, the
term

```text
3aK^2L^2f^4                                               (10)
```

has y-degree `4n+4`.  It is uniquely highest:

- the other `K^2` terms have degrees at most `3n+4` and
  `2m+2n+4<=4n+2`;
- every term linear in `K` has degree at most `m+3n+2<=4n+1`;
- every `K`-free term has still smaller degree.

If `v=lc_y(f)`, equation `(3)` therefore gives

```text
lc_y(G)^2=3aL^2v^4.                                       (11)
```

But `L(-1)=-5`, so the right side has valuation

```text
ord_a(3aL^2v^4)=1+4ord_a(v),                              (12)
```

which is odd.  It cannot be a square in the UFD `k[x]`.  This proves `(5)`.

The obstruction is global: after passing to `k(x)`, the factor `a` is a
unit and its odd divisor valuation disappears.

## 3. The equality seam is an equianharmonic norm

Now suppose `m=n>=1`.  Exactly two terms of `(9)` have degree `4n+4`, and
their sum has leading coefficient

```text
3v^2(aL^2v^2-u^2).                                        (13)
```

This coefficient cannot vanish: `u^2=aL^2v^2` would again have odd
`a`-valuation.  Hence `deg_y(G)=2n+2`, and its leading coefficient `g`
satisfies

```text
g^2=3v^2(aL^2v^2-u^2).                                    (14)
```

Valuations in the UFD show `v|g`; putting `h=g/v` proves `(6)`.  Factoring
`h^2+3u^2` over `d^2=-3` gives `(7)`, and taking `a`-valuations gives `(8)`
because `a` and `L` are coprime.

The norm condition is genuinely nonempty.  With `v=1`, take

```text
h=(a+3L^2)/2,                 u=(3L^2-a)/(2d).             (15)
```

Then

```text
h-du=a,                       h+du=3L^2.                   (16)
```

This pays the entire leading norm, but no lower coefficient of `(3)` has
been asserted.  It is the canonical hostile control against mistaking the
necessary color law for a no-go.

## 4. Exact C3 coordinate match

Let `lambda^2-lambda+1=0` and put

```text
d=2lambda-1.                                               (17)
```

Then `d^2=-3`, and the two shifted quartic colors from THM-3895 become

```text
(3+d)/2=1+lambda,              (3-d)/2=2-lambda.           (18)
```

The finite-exact carrier packet currently being developed under the reserved
THM-3898 namespace independently singles out the address `2-lambda`.  That
packet is context, not a dependency or active theorem claim here.  The exact
identities `(7)` and `(18)` nevertheless show that the proposed carrier uses
the same two-color coordinate governing entry of the first nonzero sidecar.
What is preserved is the odd `a`-valuation assigned between the colors.
What is lost by the leading view is every even common factor and every lower
y-coefficient.

## 5. Scope and next exact test

This theorem proves only a necessary leading filtration.  It does not close
the equality seam, the regimes `m>n`, or the constant seam `m=n=0`; it does
not assert that the positive payment `(15)` extends to a square.  The
cheapest decisive continuation is to substitute the two color allocations
from `(7)` into the next two y-coefficients of `(3)` and test whether their
odd `a`-valuation can be continued without introducing an additional
denominator divisor.

## 6. Reproduction

```bash
python3 04-computation/jc2_nonzero_sidecar_y_degree_tariff_thm3899.py
python3 -O 04-computation/jc2_nonzero_sidecar_y_degree_tariff_thm3899.py
```

Both runs must byte-match
`05-knowledge/results/jc2_nonzero_sidecar_y_degree_tariff_thm3899.out`.
The companion has 28 active gates.  Its finite degree-deficit table and
symbolic identities are controls; the valuation argument above is the
all-degree proof.
