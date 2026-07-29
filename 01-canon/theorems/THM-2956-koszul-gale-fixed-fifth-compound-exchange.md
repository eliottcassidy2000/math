---
id: THM-2956
title: "Koszul--Gale fixed fifth-compound exchange"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For
  universal ternary forms of degrees (2,3,4), the six quadratic
  Koszul syzygies give the exact fixed-cofactor identity
  q_200 P_fixed+c_300 P_opt=0 between a 20Q+10C+5F rank-35
  Macaulay minor and its 21Q+9C+5F exchange.  On every factorial
  specialization q_200 is positive, so nonvanishing of the fixed
  cofactor forces nonvanishing of the exchanged cofactor.  Its
  cleared degree invoice is 54M-35 rather than 55M-35.  This is a
  universal Pluecker mutation and certificate-transfer theorem, not
  by itself a new support atlas, positivity theorem, resultant
  factorization, SFC(4), or GMC(2) proof.
source: codex-gmc-koszul-gale-exchange-2026-07-29
audit: >
  An independent hostile audit rederived both complementary Gale
  minors, the opposite Pluecker sign, the generic-locus polynomial
  extension, and the factorial consequence.  Normal, optimized, and
  stored exact transcripts match after LF normalization with empty
  stderr.
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
related:
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
  - THM-2949-fixed-rank-thirty-five-cofactor-newton-atlas
  - THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary
script: 04-computation/gmc_koszul_gale_fixed_fifth_compound_exchange_thm2956.py
output: 05-knowledge/results/gmc_koszul_gale_fixed_fifth_compound_exchange_thm2956.out
script_sha256: cfd1c996f9570bdd646462788a9ce8e649efedcc11dd27290f629cd454c39b9c
output_sha256: 6dac6b8701474cc8a52d00c5219cce8f40bdad2f2e383d4e3ce72f0c1bf8de47
constructor_dependency_sha256: d2f8afeba7dd6c7950405a4845d7bf112b6c9872dd8161146446be8bbdaae0ba
hash_basis: LF-normalized bytes
---

# THM-2956 -- Koszul--Gale fixed fifth-compound exchange

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Universal statement

Let `Q,C,F` be ternary forms of degrees `2,3,4`.  Use the inherited
degree-seven Macaulay row order

```text
Q_0,...,Q_20, C_0,...,C_14, F_0,...,F_9.             (1)
```

Here `Q_i` is one degree-five monomial times `Q`, `C_j` is one
degree-four monomial times `C`, and `F_k` is one degree-three
monomial times `F`.

Choose any common five quartic rows, in the same order, and any common
`35` of the `36` target columns.  Define

```text
P_fixed
 =det(Q_0,...,Q_19,C_0,...,C_8,C_14, five F rows),

P_opt
 =det(Q_0,...,Q_20,C_0,...,C_8,      five F rows),    (2)
```

using those same target columns.  Put

```text
q_200=[x0^2]Q,                     c_300=[x0^3]C.     (3)
```

With the inherited row order,

```text
q_200 P_fixed+c_300 P_opt=0.                         (4)
```

Equation `(4)` is a polynomial identity in the universal form
coefficients.  It is not restricted to a factorial family or to the
regular locus used in its proof.

## 2. Six-syzygy Gale proof

The `21+15=36` quadratic and cubic rows have the six degree-two
Koszul relations

```text
m(C,-Q),                         m in Mon(S_2).       (5)
```

On the generic regular-sequence locus their row span has rank `30`,
and the coefficient matrix of `(5)` is a `6`-by-`36` Gale matrix.

The Gale columns indexed by the six Macaulay rows omitted from the
quadratic/cubic part of `P_opt` are

```text
C_9,...,C_14.                                          (6)
```

In the inherited monomial order their complementary determinant is
triangular and equals

```text
q_200^6.                                               (7)
```

The Gale columns indexed by the six rows omitted from the
quadratic/cubic part of `P_fixed` are

```text
Q_20,C_9,...,C_13.                                    (8)
```

Expansion in the pure `Q_20` column leaves the same triangular
multiplication matrix of size five, so this complementary determinant
is

```text
c_300 q_200^5.                                        (9)
```

Using the global indices `0,...,35`, the two selected index sums are

```text
sum(fixed indices)=450,
sum(opt indices)=435.                                (10)
```

They have opposite parity.  Gale duality therefore gives the opposite
Pluecker sign for `(7)` and `(9)`.  Wedge the two proportional
thirty-row exterior products with the same five `F` rows and take the
same target-column coordinate.  After cancelling `q_200^5`, the result
is exactly `(4)`.

The proof first holds on the nonempty generic regular-sequence,
rank-thirty locus.  Both sides of `(4)` are polynomials in the
coefficients of `Q,C,F`, so density extends the identity universally.

## 3. Factorial certificate transfer

On a translated four-slot factorial support of exact width `M`,
PROVED THM-2925 gives

```text
q_200(n)>0                       for every real n>=0. (11)
```

This also follows from positive definiteness of the unscaled Gram
quadratic and positivity of its clearing scalar.

At a physical depth, `(4)` now implies

```text
P_fixed(n)!=0
  ==> c_300(n) P_opt(n)=-q_200(n)P_fixed(n)!=0
  ==> c_300(n)!=0 and P_opt(n)!=0.                   (12)
```

Thus the exchange creates no new zero wherever the fixed cofactor is a
certificate.  No separate sign hypothesis on `c_300` is needed.

The two THM-2925 degree invoices are

```text
deg P_fixed
 <=20(M-1)+10(2M-1)+5(3M-1)=55M-35,

deg P_opt
 <=21(M-1)+ 9(2M-1)+5(3M-1)=54M-35.                 (13)
```

The universal identity therefore transfers any proved fixed-cofactor
atlas to a fixed chart whose formal degree is lower by exactly `M`.
The atlas itself remains a separate dependency: this theorem does not
turn a candidate or reserved finite bank into proved canon.

## 4. What the exchange preserves and loses

```text
source:
  one fixed 20Q+10C+5F rank-35 Pluecker coordinate;

map:
  exchange C_14 for Q_20 through the six Koszul relations;

preserved:
  nonvanishing at every point with q_200!=0;
  the five quartic rows and the target-column coordinate;

improved:
  the cleared factorial degree invoice, by M;

not preserved or supplied:
  the numerical determinant value, a canonical projective chart,
  resultant divisibility, real-ray positivity, or an initial
  nonvanishing atlas.                                  (14)
```

In particular, `(4)` is not the false claim that a `35`-minor is
universally divisible by the `(2,3,4)` resultant.  The exchange is
Pluecker-linear and chart-relative.

## 5. Exact companion

The companion symbolically constructs the full `6`-by-`36` Koszul
Gale matrix.  It verifies `(7)`, `(9)`, and `(10)` exactly.  It then
reconstructs four factorial controls and checks `(4)` by direct
`35`-by-`35` determinants at three depths each:

```text
(0,1,2,3),
(0,1,2,11),
(0,1,6,11),
(0,3,7,11).                                           (15)
```

These twelve determinant evaluations audit the inherited row order
and displayed sign; the exterior-algebra argument proves the universal
identity.

Run

```text
python 04-computation/gmc_koszul_gale_fixed_fifth_compound_exchange_thm2956.py
python -O 04-computation/gmc_koszul_gale_fixed_fifth_compound_exchange_thm2956.py
```

Both modes byte-match the stored transcript after LF normalization,
with empty stderr.

## 6. Scope

Proved:

```text
* the universal identity (4);
* the physical nonvanishing transfer (12);
* the degree reduction (13);
* the exact symbolic and factorial controls.          (16)
```

Not proved:

```text
* a new support-width closure without a proved input atlas;
* a nonzero full determinant or resultant factorization;
* positivity between integer depths;
* arbitrary-width SFC(4), SFC(5), NC2, or GMC(2).     (17)
```

**QED.**
