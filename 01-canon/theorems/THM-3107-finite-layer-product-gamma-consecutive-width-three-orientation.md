---
id: THM-3107
title: "Finite-layer product-Gamma consecutive width-three orientation"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  For every product of one through six positive-shape Gamma variables, the
  normalized support {0,1,2} has no nonzero polynomial whose first three
  moments vanish.  After the two universal THM-2824 divisibility numerators
  are expressed through q_n=product_j(n+theta_j), they factor as
  N1=-q0^4 q1 P_A and N2=-q0^5 q1^2 Q_A, where P_A and Q_A have strictly
  positive integer coefficients in the labelled shapes.  The exact
  coefficient banks contain up to 117,439 monomials each.  This is a finite
  six-layer theorem candidate, not an all-factor-count or arbitrary-support
  product-Gamma result.
source: root-product-gamma-width-three-2026-08-02
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3056-product-gamma-reciprocal-hypergeometric-and-hankel-reversal
  - THM-3100-product-gamma-response-monge-compactification-and-bad-prefix-spectrum
  - THM-3104-multiplicatively-infinitely-divisible-stieltjes-width-three-common-nullity
script: 04-computation/gmc_product_gamma_consecutive_width3_orientation_thm3107.py
output: 05-knowledge/results/gmc_product_gamma_consecutive_width3_orientation_thm3107.out
script_sha256: 5e22b40970fe2896417cd3be9d5474a6b52bc87c352f7adcdb5ede69de680ce3
output_sha256: 7dbe4ace3cead8f69ab957844aedeef2984b2e15ab4c73fa5f0036a77fa7066b
hash_basis: LF-normalized bytes
---

# THM-3107 -- finite-layer product-Gamma consecutive width-three orientation

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

THM-3100 leaves width three open because the factorial atomic orientation
already fails for the one-factor shape `theta=2`.  THM-3104 then proves that
no argument using only Stieltjes positivity, generalized Hankel total
positivity, multiplicative infinite divisibility, or alternating logarithmic
curvature can repair the gap.

The finite product-Gamma law carries one extra coordinate: its consecutive
moment ratio is a polynomial whose roots are all negative.  On the first
consecutive support this coordinate gives a direct coefficientwise
orientation.  The result below proves the first six factor counts exactly.

## 1. Product-Gamma first window

Fix

```text
1<=A<=6,                    theta_1,...,theta_A>0,      (1)

w_n=product_(j=1)^A (theta_j)_n,                       (2)
L(s^n)=w_n,                  f_n=s^n/w_n.              (3)
```

A harmless positive scale `c^n` in `(2)` is removed by replacing `s` with
`s/c`, so `(2)` loses no generality.  Probabilistically, `L` is integration
against the product of `A` independent positive-shape Gamma variables.  It
therefore has full support on `(0,infinity)`.

On the normalized support `{0,1,2}`, put

```text
U=f_1-f_0,                    V=f_2-f_1.                (4)
```

Every mean-zero polynomial on these slots is uniquely `alpha U+beta V`.
Use THM-2824's Gram and cubic tensors

```text
g11=L(U^2),       g12=L(UV),        g22=L(V^2),
t111=L(U^3),      t112=L(U^2V),
t122=L(UV^2),     t222=L(V^3).                         (5)
```

Full support makes the real binary Gram form positive definite.  Hence its
quadratic and cubic forms have a common complex projective zero exactly when
both division-free remainders vanish:

```text
I1=3t112 g11 g22-t222 g11^2-2t111 g12 g22,
I2=3t122 g11 g22-2t222 g12 g11-t111 g22^2.             (6)
```

Thus it is enough to orient either invariant strictly.

## 2. The negative-root polynomial coordinate

Let

```text
q_n=w_(n+1)/w_n=product_j(n+theta_j).                  (7)
```

For a general normalized moment vector `w_0=1,w_1,...,w_6`, exact symbolic
reduction of `(6)` gives

```text
I1=N1/(w1^5 w2^3),              I2=N2/(w1^5 w2^4),    (8)
```

where `N1,N2` have respectively `20,23` integral monomials.  Substitute

```text
w_n=product_(r=0)^(n-1) q_r.                           (9)
```

For every factor count in `(1)`, exact division gives

```text
N1=-q_0^4 q_1 P_A(theta_1,...,theta_A),
N2=-q_0^5 q_1^2 Q_A(theta_1,...,theta_A).              (10)
```

The two quotient polynomials are symmetric in the labelled Gamma layers and
obey

```text
P_A,Q_A in Z_(>0)[theta_1,...,theta_A].                (11)
```

Their exact coefficient censuses are

| `A` | terms of `P_A` | terms of `Q_A` | total degree | least coefficients `P/Q` |
|---:|---:|---:|---:|---:|
| 1 | 2 | 2 | 1 | 16 / 12 |
| 2 | 34 | 34 | 7 | 16 / 12 |
| 3 | 308 | 308 | 13 | 16 / 12 |
| 4 | 2,331 | 2,331 | 19 | 16 / 12 |
| 5 | 16,681 | 16,681 | 25 | 16 / 12 |
| 6 | 117,439 | 117,439 | 31 | 16 / 12 |

This is a finite symbolic identity, not a sample over shape values.  The
companion constructs the universal numerators first, substitutes the labelled
linear factors in `(7)`, performs exact multivariate division in
`Z[theta_1,...,theta_A]`, and inspects every coefficient.  It also checks
endpoint transposition symmetry in every bank.

Since every shape, `q_0`, and `q_1` is positive, `(8)--(11)` imply

```text
I1<0,                         I2<0.                    (12)
```

In particular the common-nullity condition `(6)` is impossible.

## 3. Detection consequence

For every law `(1)--(3)` and every nonzero

```text
H=x f_0+y f_1+z f_2,                                  (13)
```

one has

```text
L(H)!=0  or  L(H^2)!=0  or  L(H^3)!=0.                (14)
```

Indeed, after `L(H)=0`, write `H=alpha U+beta V`.  If the second and third
moments also vanished, their binary forms would share `[alpha:beta]`,
contradicting `(12)` and THM-2824's exact divisibility criterion.

At one layer the certificate is already transparent:

```text
P_1=16 theta+40,                 Q_1=12 theta+20.      (15)
```

For `theta=2`, `(12)` specializes to the exact THM-3100 control

```text
I1=-1/2,                         I2=-11/36.             (16)
```

Thus THM-3100's negative atomic `D` is only a failure of that particular
certificate; the full divisibility remainder remains strictly oriented.

## 4. Boundaries and next mechanism

1. **Factor count.**  Equation `(11)` is proved here only for `1<=A<=6`.
   The identical term support and stable coefficient floors strongly suggest
   a layer-transfer theorem, but no claim for `A>=7` is made.
2. **Support.**  The theorem treats exactly the normalized first window
   `{0,1,2}`.  A translated or nonconsecutive triple changes the powers of the
   tilting density in the second and third moments; it is not obtained by a
   silent shift of the shapes.
3. **Positive shapes.**  Positivity is load-bearing.  Algebraically,
   `P_1(-5/2)=0` and `Q_1(-5/3)=0`.  These are sign-boundary controls, not
   Stieltjes laws and not common-nullity examples.
4. **Zero layers.**  At `A=0`, `w_n=1` is point evaluation and the Gram form
   degenerates.
5. **Soft moment cones.**  THM-3104 is not contradicted: its ratio `q_n` is
   not a fixed-degree polynomial with negative real roots.  The present
   coefficient orientation retains precisely the finite Gamma-layer datum
   erased by Stieltjes/MID compactification.

The next structural target is an `A`-layer marked-permutation or network
factorization proving `(11)` uniformly in `A`, followed separately by a
tilted construction for translated and arbitrary-gap three-slot supports.

## 5. Exact evidence

Run

```text
python 04-computation/gmc_product_gamma_consecutive_width3_orientation_thm3107.py
python -O 04-computation/gmc_product_gamma_consecutive_width3_orientation_thm3107.py
```

Both transcripts must byte-match the stored output after LF normalization.
The computation uses explicit `RuntimeError` gates, so optimized mode retains
every truth-bearing check.

**QED (candidate pending independent hostile audit).**
