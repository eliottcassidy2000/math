---
id: THM-3107
title: "Finite-layer product-Gamma initial width-three orientation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  finite nonempty product of positive Gamma shapes, both binary
  quadratic--cubic divisibility invariants on normalized support {0,1,2}
  are strictly negative.  More generally every anchored arithmetic
  progression {0,d,2d} is good.  The proof is an all-layer coefficientwise
  orientation theorem: after reciprocal-shape substitution, both cleared
  invariants have zero coefficients below total degree five and, inside the
  per-layer multidegree box 0<=k_nu<=6, strictly positive coefficients at
  every total degree at least five.  A unique dominant
  response mode controls every sufficiently long layer histogram, while
  finite exact enumeration closes the two short banks.
source: root/multiscale-newton-flag/product-gamma-width3-2026-08-02
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-2853-gamma-adjacent-tensor-cycle-weighted-positivity
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3100-product-gamma-response-monge-compactification-and-bad-prefix-spectrum
script: 04-computation/gmc_product_gamma_initial_width_three_orientation_thm3107.py
output: 05-knowledge/results/gmc_product_gamma_initial_width_three_orientation_thm3107.out
script_sha256: 1ceec4ba945aa8795e7daa2d973e969b97a54060ad211392b904ab4817983b6c
output_sha256: 8ca250deb6e30e4635d720f7c7872d06cccfdce4ce54d946aefa8a813da0f03c
hash_basis: LF-normalized bytes
---

# THM-3107 -- finite-layer product-Gamma initial width-three orientation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3100 proves every finite product-Gamma family in widths one and two,
but its factorial width-three atomic determinant already has the wrong sign
for the single Gamma shape `2`.  The actual two divisibility invariants are
more rigid.  On the initial support `{0,1,2}`, each is coefficientwise
oriented in all reciprocal shape variables, for an arbitrary finite number
of Gamma layers.

The proof introduces a reusable finite signed product-kernel lemma.  A
single positive response mode dominates every negative mode coordinatewise.
It therefore wins on every sufficiently long layer word.  Only finitely many
short histograms remain, and their exact census reveals the sharp order-five
vanishing wall.  This is the same algebraic shape as a finite-state
relation/carry spectrum: one diagonal response alphabet, one dominant state,
and a bounded exceptional prefix bank.

## 1. A finite signed product-kernel lemma

Let `E` be finite.  Give each mode `e in E` a real coefficient `c_e` and a
nonnegative response vector

```text
q_e=(q_(e,1),...,q_(e,d)).                             (1)
```

Suppose one mode `e_*` has coefficient `c_*>0` and response

```text
Q=(Q_1,...,Q_d),                 Q_k>0,                (2)
```

with

```text
q_(e,k)<=Q_k                         for all e,k.       (3)
```

For every negative mode put

```text
rho_e=max_k q_(e,k)/Q_k<1.                            (4)
```

If `n=(n_1,...,n_d)` is a histogram and
`ell=sum_k n_k`, then the signed product-kernel coefficient

```text
C(n)=sum_e c_e product_k q_(e,k)^n_k                  (5)
```

satisfies

```text
C(n)
 >=product_k Q_k^n_k
   [c_*-sum_(c_e<0)|c_e|rho_e^ell].                   (6)
```

Indeed, discard every positive mode other than `e_*`, and bound each
negative product by `rho_e^ell product Q_k^n_k`.  Consequently

```text
sum_(c_e<0)|c_e|rho_e^ell<c_*       implies C(n)>0.    (7)
```

This lemma is deliberately finite and exact.  It needs no asymptotic
regularity.  Symbols with `Q_k=0` must be removed or treated separately.
An inactive symbol with every response equal to one contributes no length
and may simply be deleted.  Equality in `(6)` would require every used
negative coordinate to attain its maximal ratio and every discarded
positive mode to vanish; neither occurs beyond the explicit order-five wall
in the application below.

## 2. Product-Gamma moments and the exact obstruction

Fix

```text
A>=1,                    theta_1,...,theta_A>0,
w_n=product_(nu=1)^A(theta_nu)_n,
L_theta(s^n)=w_n,        f_n=s^n/w_n.                 (8)
```

On support `{0,1,2}` put

```text
U=f_1-f_0,                         V=f_2-f_1.          (9)
```

As in THM-2824, define

```text
g11=L(U^2),       g12=L(UV),       g22=L(V^2),
t111=L(U^3),      t112=L(U^2V),
t122=L(UV^2),     t222=L(V^3),                           (10)

I1=3t112 g11 g22-t222 g11^2-2t111 g12 g22,
I2=3t122 g11 g22-2t222 g12 g11-t111 g22^2.            (11)
```

The law of a product of independent Gamma variables is continuous and
charges every open subinterval of `(0,infinity)`.  Thus the Gram quadratic
is positive definite.  THM-2824's purely binary algebra applies unchanged:
its quadratic and cubic have a common complex projective root exactly when

```text
I1=I2=0.                                               (12)
```

We prove the stronger orientation

```text
I1<0,                         I2<0.                    (13)
```

## 3. Reciprocal-shape response coordinates

Set

```text
x_nu=1/theta_nu>0.                                    (14)
```

Multiplying every moment `w_n` by one geometric factor to the power `n`
does not change any normalized joint moment of the `f_n`.  Divide away
`(product theta_nu)^n`.  The successive moment responses become

```text
r_j=product_nu(1+j x_nu),             0<=j<=5,
r_0=1.                                                   (15)
```

Direct Laurent algebra in the five variables `r_1,...,r_5` gives

```text
F1=-r_1^2 I1,                    F2=-r_1^2 I2.         (16)
```

Both `F_j` are ordinary polynomials.  Write a response monomial as

```text
r^e=r_1^e1 ... r_5^e5,
Q_e(x)=product_(j=1)^5(1+jx)^e_j.                     (17)
```

The complete nonzero response banks `F_j=sum_e c_e r^e` are as follows.
Entries are `e:c_e`.

```text
F1 positive:
  40000:2, 32100:7, 32000:6, 23100:2, 22000:4,
  21111:1, 21110:6, 21100:5, 12200:3, 01111:1;

F1 negative:
  41000:-1, 33000:-4, 31110:-3, 31100:-3, 31000:-7,
  22200:-3, 22100:-9, 12100:-2, 11111:-2, 11110:-3.
                                                               (18)

F2 positive:
  41000:1, 33000:4, 31110:3, 31000:4, 22200:6,
  22100:9, 13200:1, 12111:2, 12100:4, 11111:2,
  02210:3;

F2 negative:
  40000:-1, 32100:-7, 32000:-4, 23100:-4, 22000:-4,
  21111:-2, 21110:-3, 21100:-2, 12210:-3, 12200:-6,
  02200:-1, 02111:-2.                                      (19)
```

The companion derives `(18)--(19)` from `(10)--(11)`; it does not import
the displayed banks.

Since `(15)` is a product over layers,

```text
r^e=product_nu Q_e(x_nu).                              (20)
```

Put

```text
q_(e,k)=[x^k]Q_e(x),                    0<=k<=6.       (21)
```

For a monomial `product_nu x_nu^k_nu`, its coefficient in `F_j` is exactly

```text
C_j(k_1,...,k_A)=sum_e c_e product_nu q_(e,k_nu).      (22)
```

The symbol `k=0` is inactive because `q_(e,0)=1` for every mode.  Therefore
`(22)` depends only on the histogram of the active symbols `1,...,6`.

## 4. Unique dominant modes and exact tails

For `F1`, the unique positive mode dominating all other response vectors is

```text
e_*=(2,1,1,1,1),        c_*=1,
Q=(1,16,100,310,499,394,120).                         (23)
```

Here and below the initial `1` is the inactive degree-zero coordinate.  For
every negative mode the maximum ratio `q_(e,k)/Q_k`, `1<=k<=6`, occurs at
`k=1`.  After collecting equal ratios, the right side of `(7)` is controlled
by

```text
B1(ell)=
 (3/8)^ell+13(9/16)^ell+6(3/4)^ell+5(1/2)^ell
 +7(5/16)^ell+2(15/16)^ell+3(5/8)^ell.                (24)
```

The first passing length for this dominant-mode certificate is `14`, and

```text
B1(14)=33354436230658959/36028797018963968<1=c_*.     (25)
```

For `F2`, the unique dominant positive mode is

```text
e_*=(1,2,1,1,1),        c_*=2,
Q=(1,17,115,395,724,668,240).                         (26)
```

Again every negative mode has its largest coordinate ratio at `k=1`.  The
collected tail is

```text
B2(ell)=
 (4/17)^ell+8(10/17)^ell+6(7/17)^ell
 +13(11/17)^ell+4(6/17)^ell
 +4(16/17)^ell+3(15/17)^ell,                          (27)
```

whose first passing length is `16`:

```text
B2(16)=94169758536521331158/48661191875666868481<2.   (28)
```

The product-kernel lemma now proves every coefficient with at least `14`
active layers for `F1`, and at least `16` for `F2`, is strictly positive.
Every displayed contraction ratio lies in `[0,1)`, so both tail bounds are
decreasing; passage at the displayed threshold propagates to every larger
active-layer count.

## 5. The finite bank and sharp order-five wall

It remains to inspect all six-symbol histograms of lengths at most `13` for
`F1` and at most `15` for `F2`.  These are exactly

```text
sum_(ell=0)^13 binom(ell+5,5)=binom(19,6)=27132,
sum_(ell=0)^15 binom(ell+5,5)=binom(21,6)=54264       (29)
```

histograms.  Substitution into `(22)` uses integers only.  The exhaustive
result is

```text
C_j(n_1,...,n_6)=0
 iff
sum_(k=1)^6 k n_k<5,

C_j(n_1,...,n_6)>0
 iff
sum_(k=1)^6 k n_k>=5.                                (30)
```

There are exactly `12` zero histograms in each bank.  The smallest positive
coefficients are

```text
[x_nu^5]F1=16,                  [x_nu^5]F2=12.        (31)
```

Equations `(24)--(28)` cover every longer histogram, so `(30)` holds for an
arbitrary finite number of layers, not merely for the enumerated lengths.
Equivalently, `F1` and `F2` have no terms below total degree five and, in the
per-layer multidegree box `0<=k_nu<=6`, every coefficient of total degree at
least five is strictly positive.  Monomials outside that box are absent.

Because `A>=1` and every `x_nu>0`, `(31)` gives

```text
F1>0,                           F2>0.                 (32)
```

Now `(16)` proves `(13)`, and `(12)` proves that `{0,1,2}` is good for every
finite product-Gamma family.

## 6. Anchored arithmetic progressions

The conclusion is stable under dilation from the origin.  Fix `d>=1` and
put `t=s^d`.  The restricted moment sequence is

```text
W_n=w_(dn)=product_nu(theta_nu)_(dn).                 (33)
```

Gauss multiplication gives the exact identity

```text
(theta)_(dn)
 =d^(dn) product_(r=0)^(d-1)((theta+r)/d)_n.          (34)
```

Thus, up to one geometric factor to the power `n`, `(W_n)` is the moment
sequence of `Ad` positive Gamma shapes

```text
(theta_nu+r)/d,                  1<=nu<=A, 0<=r<d.    (35)
```

The geometric factor cancels from normalized moments.  Applying the initial
cell theorem in the variable `t` proves

```text
{0,d,2d} is good                    for every d>=1.    (36)
```

This dilation statement is exact; it is not an asymptotic remote theorem.

## 7. Sharp boundaries and hostiles

### Empty and deleted layers

For `A=0`, all normalized moments equal evaluation at one and both invariants
vanish.  This is exactly the degree-below-five boundary of `(30)`.  Allowing
`x_nu=0` merely deletes that layer.  If at least one reciprocal shape remains
positive, `(31)` still makes both cleared invariants positive.  Physical
finite shapes have every `x_nu>0`.

### The old atomic route really fails

For the single Gamma shape `theta=2`, exact arithmetic gives

```text
I1=-1/2,                   I2=-11/36,
2t222 g12-3t122 g22=-1/12.                              (37)
```

Thus the desired invariants have the strict sign proved here even though
THM-2824's sufficient factorial atomic determinant has the wrong sign.

### Positive polynomial coefficients are insufficient

The factorization into positive linear responses is load-bearing.  Consider
the purely algebraic response

```text
P(n)=n^3+n^2+6n+1000,
r_j=P(j)/P(0).                                         (38)
```

All coefficients of `P` are positive, but its discriminant is

```text
-26896828<0,                                           (39)
```

so it is not a product of positive-root linear factors.  Substitution in the
same exact response banks gives

```text
F1=-1229376/3814697265625<0,
F2= 1539056/3814697265625>0.                           (40)
```

This is a response-cone hostile, not a Stieltjes moment or GMC
counterexample.  It shows that arbitrary positive coefficients cannot
replace the product-Gamma layer structure.

### Translation does not commute with the proof

No translated arithmetic progression is claimed.  If
`H=f_a K(s^d)`, then

```text
L_theta(H^r)
 =(w_(ra)/w_a^r) L_(theta+ra)(K^r).                   (41)
```

The quadratic uses the tilt `theta+2a`, the cubic uses `theta+3a`, while the
coefficients of `K` come from the normalization at `theta+a`.  There is no
single product-Gamma functional to which Sections 3--5 apply.  Gauss
multiplication removes dilation but not this `2a` versus `3a` tilt mismatch.

## 8. Exact companion and scope

The dependency-free companion:

1. derives the `20`- and `23`-mode banks from the moment tensors;
2. finds the two unique dominant modes and proves every negative ratio is
   maximized at degree one;
3. verifies the exact first tail thresholds `14` and `16`;
4. exhausts all `27132+54264` short histograms and the order-five boundary;
5. checks the shape-two atomic hostile and the non-real-root response
   hostile;
6. independently compares direct tensors with the cleared response banks on
   `9` rational shape tuples with one, two, and three layers; and
7. verifies `126` rational Gauss multiplication cells.

An orthogonal finite-layer companion from the reservation lineage is
preserved unchanged at

```text
04-computation/gmc_product_gamma_consecutive_width3_orientation_thm3107.py
05-knowledge/results/gmc_product_gamma_consecutive_width3_orientation_thm3107.out.
```

It works directly in the labelled shape variables and independently derives

```text
N1=-q_0^4 q_1 P_A,                 N2=-q_0^5 q_1^2 Q_A
```

with every coefficient of `P_A,Q_A` positive for `1<=A<=6`.  Its exact term
censuses are

```text
A:       1       2       3       4        5         6
terms:   2      34     308    2331    16681    117439,
```

and its least coefficients are always `16/12`.  Normal and optimized runs
both reproduce its stored fourteen-line transcript.  Its LF hashes are

```text
script  5e22b40970fe2896417cd3be9d5474a6b52bc87c352f7adcdb5ede69de680ce3
output  7dbe4ace3cead8f69ab957844aedeef2984b2e15ab4c73fa5f0036a77fa7066b.
```

That sidecar certifies only `A<=6`; the signed product-kernel argument in
Sections 4--5 is the separate step supplying the all-finite-`A` quantifier.
Accordingly its stored `A7plus_open` boundary describes the sidecar method,
not the scope of this theorem.

Run

```text
python 04-computation/gmc_product_gamma_initial_width_three_orientation_thm3107.py
python -O 04-computation/gmc_product_gamma_initial_width_three_orientation_thm3107.py
```

Both executions must byte-match the stored ten-line transcript after LF
normalization.

Two independent hostile audits rebuilt the Laurent tensors and both signed
response banks rather than trusting the displayed formulas.  They separately
checked the clearing sign in `(16)`, the unique coordinatewise dominators,
all negative contraction ratios, monotone propagation beyond the exact
`14/16` thresholds, the complete short-histogram wall, the per-layer
multidegree qualification in `(30)`, and the Gauss dilation and translated
tilt boundaries.  Normal and optimized executions of both companions
byte-match their stored transcripts, and all four LF hashes agree with the
declared values.

This theorem proves the initial cell and every anchored three-term arithmetic
progression.  It does not prove arbitrary translated or non-arithmetic
three-slot product-Gamma goodness, the conditional `O(X^(t-4))` spectrum of
THM-3100, SFC in every width, NC2, arbitrary-radial GMC(2), LRC(14), JC(2),
or DC(2).

**QED.**
