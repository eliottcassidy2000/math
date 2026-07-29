---
id: THM-2815
title: "Optimal finite Laguerre carrier and radial-selector access boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The monic
  Laguerre quotient of dimension D
  represents L(s^n)=n! exactly through degree 2D-1 and is the unique optimal
  s-generated multiplicative carrier at that horizon.  Its first alias is
  ell_D^2, with L(ell_D^2)=(D!)^2.  The inverse factorial-Hankel matrix gives
  exact finite-height coefficient selectors, but scalar moment nullity does
  not provide the required multiplier observations.  This is not arbitrary
  Wick-channel separation, HYP-8765, or a new GMC proof.
source: root/optimal-finite-laguerre-carrier-2026-07-28
audit: root/audit-2809-2026-07-28 (independent finite-difference,
  quotient-horizon, optimality, uniqueness, selector, hostile,
  normal/-O/stored/hash and docs audit: ACCEPT)
depends_on: []
related:
  - THM-1620-the-pochhammer-bridge-toral-legendre-radial-hermite
  - THM-1790-the-emp-floor-detection-depth-at-least-degree-plus-one
  - THM-2638-radial-height-graded-wick-decoder-and-laplace-forgetting-boundary
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-2812-consecutive-three-slot-factorial-moment-six-detection
script: 04-computation/gmc_optimal_finite_laguerre_carrier_thm2815.py
output: 05-knowledge/results/gmc_optimal_finite_laguerre_carrier_thm2815.out
script_sha256: 7a6d9670b16cfd6db651f978438fdc0606d368adcce6c5a794b541fee6f0c262
output_sha256: 903c60e7abaa7108b0729577ba65071350697f547988a2949397b67cec82514d
addendum_status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over R the
  optimal Laguerre quotient splits into D simple positive nodes with the
  strictly positive Gauss--Laguerre weights
  1/(x_i L_D'(x_i)^2).  This is the unique D-node atomic carrier through
  degree 2D-1, already needs D nodes through degree 2D-2, and first misses
  degree 2D by (D!)^2.  The positive splitting is cutoff-dependent and
  does not turn scalar factorial nullity into nodewise nullity or supply
  the missing multiplier observations.
addendum_audit: >
  root/audit-2809-2026-07-28 (independent positive-root, sign-change,
  Christoffel--Darboux normalization, complex-node uniqueness, cutoff,
  exact quotient-trace, normal/-O/stored/hash, and scope audit: ACCEPT)
addendum_script: 04-computation/gmc_positive_gauss_laguerre_splitting_addendum_thm2815.py
addendum_output: 05-knowledge/results/gmc_positive_gauss_laguerre_splitting_addendum_thm2815.out
addendum_script_sha256: b770897ec5836aa37076bc0753e0e7980bb9b1408b3f7c19d356b722bc4c3fd2
addendum_output_sha256: 5a41f6cfef9bcf39a5117d50c746e68a3378b44402b9e521ea5100e54ab614d9
hash_basis: LF-normalized bytes
---

# THM-2815 -- the unique optimal finite factorial carrier is Laguerre

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2810 proves that no fixed finite-dimensional carrier can represent the
whole factorial tower in characteristic zero.  This theorem identifies the
sharp positive boundary it leaves open: at every finite height there is a
unique optimal carrier, and it is the monic Laguerre quotient.

The same factorization gives exact dual coefficient selectors.  They expose,
rather than erase, the remaining Gaussian obstruction: the selectors require
new multiplier observations which the scalar moment tower does not supply.

## 1. Monic Laguerre relation

Let

```text
L:Q[s] -> Q,                         L(s^n)=n!,           (1)
```

and, for `D>=1`, define

```text
ell_D(s)=(-1)^D D! L_D(s)
 =sum_(k=0)^D (-1)^(D+k) D! binom(D,k) s^k/k!.          (2)
```

The displayed formula is the definition needed here; no analytic
orthogonal-polynomial theorem is imported.  Its leading coefficient is one.

For `r>=0`,

```text
L(s^r ell_D)
 =D! sum_(k=0)^D (-1)^(D+k)binom(D,k)(k+r)!/k!.         (3)
```

The function

```text
k |-> (k+r)!/k!=(k+1)(k+2)...(k+r)                     (4)
```

is a polynomial of degree `r`.  Equation `(3)` is `D!` times its `D`th
forward difference.  Hence

```text
L(s^r ell_D)=0,                         0<=r<D,          (5)
L(s^D ell_D)=(D!)^2.                                      (6)
```

## 2. Exact finite-height quotient

Divide every polynomial `f` of degree at most `2D-1` uniquely as

```text
f=q ell_D+r,                    deg q<D,  deg r<D.       (7)
```

By `(5)`,

```text
L(f)=L(r).                                             (8)
```

Therefore

```text
A_D=Q[s]/(ell_D)                                      (9)
```

with readout equal to `L` on the unique representative of degree below `D`
is an exact `D`-dimensional multiplicative carrier for every factorial
moment through degree `2D-1`.

The rational quotient `(9)` and identity `(8)` are the proved mechanism used
here; no node-splitting, positivity, or simplicity claim is needed.

## 3. The first failure is exact

Because `ell_D` is monic, `ell_D-s^D` has degree below `D`.  Thus `(5)` gives

```text
L(ell_D^2)
 =L(s^D ell_D)+L((ell_D-s^D)ell_D)
 =(D!)^2.                                               (10)
```

But `ell_D^2=0` in the quotient `(9)`.  Hence degree `2D` is an explicit
first alias:

```text
quotient readout of ell_D^2=0,
L(ell_D^2)=(D!)^2!=0.                                  (11)
```

The horizon `2D-1` is sharp.

## 4. Optimal dimension and unique relation

The factorial Hankel block

```text
H_(D-1)=((i+j)!)_(0<=i,j<D)                            (12)
```

has determinant

```text
det H_(D-1)=product_(j=0)^(D-1)(j!)^2!=0.              (13)
```

Any time-homogeneous linear carrier exact through degree `2D-2` factors
`H_(D-1)` through its state space.  Its dimension is therefore at least
`D`.  The carrier `(9)` attains this lower bound and is exact one degree
farther.

There is also uniqueness in the multiplicative category.  Let an
`s`-generated algebra of dimension `D` represent `L` through degree
`2D-1`.  The image of `s` has a nonzero monic relation `q` of degree at most
`D`.  Exactness gives

```text
L(s^r q)=0,                           0<=r<D.            (14)
```

If `deg q<D`, the nondegenerate block `(12)` forces `q=0`, a contradiction.
Thus `deg q=D`.  The same nondegeneracy makes the monic solution of `(14)`
unique; `(5)` identifies it as

```text
q=ell_D.                                                  (15)
```

So `(9)` is the unique optimal `s`-generated multiplicative relation at this
horizon.

## 5. Exact inverse-Hankel selectors

For `d>=0`, put

```text
H_d=((i+j)!)_(0<=i,j<=d).                              (16)
```

THM-2810's self-contained Pascal factorization gives

```text
(H_d^(-1))_(ij)
 =(-1)^(i+j)/(i!j!)
   sum_(k=max(i,j))^d binom(k,i)binom(k,j).             (17)
```

Indeed `H_d=D P P^T D`, and the inverse Pascal entries are
`(-1)^(i-j)binom(i,j)`.

Define

```text
phi_j^(d)(s)=sum_(i=0)^d (H_d^(-1))_(ij)s^i.            (18)
```

Matrix inversion gives the exact duality

```text
L(s^r phi_j^(d))=delta_(rj),             0<=r,j<=d.     (19)
```

The selector norm is

```text
L((phi_j^(d))^2)
 =(H_d^(-1))_(jj)
 =(j!)^(-2) sum_(k=j)^d binom(k,j)^2.                   (20)
```

Consequently every `G` of degree at most `d` obeys

```text
[s^j]G=L(G phi_j^(d)).                                  (21)
```

This is a complete finite-height coefficient decoder when the multiplier
observations on the right of `(21)` are available.

## 6. The multiplier-access boundary

Gaussian nullity supplies scalar equations `L(G_m)=0`; it does not supply
the family `L(G_m phi_j)`.  MISTAKE-211 is the sharp live witness.  With

```text
G_4=4ab^3s^6+4a^3cs^18,                                (22)
a=b=1,                         c=-6!/18!,               (23)
```

both coefficients in `(22)` are nonzero but

```text
L(G_4)=0.                                               (24)
```

The degree-eighteen selectors correctly recover

```text
L(G_4 phi_6^(18))=4,
L(G_4 phi_18^(18))=-1/2223046425600.                   (25)
```

Equations `(25)` do not follow from `(24)`.  Similarly, the quotient `(9)`
preserves the scalar factorial sum exactly; it does not make zero readout
imply a zero quotient class.

The connection contract is therefore:

| item | exact content |
|---|---|
| source | one finite radial shell `G_m` |
| optimal carrier | `Q[s]/(ell_D)` through height `2D-1` |
| exact dual | inverse-Hankel selectors `(18)` |
| preserved | scalar factorial readout |
| missing sidecar | lawful observations `L(G_m phi_j)` or algebra-valued nullity |
| sharp hostile | MISTAKE-211 two-height cancellation `(22)--(25)` |

The theorem sharpens THM-2810 but does not separate arbitrary Wick channels,
prove HYP-8765, or replace THM-2022.

## 7. Exact companion

Run

```bash
python 04-computation/gmc_optimal_finite_laguerre_carrier_thm2815.py
python -O 04-computation/gmc_optimal_finite_laguerre_carrier_thm2815.py
```

Both executions byte-match

```text
05-knowledge/results/gmc_optimal_finite_laguerre_carrier_thm2815.out.
```

The companion verifies the carrier and first failure for `D=1,...,12`,
all inverse-Hankel entries and selector norms through `d=10`, and the exact
degree-eighteen MISTAKE-211 recovery.  It has explicit exception gates, no
truth-bearing Python assertions, no floating point, and no scratch
dependency.  The universal statements are proved above rather than inferred
from these finite controls.

**QED.**

## Independently audited addendum: the optimal quotient splits positively

The rational quotient in `(9)` has a sharper real form.  Write

```text
ell_D(x)=(-1)^D D! L_D(x),                              (26)
```

with `L_D=L_D^(0)` in the standard normalization.

> **Positive-splitting theorem.**  The polynomial `ell_D` has `D` simple
> roots
>
> ```text
> 0<x_1<...<x_D.
> ```
>
> The evaluation map splits the optimal carrier as
>
> ```text
> A_D tensor_Q R  ->  product_(i=1)^D R,
> f                |-> (f(x_1),...,f(x_D)).             (27)
> ```
>
> Its readout through degree `2D-1` is the strictly positive atomic
> functional
>
> ```text
> L(f)=sum_(i=1)^D w_i f(x_i),
> w_i=1/[x_i L_D'(x_i)^2]
>    =x_i/[(D+1)^2 L_(D+1)(x_i)^2]>0.                  (28)
> ```
>
> This is the unique `D`-node real or complex atomic carrier through
> degree `2D-1`, up to ordering.  Fewer than `D` nodes are impossible
> already through degree `2D-2`.  At degree `2D` its moment is
>
> ```text
> sum_i w_i x_i^(2D)=(2D)!-(D!)^2,                     (29)
> ```
>
> so the first missing factorial moment is exactly `(D!)^2`.

### Positive roots

Equation `(5)` says that `ell_D` is orthogonal to every polynomial of
degree below `D` for

```text
<f,g>=integral_0^infinity f(x)g(x)e^(-x) dx=L(fg).     (30)
```

Suppose `ell_D` had fewer than `D` distinct positive odd-multiplicity
roots, and let `q` be their product.  Then `deg q<D`, while

```text
ell_D q=q^2(ell_D/q)
```

has one fixed sign on `(0,infinity)` and is not zero.  Its integral in
`(30)` cannot vanish, contradicting `L(q ell_D)=0`.  Thus all `D` roots
are distinct and positive.

### Christoffel--Darboux weights

The standard polynomials `L_0,...,L_(D-1)` are orthonormal for `(30)`.
Put

```text
K_(D-1)(x,y)=sum_(k=0)^(D-1) L_k(x)L_k(y).             (31)
```

Christoffel--Darboux at a root `x_i` gives

```text
K_(D-1)(x_i,x_i)=x_i L_D'(x_i)^2.                     (32)
```

The cardinal polynomial

```text
lambda_i(x)=K_(D-1)(x,x_i)/K_(D-1)(x_i,x_i)
```

has degree below `D` and values `lambda_i(x_j)=delta_(ij)`.  Only its
`L_0=1` component survives integration, so

```text
L(lambda_i)=1/K_(D-1)(x_i,x_i)=w_i>0.                 (33)
```

For `deg f<=2D-1`, divide `f=q ell_D+r` with `deg q,deg r<D`.
Orthogonality kills the first term and cardinal interpolation integrates
the second, proving `(28)`.  The root identity

```text
(D+1)L_(D+1)(x_i)=x_i L_D'(x_i)
```

gives the second weight formula.

### Minimality, uniqueness, and first failure

The Hankel-rank argument in Section 4 already proves the `D`-node lower
bound through degree `2D-2`.  For uniqueness, delete zero weights and
combine coincident nodes in a putative `D`-node carrier.  Minimality leaves
exactly `D` distinct effective nodes.  Their monic node polynomial `q`
satisfies

```text
L(x^r q)=0,                         0<=r<D.             (34)
```

The coefficient vector of `q-ell_D` lies in the kernel of the invertible
factorial Hankel block `H_(D-1)`, even over `C`; hence `q=ell_D`.
Vandermonde inversion then makes the weights unique.  Finally,
`ell_D(x_i)=0` and `(10)` give

```text
sum_i w_i ell_D(x_i)^2=0,
L(ell_D^2)=(D!)^2,                                      (35)
```

which is equivalent to `(29)`.

### Exact quotient-trace certificate

No algebraic root approximation is needed.  In

```text
A_D=Q[x]/(ell_D)
```

define

```text
omega_D=(D!)^2 [x ell_D'(x)^2]^(-1) mod ell_D.          (36)
```

The inverse exists because the roots are nonzero and simple, and

```text
Tr_(A_D/Q)(omega_D f)=sum_i w_i f(x_i).                 (37)
```

The addendum companion verifies `(26)--(37)` exactly for `1<=D<=12`:
Sturm positive-root counts; squarefreeness; orthogonality and norm;
Christoffel--Darboux and both weight normalizations modulo `ell_D`;
all quotient-trace moments through `2D-1`; the first miss at `2D`;
Hankel determinants; recurrence; and coprimality of successive `ell_D`.

### Sharp scalar-nullity and tower boundary

For `D=2`,

```text
ell_2=x^2-4x+2,                 omega_2=1-x/4.           (38)
```

The polynomial `f=x-1` has `L(f)=0` but is nonzero at both roots
`2-sqrt(2)` and `2+sqrt(2)`.  Thus positivity of the atomic weights does
not turn one scalar observation into nodewise zero, and it does not supply
the multipliers in `(21)`.  The algebraic nodes depend on the cutoff and
successive `ell_D` are coprime, so these splittings do not form a bounded
compatible tower.  This addendum sharpens the finite carrier boundary; it
does not prove a new NC2 or GMC statement.

Reproduce the addendum by

```text
python 04-computation/gmc_positive_gauss_laguerre_splitting_addendum_thm2815.py
python -O 04-computation/gmc_positive_gauss_laguerre_splitting_addendum_thm2815.py
```

Both modes byte-match

```text
05-knowledge/results/gmc_positive_gauss_laguerre_splitting_addendum_thm2815.out.
```

**QED (positive-splitting addendum).**
