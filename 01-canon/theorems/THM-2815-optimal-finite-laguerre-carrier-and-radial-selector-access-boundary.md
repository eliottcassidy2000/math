---
id: THM-2815
title: "Optimal finite Laguerre carrier and radial-selector access boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The monic
  Laguerre quotient of dimension D represents L(s^n)=n! exactly through
  degree 2D-1 and is, up to the unique pointed readout-preserving
  isomorphism, the optimal D-dimensional s-generated multiplicative carrier
  at that horizon.  The first readout-alias degree is 2D, canonically
  witnessed by ell_D^2 with L(ell_D^2)=(D!)^2.  At arbitrary horizon H the
  minimum dimension is floor(H/2)+1; even horizons have the affine relation
  family ell_D+c ell_(D-1), while odd horizons select ell_D uniquely.  The
  inverse factorial-Hankel matrix gives exact finite-height coefficient
  selectors, but scalar moment nullity does not provide the required
  multiplier observations.  This is not arbitrary Wick-channel separation,
  HYP-8765, or a new GMC proof.
source: root/optimal-finite-laguerre-carrier-2026-07-28
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
independent_script: 04-computation/gmc_optimal_finite_laguerre_carrier_independent_audit_thm2815.py
independent_output: 05-knowledge/results/gmc_optimal_finite_laguerre_carrier_independent_audit_thm2815.out
independent_script_sha256: dfd44fb830ccaccaa78b027684ea8b4e1707ba6cff4cc4ee4b38b4c2ed00ed6f
independent_output_sha256: e8e189e152068661e40dd02e1c989d33d974be774a875886f623f60010a42f98
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

Fix the category before stating optimality.  A pointed `s`-generated
multiplicative carrier of horizon `H` is a triple

```text
(A,u,lambda),                                           (7)
```

where `A` is a finite-dimensional unital `Q`-algebra, `A=Q[u]`, and
`lambda:A->Q` is linear, such that

```text
lambda(u^n)=n!,                         0<=n<=H.         (8)
```

An isomorphism of such carriers must send the distinguished element `u` to
the distinguished element and commute with the readout.  The algebra need
not be presented as a quotient in advance.

Divide every polynomial `f` of degree at most `2D-1` uniquely as

```text
f=q ell_D+r,                    deg q<D,  deg r<D.       (9)
```

By `(5)`,

```text
L(f)=L(r).                                             (10)
```

Therefore

```text
A_D=Q[s]/(ell_D)                                      (11)
```

with readout equal to `L` on the unique representative of degree below `D`
is an exact `D`-dimensional multiplicative carrier for every factorial
moment through degree `2D-1`.

Only the rational quotient and identity `(10)` are used.  No splitting,
root-simplicity, positivity, or quadrature-node claim is needed.

## 3. The first failure is exact

Because `ell_D` is monic, `ell_D-s^D` has degree below `D`.  Thus `(5)` gives

```text
L(ell_D^2)
 =L(s^D ell_D)+L((ell_D-s^D)ell_D)
 =(D!)^2.                                               (12)
```

Call `h in (ell_D)` a **readout alias** when `L(h)!=0`: the quotient sends
`h` to zero but the intended readout does not.  Equation `(10)` shows that
no readout alias has degree at most `2D-1`.  But `ell_D^2=0` in the quotient
`(11)`, while `(12)` is nonzero.  Hence degree `2D` is the exact first
readout-alias degree:

```text
quotient readout of ell_D^2=0,
L(ell_D^2)=(D!)^2!=0.                                  (13)
```

In fact the whole first alias layer is rigid.  If
`h=q ell_D` has degree exactly `2D` and `a` is the leading coefficient of
`q`, then `q=a ell_D+r` with `deg r<D`, so

```text
L(h)=a(D!)^2!=0.                                       (13a)
```

Thus every ideal element of exact degree `2D` is a readout alias, while
none exists below it.  The horizon `2D-1` is sharp.

## 4. Optimal dimension and unique relation

The factorial Hankel block

```text
H_(D-1)=((i+j)!)_(0<=i,j<D)                            (14)
```

has determinant

```text
det H_(D-1)=product_(j=0)^(D-1)(j!)^2!=0.              (15)
```

For completeness, put
`Delta=diag(0!,1!,..., (D-1)!)` and
`P_(i,k)=binom(i,k)` for `k<=i`.  Vandermonde's identity gives

```text
H_(D-1)=Delta P P^T Delta.
```

The Pascal matrix `P` is unit lower triangular, proving `(15)` without an
analytic positivity input.

Indeed, a time-homogeneous linear carrier of dimension `R` has data
`(V,T,v,alpha)` with `dim V=R` and

```text
alpha(T^n v)=n!,                         0<=n<=2D-2.    (16)
```

For `0<=i,j<D`,

```text
(i+j)!=alpha(T^i(T^jv)),                                (17)
```

so `(14)` factors through the `R`-dimensional space `V`.  Its rank is `D`
by `(15)`, hence `R>=D`.  This proves representation-free optimality; the
carrier `(11)` attains the lower bound and is exact one degree farther.

There is also uniqueness in the precise multiplicative category of Section
2.  Let `(A,u,lambda)` have dimension `D` and horizon `2D-1`.  Evaluation
gives a surjection

```text
pi:Q[s]->A,                         pi(s)=u.             (18)
```

Because `A` is finite-dimensional, `ker(pi)` contains a nonzero monic
polynomial `q` of least degree `e<=D`.  For `0<=r<D`, exactness gives

```text
L(s^r q)=lambda(u^r q(u))=0.                            (19)
```

All degrees used in `(19)` are at most `2D-1`.  If `e<D`, the nondegenerate
block `(14)` forces the coefficient vector of `q`, padded to length `D`, to
vanish, a contradiction.  Thus `e=D`.  The same nondegeneracy makes the
monic degree-`D` solution of `(19)`
unique; `(5)` identifies it as

```text
q=ell_D.                                                (20)
```

It follows simultaneously that `1,u,...,u^(D-1)` are linearly independent.
They are therefore a basis of `A`, and `(18)` induces an algebra isomorphism

```text
Q[s]/(ell_D) ~= A                                      (21)
```

sending the distinguished class of `s` to `u`.  The values
`lambda(u^n)=n!` for `n<D` force this isomorphism to commute with the
readouts.  Thus `(11)` is unique up to the unique pointed
readout-preserving algebra isomorphism, not merely unique as a displayed
relation.

The last moment in the horizon is exactly what buys uniqueness.  Set
`ell_0=1`.  At the shortened horizon `2D-2`, every

```text
q_c=ell_D+c ell_(D-1),                         c in Q,   (21a)
```

is monic of degree `D` and satisfies

```text
L(s^r q_c)=0,                              0<=r<=D-2.  (21b)
```

Therefore `Q[s]/(q_c)`, with the same remainder readout, is exact through
degree `2D-2`.  These are all the possible monic degree-`D` relations at
that horizon: the `D-1` orthogonality equations have rank `D-1`, because
their first `D-1` coefficient columns form the invertible block
`H_(D-2)` (with the statement vacuous at `D=1`).  Finally,

```text
L(s^(D-1)q_c)=c((D-1)!)^2.                            (21c)
```

Thus the degree-`2D-1` observation kills exactly the residual `c` direction.
Dimension optimality begins at `2D-2`, uniqueness begins one moment later,
and the first readout failure of the unique carrier begins one moment after
that at `2D`.

This gives the complete horizon-by-horizon classification.  If a
dimension-`D` pointed `s`-generated carrier is exact only through `2D-2`,
a minimal relation of degree below `D` would still be orthogonal to
`1,s,...,s^(D-1)` (all resulting degrees are at most `2D-2`), contradicting
the rank of `(14)`.  Its minimal relation therefore has degree `D`, so the
preceding rank-`D-1` calculation makes it exactly one of the `q_c`.
Conversely every `q_c` works.  Also, for `D>=2`, the dimension-`D-1`
Laguerre quotient is exact through `2D-3`.  Hence for every horizon `H>=0`,

```text
minimum dimension in either carrier category=floor(H/2)+1;

H=2D-2:  the optimal pointed carriers are Q[s]/(ell_D+c ell_(D-1));
H=2D-1:  the unique optimal pointed carrier is Q[s]/(ell_D).       (21d)
```

The even/odd alternation is intrinsic: the even horizon forces the new
state dimension, and the following odd horizon fixes its last relation
parameter.

## 5. Exact inverse-Hankel selectors

For `d>=0`, put

```text
H_d=((i+j)!)_(0<=i,j<=d).                              (22)
```

The same Pascal factorization used for `(15)` gives

```text
(H_d^(-1))_(ij)
 =(-1)^(i+j)/(i!j!)
   sum_(k=max(i,j))^d binom(k,i)binom(k,j).             (23)
```

Indeed `H_d=Delta_d P P^T Delta_d`, where
`Delta_d=diag(0!,1!,...,d!)`, and the inverse Pascal entries are
`(-1)^(i-j)binom(i,j)`.

Define

```text
phi_j^(d)(s)=sum_(i=0)^d (H_d^(-1))_(ij)s^i.            (24)
```

Writing `L_k(s)=sum_i (-1)^i binom(k,i)s^i/i!`, the same coefficient
formula has the compact Laguerre-basis form

```text
phi_j^(d)(s)
 =(-1)^j/j! sum_(k=j)^d binom(k,j)L_k(s).
```

Expanding the right side and interchanging the `i,k` sums gives exactly
`(23)`.  Thus the optimal carrier relation and its coefficient duals belong
to one finite Laguerre system rather than being two unrelated Hankel
constructions.

Matrix inversion gives the exact duality

```text
L(s^r phi_j^(d))=delta_(rj),             0<=r,j<=d.     (25)
```

The selector norm is

```text
L((phi_j^(d))^2)
 =(H_d^(-1))_(jj)
 =(j!)^(-2) sum_(k=j)^d binom(k,j)^2.                   (26)
```

Consequently every `G` of degree at most `d` obeys

```text
[s^j]G=L(G phi_j^(d)).                                  (27)
```

In monomial-observation coordinates this says explicitly

```text
[s^j]G
 =sum_(i=0)^d (H_d^(-1))_(ij)L(s^iG).                  (27a)
```

For a radial shell with `s=ZW`, the required sidecar is therefore the
inserted bank `E[(ZW)^iP^m]`, `0<=i<=d`; the uninserted scalar moment is only
its `i=0` member.  Formula `(27a)` identifies the missing observations
without pretending that Gaussian nullity supplies them.

This is a complete finite-height coefficient decoder when the multiplier
observations on the right of `(27)` are available.

## 6. The multiplier-access boundary

Gaussian nullity supplies scalar equations `L(G_m)=0`; it does not supply
the family `L(G_m phi_j)`.  MISTAKE-211 is the sharp live witness.  With

```text
G_4=4ab^3s^6+4a^3cs^18,                                (28)
a=b=1,                         c=-6!/18!,               (29)
```

both coefficients in `(28)` are nonzero but

```text
L(G_4)=0.                                               (30)
```

The degree-eighteen selectors correctly recover

```text
L(G_4 phi_6^(18))=4,
L(G_4 phi_18^(18))=-1/2223046425600.                   (31)
```

Equations `(31)` do not follow from `(30)`.  Similarly, the quotient `(11)`
preserves the scalar factorial sum exactly; it does not make zero readout
imply a zero quotient class.

The connection contract is therefore:

| item | exact content |
|---|---|
| source | one finite radial shell `G_m` |
| optimal carrier | `Q[s]/(ell_D)` through height `2D-1` |
| exact dual | inverse-Hankel selectors `(24)` |
| preserved | scalar factorial readout |
| missing sidecar | inserted bank `L(s^iG_m)=E[(ZW)^iP^m]`, or algebra-valued nullity |
| sharp hostile | MISTAKE-211 two-height cancellation `(28)--(31)` |

The theorem sharpens THM-2810 but does not separate arbitrary Wick channels,
prove HYP-8765, or replace THM-2022.

## 7. Exact companion

Run

```bash
python 04-computation/gmc_optimal_finite_laguerre_carrier_thm2815.py
python -O 04-computation/gmc_optimal_finite_laguerre_carrier_thm2815.py
python 04-computation/gmc_optimal_finite_laguerre_carrier_independent_audit_thm2815.py
python -O 04-computation/gmc_optimal_finite_laguerre_carrier_independent_audit_thm2815.py
```

The primary pair byte-matches

```text
05-knowledge/results/gmc_optimal_finite_laguerre_carrier_thm2815.out.
```

The companion verifies the carrier and first failure for `D=1,...,12`,
all inverse-Hankel entries and selector norms through `d=10`, and the exact
degree-eighteen MISTAKE-211 recovery.  It has explicit exception gates, no
truth-bearing Python assertions, no floating point, and no scratch
dependency.  The universal statements are proved above rather than inferred
from these finite controls.

The independent hostile audit uses only Python integers and
`fractions.Fraction`, with no SymPy or other CAS.  It independently rebuilds
the Laguerre relation through `D=14`, solves the monic orthogonality system
from the closed inverse, verifies the exact quotient horizon and first
readout-alias layer, checks the complete short-horizon family
`ell_D+c ell_(D-1)` on hostile rational parameters, checks both orders of
the inverse identity and an integer Bareiss determinant through `d=12`,
checks the compact Laguerre-basis selector formula, and recovers both
nonzero MISTAKE-211 coefficients from the scalar-null degree-eighteen shell.
Its normal and optimized executions byte-match

```text
05-knowledge/results/gmc_optimal_finite_laguerre_carrier_independent_audit_thm2815.out.
```

The audit specifically challenged the category behind “unique,” the meaning
of “first alias,” and the unused node-diagonalization sentence.  The first
two are now typed in Sections 2--4; the third was removed rather than
importing an unnecessary squarefreeness theorem.  It also rechecked that
`(31)` uses observations absent from `(30)`.  The MISTAKE-211,
HYP-8765, arbitrary-Wick-channel, and no-new-GMC scope boundaries therefore
remain load-bearing.

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
