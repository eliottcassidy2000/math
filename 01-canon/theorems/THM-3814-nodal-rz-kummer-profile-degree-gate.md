---
id: THM-3814
title: "Nodal RZ Kummer profiles have a degree-thirteen entry gate"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
  INDEPENDENT HOSTILE AUDIT.  On the c=1 cubic pseudo-plane, a Darboux pair
  with the normalized nodal arm jet and no canonical terms beyond
  r g(e), z^2 f(e), and rz p(e) must enter a nonconstant 5/4 Kummer tower:
  p=alpha e^3v^5, q=lambda p, and kappa-lambda f=beta e^2v^4 with
  deg(v)>=2.  Thus deg(p)>=13 and deg(kappa-lambda f)>=10.  The asymmetric,
  constant-v, and linear-v branches are impossible; towers of degree at
  least two remain OPEN.  This is an exact profile gate, not a no-go for
  arbitrary higher canonical coefficients and not a counterexample claim.
source: jc_zero_debt_lift / nodal rz Kummer-tower lane, 2026-08-23
audit: >
  PROVISIONAL EXACT CANDIDATE.  The deterministic companion has 40 active
  gates checking the Poisson signs, unique monic z-normal reduction, the
  arm, Wronskian, and r^2z^2 buckets, both Kummer valuation patterns, the
  p=0 boundary, the constant-v polynomial divisions, the actual r^3 and
  z^2 source coefficients, quotient-ring reduction modulo 4lambda^2+3,
  the final residual coefficient 1/4, both linear-v arm remainders, the
  uniform rank-five source equation, its generic and t=0 compatibility,
  and the resulting exact unit ideal.  Normal and optimized replay,
  frozen-output equality, hashes, and documentation checks are required
  before audit promotion.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3812-nodal-arm-coefficient-second-normal-profile-nonentry
related:
  - THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
  - THM-3811-nodal-arm-bezout-law-for-cubic-pseudoplane-darboux-pairs
  - THM-3813-quartic-r-repairs-of-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_nodal_rz_kummer_gate_thm3814.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_nodal_rz_kummer_gate_thm3814.out
script_sha256: c8384518b8af8d5ef204e5954e3cf60b796fde82b5e0c34042293b90e36ae445
output_sha256: 03df8eed641f81967781d4641ba6bdc31d91eb6306f700ca8ab21cdbd1a75ee9
semantic_sha256: 2b2e0156b8d0dfb6cfd14c1f0cbd48d13ac79addff056e419f50c41ff5cda77c
hash_basis: raw LF bytes
---

# THM-3814 -- the first mixed nodal profile has a degree-thirteen gate

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
INDEPENDENT HOSTILE AUDIT.**  Let `k` be an algebraically closed field of
characteristic zero and put

```text
B=k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary `f,g,h,kappa,p,q in k[e]`, consider the first mixed extension
of the THM-3812 nodal profile,

```text
A=e^2-z/3+r g(e)+z^2 f(e)+rz p(e),
C=e^3-e-ez/2+r h(e)+z^2 kappa(e)+rz q(e).             (2)
```

If

```text
{A,C}=1,                                                (3)
```

then there are `lambda,alpha,beta in k`, with `alpha beta!=0`, and a
polynomial `v in k[e]` of degree at least two such that

```text
q=lambda p,
p=alpha e^3 v^5,
kappa-lambda f=beta e^2 v^4.                           (4)
```

In particular every pair `(2)` satisfying `(3)` would have

```text
deg p >= 13,            deg(kappa-lambda f) >= 10.      (5)
```

The theorem is a necessary entry gate.  It does not assert that the
nonconstant towers in `(4)` satisfy the remaining bracket equations.

## 1. Unique canonical reduction and three decisive buckets

The surface algebra has the monic presentation

```text
B=k[r,e][z]/(z^3-r^2e-r).                              (6)
```

It is therefore free over `k[r,e]` with basis `(1,z,z^2)`.  Reducing the
bracket from `(2)` by `(6)` loses no quotient information: `(3)` holds if
and only if every nonscalar coefficient of that unique normal form vanishes.

Three coefficients of `{A,C}-1` are

```text
[z] = (36e^2f-24e kappa-12f+1)/2,                     (7)

[r^3z] = 15e(pq'-qp'),                                (8)

[r^2z^2]
 = -3(-4efq'+4e kappa p'-5ep kappa'
      +5eqf'+2fq-2kappa p).                           (9)
```

Thus the arm bucket `(7)` gives

```text
(3e^2-1)f-2e kappa=-1/12,              hence f(0)=1/12. (10)
```

The factor `e` in `(8)` may be cancelled because `k[e]` is a domain.
Consequently

```text
pq'-qp'=0.                                             (11)
```

These equations are necessary identities in `k[e]`; no division by a
possibly vanishing profile has yet occurred.

## 2. The asymmetric branch is empty

First suppose `p=0`.  If also `q=0`, `(2)` is exactly the `c=1` profile
excluded by THM-3812.  Suppose instead that `q!=0`.  Equation `(9)` becomes

```text
5e q f'-4e f q'+2fq=0.                                (12)
```

Here `f!=0` by `(10)`.  In `k(e)`, equation `(12)` says

```text
(q^4/(e^2f^5))'=0.                                    (13)
```

The constant field of `k(e)` is `k`, so for some `gamma in k*`,

```text
q^4=gamma e^2f^5.                                     (14)
```

Unique factorization in `k[e]` now forces

```text
q=alpha e^3v^5,              f=beta e^2v^4             (15)
```

for nonzero constants `alpha,beta` and a polynomial `v`.  In particular
`f(0)=0`, contradicting `(10)`.  Therefore every survivor has

```text
p!=0.                                                  (16)
```

## 3. The symmetric branch is a 5/4 Kummer tower

By `(11)` and the characteristic-zero constant-field fact,

```text
q=lambda p                                             (17)
```

for some `lambda in k`.  Set

```text
W=kappa-lambda f.                                      (18)
```

The polynomial `W` is nonzero.  Otherwise `(10)` would say

```text
(3e^2-2lambda e-1)f=-1/12,                             (19)
```

which is impossible because the displayed quadratic is not a unit in
`k[e]`.  Substituting `(17)--(18)` into `(9)` gives the exact equation

```text
4eW p'-5epW'-2pW=0.                                   (20)
```

After division in `k(e)`, this is

```text
(p^4/(e^2W^5))'=0,
p^4=gamma e^2W^5                         (gamma in k*). (21)
```

For completeness, the UFD calculation behind `(4)` is local at every
irreducible factor.  If `a=ord_pi(p)` and `b=ord_pi(W)`, then `(21)` gives

```text
4a=5b                         if pi!=e,
4a=2+5b                       if pi=e.                 (22)
```

The nonnegative solutions are respectively

```text
(a,b)=(5t,4t),
(a,b)=(3+5t,2+4t),                 t>=0.               (23)
```

Collecting the irreducible factors with multiplicity `t` into `v` proves
the representation `(4)` (with the constant relation
`alpha^4=gamma beta^5`).

If `v` is nonconstant, `(4)` already gives the preliminary bounds
`deg p>=8` and `deg W>=6`.  Sections 4--5 rule out degrees zero and one,
which will sharpen these to `(5)`.  We begin with the constant case.

## 4. The constant tower collapses to one hostile normal form

Absorb the nonzero constant `v` into `alpha,beta`, so that

```text
p=alpha e^3,                 W=beta e^2.               (24)
```

Put `D=3e^2-2lambda e-1`.  Equation `(10)` becomes

```text
D f=2beta e^3-1/12.                                   (25)
```

The Euclidean remainder on division of the right side by `D` is

```text
4beta lambda/9-1/12
 + e(8beta lambda^2/9+2beta/3).                       (26)
```

Its two coefficients must vanish.  Since `beta!=0`, this forces

```text
4lambda^2+3=0,              beta=-lambda/4,            (27)

f=1/12-lambda e/6,
kappa=lambda/12+e/8-lambda e^2/4.                     (28)
```

Write `H=h-lambda g`.  Two further coefficients in the same unique normal
form now have a useful triangular order.  The pure `r^3` coefficient gives

```text
eH'-2H=(3e+4lambda)/(120alpha),                        (29)
```

and hence

```text
H=delta e^2-e/(40alpha)-lambda/(60alpha)               (30)
```

for a constant `delta`.  The pure `z^2` coefficient then gives

```text
D g=(1-2lambda e+48eH)/24.                            (31)
```

Reducing the Euclidean remainder in `(31)` modulo the already necessary
law `4lambda^2+3=0` yields

```text
(160alpha delta lambda+15alpha-6)/(360alpha)
 - e lambda(5alpha+4)/(60alpha).                      (32)
```

Because `(27)` makes `lambda!=0`, vanishing of `(32)` determines

```text
alpha=-4/5,                 delta=3lambda/16,           (33)

g=(3lambda e-1)/24,
h=(9lambda e^2-3e-lambda)/48.                         (34)
```

Thus every constant-`v` possibility has been forced to the single normal
form `(27)--(28),(33)--(34)`.  Substitute it into the *full* bracket and
reduce by `(6)`.  In the quotient `k[lambda]/(4lambda^2+3)`, the coefficient
of the pure monomial `r e^0` is

```text
[r z^0 e^0]({A,C}-1)=1/4.                             (35)
```

Characteristic zero makes `(35)` nonzero, contradicting `(3)`.  Hence `v`
cannot be constant.

## 5. The linear tower is empty

It remains possible at this stage that `deg v=1`.  Absorb its leading
coefficient into `alpha,beta` and write

```text
v=e-t,
p=alpha p_0,       p_0=e^3(e-t)^5,
W=beta e^2(e-t)^4.                                   (36)
```

The arm equation is again

```text
D f=2eW-1/12,             D=3e^2-2lambda e-1.         (37)
```

Exact Euclidean division gives

```text
rem_D(2eW-1/12)=A_0/2916 + e(2beta A_1/729),          (38)
```

where

```text
A_0=
 256beta lambda^5-1536beta lambda^4t
 +3456beta lambda^3t^2+768beta lambda^3
 -3456beta lambda^2t^3-3456beta lambda^2t
 +1296beta lambda t^4+5184beta lambda t^2+432beta lambda
 -2592beta t^3-864beta t-243,                         (39)

A_1=
 64lambda^6-384lambda^5t+864lambda^4t^2+240lambda^4
 -864lambda^3t^3-1152lambda^3t
 +324lambda^2t^4+1944lambda^2t^2+216lambda^2
 -1296lambda t^3-648lambda t+243t^4+486t^2+27.        (40)
```

Thus polynomiality of `f` requires `A_0=A_1=0`.

There is a small exact compatibility behind the remaining equations.  Put
`H=h-lambda g` and `X=alpha H`.  After division by its nonzero scalar factor
`-3`, the pure `r^3` bucket of the full canonical bracket is equivalent to

```text
4e^2(Wf'-fW')
 +3eX p_0'-5e p_0X'+X p_0=0.                         (41)
```

The forcing term has degree at most twelve.  Since `p_0` is monic of degree
eight, a leading term `x_de^d` in `X` contributes
`(25-5d)x_de^(d+8)` to the homogeneous part of `(41)`.  It follows that
`deg X<=5` (degree five is the unique resonance).

On the six-dimensional space `deg X<=5`, the coefficient matrix of the
homogeneous operator in `(41)` has rank exactly five for every `t`.  Indeed
`e^2(e-t)^3` spans a visible kernel, while the minor in coefficient rows
`e^8,...,e^12` and columns `1,e,...,e^4` is the nonzero constant
`375000`.  A left-null compatibility for `t!=0`, together with the separate
`t=0` coefficient, gives in both cases

```text
C=(-2lambda+3t)
  (16lambda^4-48lambda^3t+36lambda^2t^2+48lambda^2
   -90lambda t+27t^2+27)=0.                           (42)
```

More explicitly, for `t!=0` the functional with nonzero coefficient-row
entries

```text
ell_3=-1176/t^6, ell_4=-224/t^5, ell_6=20/t^3,
ell_7=10/t^2,    ell_8=4/t,       ell_9=1
```

annihilates the rank-five operator, and its value on the right side of
`(41)` is `32beta^2 C/(729t^2)`.  When `t=0`, the `e^7` coefficient of
`(41)` is instead `32beta^2 C/486`.  Since `beta!=0`, no denominator case
is lost and `(42)` is necessary uniformly.

The three necessary equations are already inconsistent.  Introducing
`beta_inv` only to encode `beta!=0`, exact characteristic-zero Buchberger
reduction in the displayed variable order gives

```text
<A_0,A_1,C,beta beta_inv-1>
   = <1> in Q[beta_inv,beta,lambda,t].                 (43)
```

The companion constructs `(38)--(43)` from the bracket, rather than taking
the packet as input, and freezes the three-polynomial packet by a separate
hash.  Therefore `deg v=1` is impossible.  Together with Section 4, this
proves `deg v>=2` and hence `(4)--(5)`.

## 6. Exact boundary and next construction target

This theorem closes the first mixed layer only up to a sharp structural
boundary.  It proves that neither one-sided `rz` repair nor a constant or
linear Kummer profile can work.  It does **not** close

```text
p=alpha e^3v^5,
q=lambda p,
kappa-lambda f=beta e^2v^4,       deg v>=2.            (44)
```

and it does not allow still higher canonical coefficients to be discarded.
The cheapest remaining internal test is therefore `deg v=2`, beginning at
degrees `(deg p,deg W)=(13,10)`, while the orthogonal escape is to add the
next canonical `r,z` monomial layer.  Notice also that `rz in (r,z)^4` on
this surface, as typed in THM-3812; it is the first omitted *mixed canonical
coefficient*, not an `I^3` correction.

The exact companion named in the metadata derives the entire arbitrary-
function bracket, verifies `(7)--(9)`, checks both valuation patterns in
`(23)`, freezes the divisions `(26),(32),(38)`, recomputes `(35)` from the
full normal form, and proves the uniform rank and compatibility in
`(41)--(43)`.  Quotient-ring reductions reject any `lambda`-dependent
denominator before reducing modulo `(27)`.  There is no finite-field or
bounded-degree inference in the proof.  **QED.**
