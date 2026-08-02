---
id: THM-3089
title: "Square-root moving-gap cluster cone and uniform alternant conditioning"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.  The
  fixed-gap cluster theorem remains physical uniformly for arbitrary
  distinct moving gaps of diameter o(sqrt C), and in a sufficiently small
  fixed square-root cone.  A total-positivity compactness lemma and the
  regular-bipartite character of every resultant monomial convert the exact
  coefficient error O(H^2/C) into a relative resultant error with no
  exponential condition loss.  At H asymptotic to sqrt C a Gaussian
  multi-index deformation is the exact new face.  Every covered multi-normal
  suspension is eventually positive.  This is fixed-width and
  does not prove positivity on that boundary or for arbitrary supports.
source: root-gmc-moving-gap-cluster-2026-08-01
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
  - THM-3085-multi-normal-fixed-gap-cluster-and-unconditional-all-width-tail
  - THM-3086-arbitrary-cluster-composition-chambers-and-alternant-clutch-holotopy
related:
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
script: 04-computation/gmc_logarithmic_moving_gap_cluster_thm3089.py
output: 05-knowledge/results/gmc_logarithmic_moving_gap_cluster_thm3089.out
script_sha256: 765a414dfa5672b9391325d14dfe6de6ab51e41bc3cccaeb7552f5db8c31bb77
output_sha256: cf5b3db020e61f11803cb5eaa46344fbf23cb025fbe852746d0fa79d31e89290
hash_basis: LF-normalized bytes
---

# THM-3089 -- square-root moving-gap cluster cone

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3085 fixes every internal gap of the remote cluster.  THM-3086 retains
the resulting generalized alternants only as fixed-gap physical sidecars and
as growing-gap algebraic symbols.  The fixed-gap hypothesis is not the real
boundary.  The exact factorial distortion is quadratic in the cluster
diameter, while the positive generalized alternant has a uniformly
controlled absolute envelope.  This makes a full square-root moving-gap
cone physical.

## 1. Statement

Fix a physical lower support of width `m>=3` with normalized first-window
resultant

```text
S_m!=0.                                                     (1)
```

Fix `q>=2`, put

```text
k=m+q,                    p=m+1,                           (2)
```

and for every sufficiently large integer `C` choose arbitrary distinct
integer gaps

```text
0=h_1(C)<...<h_q(C),      H(C)=h_q(C).                    (3)
```

Append the physical offsets

```text
C+h_1(C),...,C+h_q(C).                                    (4)
```

Let `R_k(C)` be the normalized physical resultant, let `U_(r,C)>0` be the
exact THM-3085 carrier for `p<=r<=k`, and let `E_C` be the exact intrinsic
all-high `q`-normal resultant.  Put

```text
D_z=p(p+1)...k=k!/m!,
rho_m=(m/(m+1))^m.                                        (5)
```

There is a constant `kappa_(m,q)>0`, depending only on the fixed width, such
that if

```text
H(C)<=kappa_(m,q) sqrt(C),                                (6)
```

then, uniformly over every distinct gap vector `(3)`,

```text
R_k(C)=S_m^D_z E_C^m!
       product_(r=p)^k U_(r,C)^[k!/r]
       [1+O(poly(C,H)e^(c_(m,q)H) rho_m^C)],              (7)
```

and the exact normal factor has the sign of its generalized-alternant model.
Consequently

```text
R_k(C)>0.                                                 (7a)
```

for all sufficiently large `C`, uniformly in `(3)`.  Indeed `q>=2` makes
`D_z` even, `m!` is even, and every `U` factor is positive, so both `S_m^D_z`
and `E_C^m!` are positive regardless of the signs of `S_m,E_C`.

The asymptotic version is sharper.  If

```text
H(C)=o(sqrt(C)),                                          (8)
```

then, writing

```text
M_C=(r^h_j(C))_(p<=r<=k,1<=j<=q),                        (9)
```

one has the relative estimate

```text
E_C=(-1)^[qD_z] det(M_C)^D_z
       [1+O_(m,q)(H(C)^2/C)].                            (10)
```

In particular `(10)` tends relatively to one.  The error in `(7)` also
tends to zero under either `(6)` or `(8)` because
`exp(O(sqrt C))rho_m^C->0`.

When `m=3`, THM-2824 gives `S_3>0` for every arbitrary three-slot physical
base.  Therefore every fixed width `K>=5`, every such base, and every
arbitrary distinct remote gap vector with diameter `o(sqrt C)` have positive
first-window resultant for all sufficiently large `C`.  This is an
unconditional physical moving-cluster tail, not merely a formal symbol.

## 2. Exact factorial distortion

Put `N=n+C` and use the rising factorial

```text
(x)_s=x(x+1)...(x+s-1).                                  (11)
```

In the direct high-function coordinates, a coefficient of the all-high
degree-`r` form is indexed by

```text
alpha=(alpha_1,...,alpha_q),
sum alpha_j=r,
s_alpha=sum alpha_j h_j.                                 (12)
```

After dividing by `U_(r,C)`, its exact/model ratio is

```text
Q_alpha/r^s_alpha,
Q_alpha=(rN+1)_(s_alpha) /
        product_j (N+1)_(h_j)^alpha_j.                   (13)
```

Indeed `(13)` is the literal ratio of the factorial moments after pulling
out `(rN)^s_alpha/N^s_alpha=r^s_alpha`.  Since
`log(1+x)<=x`, one has uniformly in `alpha`, `r`, and the gap vector

```text
|log(Q_alpha/r^s_alpha)|
 <=sum_(u=1)^s_alpha u/(rN)
   +sum_j alpha_j sum_(u=1)^h_j u/N
 <=r H(H+1)/N.                                          (14)
```

Thus every direct-basis coefficient differs relatively from its line-power
model by

```text
O_(m,q)(H^2/C).                                          (15)
```

Make that transformation explicit.  With direct variables
`v=(v_1,...,v_q)` and normal variables `u=(w,z_1,...,z_(q-1))`, set

```text
v_j=z_j                  (j<q),
v_q=-w-sum_(j<q)z_j.                                    (15a)
```

Its determinant is `(-1)^q`, independently of `C` and the gaps, and

```text
sum_j r^h_j v_j
 =-r^h_q w+sum_(j<q)(r^h_j-r^h_q)z_j.                  (15b)
```

which is exactly THM-3085's row.  Both exact and model systems undergo
`(15a)` as whole forms, so resultant covariance contributes
`(-1)^[qD_z]`.  No termwise relative estimate is made after this subtractive
change; cancellations there can create zero model coefficients.

## 3. Uniform generalized-alternant condition lemma

Let fixed positive rows satisfy

```text
0<r_1<...<r_q,                                            (16)
```

and for an arbitrary integer gap vector `0=h_1<...<h_q` put

```text
M_h=(r_i^h_j),
per(M_h)=sum_pi product_i r_i^h_(pi(i)).                  (17)
```

Then there is a finite constant `C_(r,q)` such that

```text
1<=per(M_h)/det(M_h)
  <=q! product_i r_i^(i-1) /
        product_(i<j)(r_j-r_i)=C_(r,q)                   (18)
```

for **every** gap vector.  The determinant is positive by strict total
positivity and is a nonzero integer, but the uniform upper bound is the
load-bearing statement.

Write

```text
d_i=h_i-(i-1),
V(r)=product_(i<j)(r_j-r_i).                              (19)
```

The `d_i` are nonnegative and increasing.  The generalized-Vandermonde
factorization is

```text
det(M_h)=V(r)s_lambda(r_1,...,r_q),                       (20)
```

where `lambda` is the reverse of `(d_1,...,d_q)`.  Schur positivity includes
the monomial `product_i r_i^d_i` with coefficient at least one.  Therefore

```text
det(M_h)>=V(r) product_i r_i^[h_i-(i-1)].                 (21)
```

Rearrangement says every permanent monomial is at most the increasing
matching `product_i r_i^h_i`, so

```text
per(M_h)<=q! product_i r_i^h_i.                           (22)
```

Dividing `(22)` by `(21)` proves the explicit upper bound `(18)`.

The lemma is stronger than an integer lower bound: it controls all absolute
monomials by the signed alternant uniformly even when the gaps diverge.

## 4. The regular-bipartite resultant character

Let `P` be the universal homogeneous resultant of `q` forms of degrees
`p,...,k` in `q` variables.  Its degree in the coefficients of the
degree-`r` form is

```text
mu_r=D_z/r.                                               (23)
```

Substitute the positive line powers associated to `M_h` and fully expand one
resultant monomial in the entries of `M_h`.  Let `e_(r,j)` be its exponent
matrix.  Each row sum is

```text
r mu_r=D_z,                                               (24)
```

while diagonal variable covariance of the resultant makes every column sum
also `D_z`.  Thus `e` is the incidence matrix of a `D_z`-regular bipartite
multigraph.  Konig's line-colouring theorem, equivalently the integral
Birkhoff theorem, decomposes it into `D_z` permutation matrices.

Every expanded monomial consequently satisfies

```text
product_(r,j) M_(r,j)^e_(r,j)
 <=per(M_h)^D_z.                                         (25)
```

The universal polynomial has a finite coefficient bank depending only on
`m,q`.  If `P_abs` denotes the sum with absolute coefficients, `(18),(25)`
give

```text
P_abs(M_h)<=C_(m,q) det(M_h)^D_z.                        (26)
```

Let

```text
T_z=sum_(r=p)^k mu_r                                     (27)
```

be the total coefficient degree.  Telescoping each coefficient product and
using `(15)` and `(26)` gives, if `epsilon=O_(m,q)(H^2/C)`, the more explicit
relative envelope

```text
C_0[exp(T_z epsilon)-1]
   [per(M_C)/det(M_C)]^D_z.                              (28)
```

The uniform condition lemma bounds the last factor.  Hence

```text
|E_C-(-1)^[qD_z]det(M_C)^D_z|
 <=C'_(m,q)(H^2/C) det(M_C)^D_z.                         (29)
```

This is `(10)`.  Choose `kappa_(m,q)>0` so that the exact envelope `(28)` is
less than one half whenever `(6)` holds.  This proves the
uniform sign claim without an exponential condition-number loss.

## 5. The physical contraction and composition extension

THM-3085's nonsurviving coefficient layers still have exponential base at
most `rho_m^C`.  Here is a uniform invoice with no hidden `C^H`.  In a layer
with `j` high entries, let `A` be the total offset supplied by the remaining
slower entries and let `S<=jH` be the sum of the moving internal gaps.  The
exact extra factorial ratio relative to replacing those high entries by
their cluster base has the form

```text
(rn+1+A+jC)_S /
 product_(high entries ell)(n+C+1)_(h_ell).              (29a)
```

At a local node `A+jC<=rC`.  Once `C>=max(n,H)`, every numerator factor in
`(29a)` is at most `2r(n+C)` and every denominator factor is at least
`n+C`.  Hence `(29a)<=(2r)^S<=(2k)^(kH)`.  After the finite signed inclusion
sum, every nonsurviving layer is therefore bounded by

```text
poly(C,H)exp(c_(m,q)H)rho_m^C.                           (29b)
```

Under `(6)` this is subexponential in `C`, so the strict `rho_m` gap
survives.  The comparison is made before any signed layer cancellation.

The exact lower-error quotient is unchanged.  THM-3073 gives

```text
Res=S_m^D_z E_C^m!,                                      (30)
```

and THM-3085's variable/equation covariance supplies the `U` product in
`(7)`.  Since the exact `E_C` is retained in the carrier, `(29)--(30)` and
coefficientwise resultant continuity prove `(7)`.

The same argument crosses THM-3086's former growing-gap boundary.  At each
node of a fixed cluster composition, permit arbitrary distinct internal gaps
of diameter

```text
H_i(T)=o(sqrt(L_i(T))),                                   (31)
```

or choose a sufficiently small nodewise square-root constant.  Keep the
same strict entropy chamber

```text
J_(m_i)(delta_i)<Gamma_i.                                (32)
```

The new factors in `(31)` are `exp(o(L_i))`, so they do not change the
`Gamma_i-J` margin.  Formula THM-3086 `(18)` remains valid with the exact
moving `E_i(T)`, every such factor has its generalized-alternant sign, and
the final resultant is positive.  Taking the minimum of the finitely many
prefix constants preserves THM-3086's simultaneous positivity at every
individual support width.

## 6. The square-root Gaussian face

The exponent `1/2` is the exact boundary of the pure-line approximation.
Suppose

```text
h_j(C)=t_j sqrt(C)+o(sqrt(C)),                            (33)
```

with fixed `t_j>=0`.  Expanding `(13)` one order farther gives

```text
log(Q_alpha/r^s_alpha)
 ->(sum_j alpha_j t_j)^2/(2r)
   -(1/2)sum_j alpha_j t_j^2

  =-(1/(2r))sum_(i<j) alpha_i alpha_j(t_i-t_j)^2.        (34)
```

If two distinct slopes occur in the support of `alpha`, `(34)` is strictly
negative.  The limiting all-high form is therefore a genuine
Gaussian/Hadamard deformation of the pure line power.  Its logarithm is
quadratic, not linear, in `alpha`, so it cannot be absorbed into a single
linear coordinate change.

This is the first new physical face.  It proves sharp failure of the
coefficientwise pure-line mechanism at square-root scale; it does **not**
prove that the deformed resultant vanishes or changes sign.  A mesoscopic
normal-cone or heat-semigroup analysis might extend positivity across this
face.

## 7. Scope and exact evidence

The theorem fixes `m,q` and the lower support.  The constant in `(6)` is not
uniform in width.  Repeated gaps kill `det(M_C)`.  Diameters comparable to
`sqrt C` with arbitrary constant, or larger diameters, are not covered.
The result does not prove arbitrary equal-scale supports, a maximal Newton
fan, an all-order Stieltjes tail, arbitrary-radial GMC(2), NC2, LRC(14),
JC(2), or DC(2).

The initially reserved logarithmic route used only `|det M|>=1` and paid an
exponential absolute-condition invoice.  Sections 3--4 replace it by the
structured permanent bound; generic `N-(N-1)` cancellation is therefore not
a physical hostile for this resultant family.  The true boundary is the
Gaussian face `(34)`.

The exact companion verifies:

1. `3,850` exact rising-factorial ratios against the rational logarithmic
   bound `(14)`;
2. `31,930` positive integer alternants through rank five and gap height
   thirty, with permanent/determinant maxima
   `9,375,168268/3,241529660/9` at consecutive gaps and the explicit Schur
   condition bound `(18)` in every cell;
3. `1,361` regular exponent matrices, integral permutation decompositions,
   and the monomial bound `(25)`;
4. sixty exact direct-to-normal transformations, including determinant
   `(-1)^q` and every displayed row in `(15b)`;
5. twenty-four resultant degree invoices and prefix monotonicity controls;
6. the width-twelve value `T_z=101378880`; and
7. `387` exact identities `(34)`, including `351` strict mixed-sheet
   Gaussian cells.

Run

```text
python 04-computation/gmc_logarithmic_moving_gap_cluster_thm3089.py
python -O 04-computation/gmc_logarithmic_moving_gap_cluster_thm3089.py
```

Both modes must equal the stored transcript after LF normalization.

**QED, pending independent audit and status promotion.**
