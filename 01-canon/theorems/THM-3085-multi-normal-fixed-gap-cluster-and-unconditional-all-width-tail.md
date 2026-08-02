---
id: THM-3085
title: "Multi-normal fixed-gap cluster and unconditional all-width tail"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.  Appending
  any finite set of distinct fixed gaps at one remote scale to a lower
  physical support gives a multi-normal suspension whose only obstruction is
  the lower resultant.  Its exact normal symbol is a positive generalized
  Vandermonde, and its physical contraction has gap (m/(m+1))^m.  Since every
  arbitrary three-slot factorial support has positive resultant, every such
  base plus an arbitrary fixed remote cluster has eventual positive
  first-window resultant in every width, uniformly with gap 27/64 for each
  fixed width.  This is not an arbitrary-support or moving-gap theorem.
source: root-gmc-multinormal-cluster-2026-08-01
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3063-terminal-suspension-transverse-resultant-and-five-slot-tail-holotopy
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
related:
  - THM-3082-admissible-suspension-word-simultaneous-chambers-and-scale-tree-holotopy
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
script: 04-computation/gmc_multinormal_fixed_gap_cluster_thm3085.py
output: 05-knowledge/results/gmc_multinormal_fixed_gap_cluster_thm3085.out
script_sha256: ce2215af6c9f8b48edeb2e8ce5b0333a5cccd95a9e1b3a1a258b35f9d152cb15
output_sha256: b0a024d52d930764bbfb9b8c8018a31b94cb14770c41fe70c019b170a2e3b30d
hash_basis: LF-normalized bytes
---

# THM-3085 -- multi-normal fixed-gap cluster and unconditional all-width tail

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3063 leaves arbitrary six-slot SFC blocked when a remote **pair** is
attached to an arbitrary four-slot base.  Change the cut: keep only an
arbitrary three-slot base and put every remaining slot into one fixed-gap
remote cluster.  The normal space can have any finite dimension, and its
symbol is an exact generalized Vandermonde.  This yields an unconditional
all-width physical tail.

## 1. Conditional multi-normal suspension

Fix `n>=1` and a physical lower support of width `m>=3`, with normalized
first-window resultant

```text
S_m!=0.                                                     (1)
```

Put `q>=1`, `k=m+q`, `p=m+1`, and fix distinct integer gaps

```text
0=h_1<h_2<...<h_q.                                         (2)
```

Append the offsets

```text
C+h_1,...,C+h_q                                            (3)
```

and let `R_k(C)` be the normalized physical resultant of the degree
`2,...,k` first-window system.  The lower support and all gaps are fixed while
`C->infinity`.

For `p<=r<=k`, use the positive carriers

```text
U_(r,C)=L(f_(n+C)^r)/L(f_n^r)>0.                           (4)
```

There is an exact intrinsic `q`-variable all-high block `A_(r,C)` in the
degree-`r` equation, normalized by `U_(r,C)`.  Define its normal resultant

```text
E_C=Res_z(A_(p,C),...,A_(k,C)).                             (5)
```

Then

```text
R_k(C)=S_m^[k!/m!] E_C^[m!]
       product_(r=p)^k U_(r,C)^[k!/r]
       [1+O(poly(C) rho_m^C)],

rho_m=(m/(m+1))^m<1.                                      (6)
```

In particular `R_k(C)` is eventually nonzero.  Its eventual sign is
`sign(S_m)^[k!/m!]`; the normal factor in `(6)` is positive after the even
power `m!`.

## 2. Fixed-pivot coordinates and the normal symbol

Choose one fixed lower offset `c` and use the determinant-one difference
coordinates

```text
(y_1,...,y_(m-1),z_1,...,z_q),                             (7)
```

pivoted on `c`, never on a moving cluster point.  Restriction to `z=0`
recovers the lower physical forms.  In the all-high layer, the normalized
degree-`r` normal form tends

```text
A_(r,C) -> L_r(z)^r,

L_r=(-r^h_q,
     r^h_1-r^h_q,...,r^h_(q-1)-r^h_q),
                  p<=r<=k.                                (8)
```

Let `M=(r^h_j)` have rows `r=p,...,k` and columns `j=1,...,q`.
Elementary column operations give

```text
det(L_r)=(-1)^q det M.                                     (9)
```

Because the rows and integer exponents in `(2)` are strictly increasing,

```text
det M=product_(i<j)(r_j-r_i) s_lambda(r_p,...,r_k)>0,       (10)
```

where `s_lambda` is the Schur polynomial associated to the shifted gap
partition.  This is the generalized-Vandermonde positivity theorem; all its
monomial coefficients are nonnegative and it is nonzero on positive inputs.

The resultant of powers of independent linear forms is

```text
Res(L_p^p,...,L_k^k)=det(L)^[D_z],
D_z=product_(r=p)^k r=k!/m!.                              (11)
```

This follows from common linear covariance and calibration on coordinate
powers.  Consequently

```text
E_C=(-1)^[qD_z](det M)^D_z[1+O(C^-1)]!=0.                 (12)
```

If `q>=2`, `D_z` contains two consecutive integers and is even, so `E_C` is
itself eventually positive.  If `q=1`, its sign can be `(-1)^k`; the exponent
`m!` in `(6)` is even and removes it.  Repeated gaps give repeated columns in
`M` and are the exact leading-face-zero boundary.

## 3. The coefficientwise contraction

Scale every normal coordinate by

```text
lambda_C=U_(p,C)^(-1/p),                                  (13)
```

and multiply each upper equation of degree `r>=p` by

```text
s_(r,C)=U_(p,C)^(r/p)/U_(r,C).                             (14)
```

The variable covariance contributes `lambda_C^(q k!)`; equation-block
homogeneity contributes `product s_(r,C)^(k!/r)`.  The powers of `U_p`
cancel exactly:

```text
lambda_C^(qk!) product_(r=p)^k s_(r,C)^(k!/r)
 =product_(r=p)^k U_(r,C)^(-k!/r).                        (15)
```

Thus the transformed resultant is exactly `R_k(C)` divided by the carrier
product in `(6)`.

Consider one coefficient layer containing `ell` normal variables and `j`
actual high-cluster factors.  Necessarily `ell>=j`; fixed-pivot terms can have
`j=0<ell`.  Ignoring only polynomial and fixed-gap constants, its exponential
base after `(13)--(14)` is

```text
j^j/p^ell                         for degree r<=m,
j^j p^(r-ell)/r^r                 for degree r>=p,          (16)
```

with `0^0=1`.

The only base-one survivors are

```text
(r<=m,ell=j=0):            H_r,
(r=p,ell=j=0):             H_p,
(r>=p,ell=j=r):            A_(r,C).                        (17)
```

Every other cell has base at most `rho_m`.

For the lower cells, `ell>=j` gives `(j/p)^j`; convexity reduces its maximum
on `1<=j<=m` to an endpoint, and

```text
1/p <=(m/p)^m=rho_m.                                      (18)
```

The `j=0,ell>=1` cells satisfy the same bound.  For upper cells with `j=0`,
either one has the declared `r=p,ell=0` survivor or the base is bounded by
`rho_m`.  For `j>0`, use `ell>=j`; convexity in `j` reduces to the zero end or
`j=r-1`, where

```text
(p/r)((r-1)/r)^(r-1)<=rho_m,                              (19)
```

with equality at `r=p,j=ell=p-1`.  Monotonicity in `r` proves the remaining
cases.  The lower equality occurs at `r=m,j=ell=m`.  These are the sharp
`rho_m` cells.

Therefore the entire transformed physical system contracts coefficientwise
to

```text
H_2,...,H_m, H_p+A_(p,C), A_(p+1,C),...,A_(k,C)
 +O(poly(C)rho_m^C).                                      (20)
```

This is a whole-form contraction on one fixed projective space, not a
termwise sign comparison of Macaulay monomials.

## 4. Exact separated resultant

The lower block in `(20)` is pure in `y`, while the upper block is arbitrary
in `(y,z)` and restricts at `y=0` to the exact normal system in `(5)`.
THM-3073 applies with

```text
D_y=2*3*...*m=m!,
D_z=p*(p+1)*...*k=k!/m!.                                  (21)
```

It gives, with canonical positive normalization,

```text
Res(20)=S_m^D_z E_C^D_y.                                  (22)
```

Equations `(12)`, `(15)`, coefficientwise continuity of the fixed-degree
resultant, and `(20)--(22)` prove `(6)`.  Exact quotienting is what makes the
arbitrary child form `H_p` harmless; it is not required to converge to zero.

## 5. Unconditional all-width remote clusters

Now set `m=3`.  For every arbitrary triple `(a,b,c)`, THM-2824 gives

```text
S_3>0.                                                      (23)
```

Hence for every fixed `K>=4` and every fixed gap set

```text
0=h_1<...<h_(K-3),                                        (24)
```

the actual support

```text
(a,b,c,C+h_1,...,C+h_(K-3))                               (25)
```

has positive first-window resultant for all sufficiently large `C`.  More
precisely,

```text
R_K(C)=S_3^[K!/6] E_(4..K,C)^6
       product_(r=4)^K U_(r,C)^[K!/r]
       [1+O(poly(C)(27/64)^C)]>0.                         (26)
```

The exponential gap `27/64` is independent of `K`; the implicit polynomial,
threshold, and normal determinant still depend on the fixed width and gaps.
Thus `(26)` is not a width-uniform numerical threshold.

The conceptual gain is substantial at width six.  A remote pair at `k=6`
passes an arbitrary four-slot obstruction into THM-3063.  A remote **triple**
instead cuts transversely back to the universally positive three-slot cell,
so no arbitrary four-slot hypothesis appears.

## 6. Full asymptotic and the six-slot triple

Gauss multiplication gives

```text
U_(r,C)~kappa_(r,n) r^(rC) C^(-(r-1)/2),
kappa_(r,n)>0.                                             (27)
```

Let

```text
alpha_K=(K!/2) sum_(r=4)^K (r-1)/r.                       (28)
```

Since `(12)` raised to the sixth power contributes `(det M)^(K!)`, equation
`(26)` yields

```text
R_K(C)~A_(K,h) (K!/6)^(K! C) C^(-alpha_K),
A_(K,h)>0.                                                 (29)
```

For `K=6` and gaps `(0,h,g)`, `0<h<g`, the normal symbol is

```text
det [[1,4^h,4^g],
     [1,5^h,5^g],
     [1,6^h,6^g]]^720>0.                                  (30)
```

At the consecutive gaps `(0,1,2)`, the determinant is exactly two, so

```text
R_6(C)~A S_3^120 2^720 120^(720C) C^-858,
A>0.                                                       (31)
```

The exact carrier is `(26)`; `(29)--(31)` replace it only by its leading
inverse-power symbol.

## 7. Boundaries and exact evidence

The theorem fixes the child support, width, and all high gaps.  Repeated gaps
kill the generalized Vandermonde; moving gaps or multiple positive scales
require THM-3082's entropy chambers; a zero lower resultant removes the
carrier; arbitrary high offsets need not lie in any fixed cluster.  The
theorem does not prove arbitrary support SFC, a uniform threshold in width or
gaps, a maximal fan, an all-order Stieltjes tail, arbitrary-radial GMC(2),
NC2, LRC(14), JC(2), or DC(2).

The exact companion verifies:

1. all `5,296` nonsurviving lower/upper layer cells for `3<=m<=10` and
   degrees through fourteen, the `112` declared survivors, and all `16`
   sharp `rho_m` equalities;
2. `381` generalized-Vandermonde and normal-sign controls for ranks two
   through six and anchored gap subsets of `{0,...,9}`;
3. sixty conditional covariance/exponent ledgers and every unconditional
   width `4<=K<=14`;
4. the carrier exponents in `(6)` and `(26)`; and
5. the exact `K=6` determinant, exponential base, symbol exponent, and
   inverse-power exponent in `(31)`.

Run

```text
python 04-computation/gmc_multinormal_fixed_gap_cluster_thm3085.py
python -O 04-computation/gmc_multinormal_fixed_gap_cluster_thm3085.py
```

Both modes must equal the stored transcript after LF normalization.

**QED, pending independent audit and status promotion.**
