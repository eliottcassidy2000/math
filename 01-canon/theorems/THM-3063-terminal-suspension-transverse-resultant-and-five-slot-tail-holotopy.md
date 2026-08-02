---
id: THM-3063
title: "Terminal suspension, transverse resultant, and the five-slot physical tail"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every k>=3,
  a fixed-low physical factorial cell with one
  moving terminal pair has an eventual positive resultant whenever its
  transverse (k-2)-slot resultant is nonzero.  The exact leading suspension
  is S^[k(k-1)] times the two top normal forms, and every fixed affine finite
  generalized-Hankel bank is uniformly positive along a carrier-to-physical
  holotopy.  For k=5 the transverse cell is an arbitrary three-slot cell, so
  THM-2824 proves the hypothesis unconditionally.  The first new obstruction
  is k=6, exactly the arbitrary four-slot first-window resultant.
source: kind-pasteur-2026-08-01-terminal-suspension
audit: >
  Two independent hostile audits rederived the determinant-one coordinate
  covariance, every scaled layer base including fixed-layer leakage and the
  k=3 endpoint, independence from H_(k-1), the separated resultant sign and
  exponents, the five-slot constants, and the carrier holotopy.  Both replayed
  normal, optimized, and stored companions, matched the declared LF hashes,
  and cross-checked the executable all-high block against the literal physical
  top layers.  Their audit replaced the earlier invalid separate-line
  coercivity argument by the current isotropic whole-form contraction.
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-3040-formal-corner-resultant-width-quotient-and-all-order-bernoulli-law
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3054-affine-moving-lower-tropical-beta-gamma-tail-holotopy
  - THM-3062-four-slot-physical-transverse-augmentation-and-affine-tail-holotopy
related:
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-3060-three-slot-physical-terminal-face-and-affine-tail-holotopy
script: 04-computation/gmc_terminal_suspension_five_slot_tail_thm3063.py
output: 05-knowledge/results/gmc_terminal_suspension_five_slot_tail_thm3063.out
script_sha256: 6105bddcc05a776272b8a6b149017e589de8c06a19bcdc07f29075fb749a4543
output_sha256: b08563d21d691ba5bdeeaf019e678f28ffeaffb29de98025ae128625e166ac1b
hash_basis: LF-normalized bytes
---

# THM-3063 -- terminal suspension and the five-slot physical tail

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3060 and THM-3062 expose the same phenomenon at consecutive transverse
depths.  The all-large physical face always has only two normal directions.
The information apparently lost in the other coordinates is not arbitrary:
it is exactly the lower-slot factorial resultant on the transverse plane.
This gives a suspension theorem.  It proves the full fixed-low terminal tail
through five slots and identifies, without termwise factorial comparison,
the first genuinely new obstruction at six slots.

## 1. Physical family and transverse sidecar

Let `L(s^j)=j!`, `f_j=s^j/j!`, and fix

```text
k>=3, n>=1, 0<=a_0<...<a_(k-3)<C, h>=1, M=C+h.       (1)
```

Use homogeneous coordinates

```text
(x_0,...,x_(k-3),z) in P^(k-2),
w=sum_(i=0)^(k-3)x_i,                                 (2)
```

and put

```text
Y_C=sum_i x_i(f_(n+a_i)-f_(n+M))
       +z(f_(n+C)-f_(n+M)),                           (3)

F_(r,C)=L(Y_C^r)/L(f_n^r),                 2<=r<=k,
R_C=Res(F_(2,C),...,F_(k,C)).                         (4)
```

This is the normalized intrinsic physical resultant, not a formal prime
corner or one selected Macaulay chart.

For `k>=4`, the transverse space is

```text
K={w=z=0} ~= P^(k-4).                                 (5)
```

It has vector dimension `k-3`.  Restrict the fixed layers:

```text
H_r=F_(r,0,C)|_K,                         2<=r<=k-1,
S=Res_K(H_2,...,H_(k-2)).                           (6)
```

All `H_r` are independent of `C`.  The resultant in `(6)` uses `k-3`
forms in `k-3` homogeneous variables; the extra fixed form `H_(k-1)` will
enter the normal contraction but not its obstruction.  The normalization of
`S` uses the determinant-one `y` coordinates fixed in `(16)`.  At the
endpoints use

```text
k=3: K=P^(-1), S=1;
k=4: K=P^0,    S=H_2(K).                              (7)
```

The theorem assumes only

```text
S!=0.                                                  (8)
```

Thus the fixed-low scalar cell is the complete obstruction passed into the
moving tail.  No sign hypothesis on `S` is needed.

## 2. Normal filtration and the exact top pair

Split each physical form by the number of factors drawn from the moving pair
`{f_(n+C),f_(n+M)}`:

```text
F_(r,C)=sum_(q=0)^r F_(r,q,C).                        (9)
```

Every moving factor is linear in `(w,z)`, so exactly

```text
F_(r,q,C) belongs to (w,z)^q.                        (10)
```

For `r=k-1,k`, define

```text
U_(r,C)=L(f_(n+C)^r)/L(f_n^r)
       =(rn+1)_(rC)/(n+1)_C^r,
g_(r,C)=F_(r,C)/U_(r,C).                              (11)
```

The all-high binary normal blocks have full coefficientwise inverse-`C`
expansions and leading forms

```text
F_(r,r,C)/U_(r,C) -> L_r(w,z)^r,
L_r=-r^h w+(1-r^h)z.                                  (12)
```

Put

```text
p=k-1, q=k, Delta=k^h-(k-1)^h>0.                     (13)
```

The two limiting lines are independent.  Define the exact all-high binary
normal blocks and their resultant by

```text
A_(r,C)=F_(r,r,C)/U_(r,C),
E_C=Res_bin(A_(k-1,C),A_(k,C)).                       (14)
```

Equation `(12)` gives the full symbol

```text
E_C=Delta^[k(k-1)](1+O(C^-1))!=0                     (15)
```

for every sufficiently large `C`.  Unlike a boundary factor from an
auxiliary projective flag, `(14)` is intrinsic to the fixed quotient
coordinates `(w,z)`.

## 3. Isotropic normal contraction

For `k>=4`, take the unimodular coordinates

```text
y_i=x_i,                         0<=i<=k-4,
x_(k-3)=w-sum_(i=0)^(k-4)y_i,                           (16)
```

together with `z`.  Thus `K={w=z=0}` with homogeneous coordinates `y`, and
the coordinate determinant is one.  Put

```text
lambda_C=U_(k-1,C)^[-1/(k-1)],
s_C=U_(k-1,C)^[k/(k-1)]/U_(k,C).                      (17)
```

Substitute `(w,z)->lambda_C(w,z)` in every form and then multiply only the
order-`k` form by `s_C`.  If `R~_C` is the resultant of the transformed
system, variable covariance and coefficient-block homogeneity give exactly

```text
mu_(k-1)=k!/(k-1),                 mu_k=(k-1)!,

R~_C=lambda_C^(2k!) s_C^mu_k R_C
    =R_C/[U_(k-1,C)^mu_(k-1) U_(k,C)^mu_k].           (18)
```

This is an isotropic normal scaling.  It does not make the invalid
anisotropic replacement of the two exact normal blocks by their limiting
lines.

The layer estimate behind the contraction is coefficientwise.  For
`2<=r<=k-2` and `1<=j<=r`, or for `r=k-1` and `1<=j<k-1`, a `j`-high layer
after `(17)` costs at most

```text
poly(C)(j/(k-1))^(jC).                                (19a)
```

The omitted `j=0` fixed layer has a separate elementary estimate:

```text
F_(r,0,C)(y,lambda_C w)-H_r(y)
 =O(lambda_C)=O(poly(C)(k-1)^(-C)),       2<=r<=k-1.  (19c)
```

In the scaled order-`k` form, a `j`-high layer with `0<=j<k` costs at most

```text
poly(C)[j^j (k-1)^(k-j)/k^k]^C.                       (19b)
```

Convexity and the endpoint comparison give the common safe gap

```text
rho_k=((k-2)/(k-1))^(k-2)<1,                         (20)
```

which dominates every base in `(19a)--(19c)`; the largest lower base is the
`j=r=k-2` boundary in `(19a)`.  Therefore, coefficientwise on the fixed
projective space,

```text
(F~_2,...,F~_k)
 =(H_2,...,H_(k-2),H_(k-1)+A_(k-1,C),A_(k,C))
   +O(poly(C)rho_k^C).                                (20a)
```

For `k=3`, there are no `y` variables and the transformed pair is exactly
the full normalized pair `(g_(2,C),g_(3,C))`.  The same layer estimates give

```text
(g_(2,C),g_(3,C))
 =(A_(2,C),A_(3,C))+O(poly(C)2^(-C)).                 (21)
```

Thus binary-resultant continuity supplies the same conclusion.  This
whole-form contraction, rather than termwise sorting of equal-base Macaulay
monomials, is the load-bearing separation.

## 4. Separated resultant and the positive face

For fixed `C`, the model system in `(20a)` has resultant

```text
Res(H_2,...,H_(k-2),H_(k-1)+A_(k-1,C),A_(k,C))
 =S^[k(k-1)] E_C^[(k-2)!].                            (22)
```

Indeed `S!=0` makes the homogeneous quotient by
`H_2,...,H_(k-2)` finite of length

```text
D=2*3*...*(k-2)=(k-2)!.                              (23)
```

For every sufficiently large `C`, `(15)` also gives `E_C!=0`.  The lower
equations then force `y=0` at any common cone zero, and the two normal forms
have no common projective zero.  Consequently the model resultant is
nonzero for *every* choice of the coefficients of `H_(k-1)`.  Over `C`, a
nonconstant polynomial in those coefficients has a zero, so the resultant
is independent of `H_(k-1)`.  Set that form to zero.  The standard separated
product formula now repeats `S` through the normal complete intersection of
length `k(k-1)` and repeats `E_C` through the transverse quotient of length
`D`.  Its sign is positive in the standard resultant normalization, as is
seen on the coordinate-power monomial system.  This proves `(22)`; the
exponent identity is

```text
k(k-1)(k-2)!=k!.                                      (24)
```

Since the right side of `(22)` stays away from zero, polynomial continuity
of the fixed-degree resultant applied to `(20a)` and then `(18)` proves

```text
R_C=S^[k(k-1)] E_C^[(k-2)!]
    U_(k-1,C)^[k!/(k-1)] U_(k,C)^[(k-1)!]
    (1+O(poly(C)rho_k^C)).                            (25)
```

Combining `(15)` and `(25)` gives the intrinsic leading face

```text
R_C~S^[k(k-1)] Delta^(k!)
     U_(k-1,C)^[k!/(k-1)] U_(k,C)^[(k-1)!].           (26)
```

Because `k(k-1)` and `k!` are even and every `U` is positive, `(8)` implies

```text
R_C>0 for every sufficiently large C.                 (27)
```

For `k=3`, `(21)` and binary-resultant continuity prove `(25)` with the
empty lower product.  For `k=4`, `(25)` recovers the leading
`G_2^12 E_C^2` mechanism of THM-3062; that theorem's exact Poisson boundary
factor uses its full boundary restriction rather than only `(14)`.

## 5. Exact exponential base and power

Gauss multiplication gives

```text
U_(r,C)=r^(rC) product_(j=1)^(r-1)
         [(n+j/r)_C/(n+1)_C],                         (28)

U_(r,C)~kappa_(r,n) r^(rC) C^[-(r-1)/2],
kappa_(r,n)=Gamma(n+1)^(r-1)/
             product_(j=1)^(r-1)Gamma(n+j/r)>0.       (29)
```

Put

```text
alpha_k=(k-2)! [2(k-1)^2-1]/2.                       (30)
```

Then `(26)` becomes

```text
R_C=A [k(k-1)]^(k! C) C^(-alpha_k)(1+O(C^-1)),       (31)

A=S^[k(k-1)] Delta^(k!)
  kappa_(k-1,n)^[k!/(k-1)]
  kappa_(k,n)^[(k-1)!]>0.                             (32)
```

All normalized coefficients in `(14)` are finite Gamma-ratio expressions,
and the remaining error in `(25)` is exponentially small.  Thus `(31)` has
a complete inverse-power expansion.

## 6. Unconditional five-slot theorem

Take `k=5`, so the fixed lows are arbitrary

```text
0<=a<b<c<C.                                           (33)
```

Now `K=P^1`, and

```text
S_3=Res_K(H_2,H_3).                                   (34)
```

THM-2824 proves that the positive-definite binary quadratic `H_2` and the
real binary cubic `H_3` have no common projective root for every three-slot
factorial support.  In the standard orientation

```text
S_3=leading(H_2)^3 |H_3(alpha)|^2>0.                 (35)
```

Therefore `(8)` is automatic and

```text
R_C=S_3^20 E_C^6 U_(4,C)^30 U_(5,C)^24
    (1+O(poly(C)(27/64)^C)),                          (36)

R_C~S_3^20(5^h-4^h)^120
     kappa_(4,n)^30 kappa_(5,n)^24
     20^(120C) C^(-93)>0.                             (37)
```

This proves the actual normalized five-slot physical terminal resultant is
eventually positive for every fixed triple of lower offsets and fixed gap.

The next case is exact rather than rhetorical.  At `k=6`, `(6)` is the
arbitrary four-slot resultant of degrees `2,3,4` on `P^2`.  THM-2843 gives
its nonnegative norm reduction but not universal nonvanishing.  Hence six
slots are blocked precisely by arbitrary four-slot first-window SFC; no new
terminal-tail cancellation is hiding before that sidecar.

## 7. Exact carrier and affine tail holotopy

Let

```text
A_k=k!(Harmonic(k)-1),
B_k=k!(k+1-2Harmonic(k)),
I_k=A_k+B_k=k!(k-Harmonic(k)).                        (38)
```

THM-3040's exact physical width flag at `t=1/n` is

```text
W_M=(n)_M^A_k (n+1)_M^B_k/n^(I_k M).                 (39)
```

Define the positive exact carrier

```text
B_C=S^[k(k-1)] Delta^(k!)
    U_(k-1,C)^[k!/(k-1)] U_(k,C)^[(k-1)!],
G_C=W_(C+h)B_C,                    Q_C=W_(C+h)R_C.    (40)
```

Formula `(28)` makes `B_C` a scaled product of Beta moments, while `(39)`
is a full-support product-Gamma moment.  Hence `(G_C)` is a strict Stieltjes
product.

Fix integers `s>=1,b_0>=0`, put `C_j=sj+b_0`, and retain the tail on which
all inequalities in `(1)` hold.  For

```text
T_j^(theta)=(1-theta)G_(C_j)+theta Q_(C_j),
0<=theta<=1,                                          (41)
```

every fixed finite generalized-Hankel bank is positive for all sufficiently
large `j`, uniformly in `theta`.  The actual/carrier ratio has a common full
symbol by `(25),(31)`.  The width flag supplies the first confluent
alternant

```text
(I_k s)^N V(u)V(v)/product_(ell=0)^(r-1)ell!>0,
N=r(r-1)/2,                                           (42)
```

and the fixed bank makes the remainder uniform.

For five slots,

```text
(A_5,B_5,I_5)=(154,172,326),                         (43)
```

and the carrier has thirty copies of each
`Beta(n+j/4,1-j/4)`, `j=1,2,3`, plus twenty-four copies of each
`Beta(n+j/5,1-j/5)`, `j=1,2,3,4`.  Its affine alternant order is `326s`.

## 8. Exact evidence and scope

The exact companion checks the normal-ideal filtration, both isotropically
scaled layer-gap families including fixed-layer leakage, the separated
resultant exponent ledger, five-slot transverse binary resultants, exact
all-high normal resultants, Gauss/Beta products, `(154,172,326)` width flag,
and finite carrier Hankel minors.  It also builds the literal degree-11
physical Macaulay map over `F_1000003` on `36` mixed cells.  All maps have
shape `589x364` and full column rank.  This is a finite nonvanishing control,
not an intrinsic resultant value or a tail proof.  No unproved selected
`364x364` determinant chart is used as evidence.

Run

```text
python 04-computation/gmc_terminal_suspension_five_slot_tail_thm3063.py
python -O 04-computation/gmc_terminal_suspension_five_slot_tail_thm3063.py
```

Both modes must equal the stored transcript after LF normalization.

This theorem is fixed-low, fixed-gap, and affine-clock.  It proves eventual
physical resultant positivity and every fixed finite Hankel bank, not one
all-order Stieltjes tail.  It does not handle moving lower offsets, varying
gap, nonlinear clocks, raw-chart persistence, or a wall-stripped core.  It
does not prove arbitrary four-slot SFC, the six-slot tail without `(8)`, or
any new NC2, GMC, LRC, or Jacobian statement.

**QED.**
