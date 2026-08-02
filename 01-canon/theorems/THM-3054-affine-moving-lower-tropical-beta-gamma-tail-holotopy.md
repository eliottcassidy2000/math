---
id: THM-3054
title: "Affine moving-lower tropical Beta-Gamma tail holotopy"
status: >
  PROVED + VERIFIED-EXACT + THREE INDEPENDENTLY HOSTILE-AUDITED.  Along every
  affine penultimate-lower path with fixed terminal
  gap, the unique top tropical face of the formal factorial resultant is an
  explicit strict Beta-Gamma moment sequence.  Every other face is smaller
  by the explicit exponential gap ((k-1)^(k-1)/k^k)^C.  Conditional on the
  frozen lower-face resultant, every fixed finite bank of generalized Hankel
  minors is eventually positive, uniformly along the straight
  carrier-to-resultant homotopy.  The sidecar is automatic for k<=4 and at
  every sufficiently large depth for fixed general k.  Quantifiers do not
  commute: no one infinite Stieltjes tail or physical-width claim follows.
source: kind-pasteur-2026-08-01-affine-moving-lower-tail
audit: >
  A third independent hostile audit ACCEPTED the multihomogeneous top-face
  selection, Poisson specialization, Beta-Gamma factorization, affine
  pushforward, fixed-bank alternant asymptotics, perturbative homotopy,
  k<=4 sidecar, and all scope hostiles.  Fresh ordinary and optimized runs
  matched the stored eleven-line LF transcript and both declared hashes.
  The reconciliation clarifies that rho_k is the relaxed redistribution bound,
  attained by the exact three-slot correction, rather than a claimed census
  of every resultant Newton monomial.
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-3040-formal-corner-resultant-width-quotient-and-all-order-bernoulli-law
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3051-stieltjes-multiplier-gamma-flow-and-moving-lower-hankel-boundary
  - THM-3053-beta-gamma-prefix-transport-and-multiplicative-holotopy-cone
related:
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
script: 04-computation/gmc_affine_moving_lower_tail_holotopy_thm3054.py
output: 05-knowledge/results/gmc_affine_moving_lower_tail_holotopy_thm3054.out
script_sha256: 9a6ddbf51b7eb35554b492bda7e466ff4c7e4b4de37571d71ac0d294e6b60a3b
output_sha256: ea06c3ab193be66d6c30f557e1c34597b01f890a3d47778e0257c9c52c3c12a7
hash_basis: LF-normalized bytes
---

# THM-3054 -- affine moving-lower corners return to every finite Hankel chamber

**PROVED + VERIFIED-EXACT + THREE INDEPENDENTLY HOSTILE-AUDITED.**

THM-3051 proves that literal moving lower offsets can leave the Stieltjes
cone at the first nontrivial widths.  THM-3053 identifies the right positive
carriers: ordered Beta-Gamma transport, not coordinatewise Gamma inventory.
The two facts fit together through the tropical Newton face of the resultant.

The top moving face is exactly Beta-Gamma.  Every other face is exponentially
smaller.  Therefore the full moving corner need not be Stieltjes at the
beginning, but it eventually re-enters every **fixed finite** total-positivity
chamber.  The straight path from the carrier to the full corner stays in that
same finite chamber after an order-dependent threshold.  This is the tail
holotopy promised by the early hostiles.

## 1. Affine moving-lower system

Fix:

```text
slot number k>=2,          integer depth n>=1,
fixed lower offsets 0<=d_0<...<d_(k-3),
affine moving offset C_c=s c+b, integers s>=1, b>=0,
fixed terminal gap h>=1,  M_c=C_c+h.                  (1)
```

For `k=2`, the fixed tuple is empty.  We only consider tail indices for which
the displayed offsets are distinct and ordered.

Take the coefficientwise formal terminal corner of THM-3040.  Its surviving
factorial forms have degrees `r=2,...,k` in the `k-1` lower variables.  Denote
them by `G_(r,C)` and their intrinsic multivariate resultant by

```text
R_C=Res(G_(2,C),...,G_(k,C)).                          (2)
```

Set the moving variable to zero in the first `k-2` forms and write

```text
S=Res(G_2^0,...,G_(k-1)^0),                           (3)
```

with the conventions `S=1` at `k=2` and the evident one-form resultant at
`k=3`.  This is the **frozen lower-face sidecar**.  For now assume `S!=0`.

The coefficient of the pure moving monomial in the degree-`k` form is

```text
U_C=(kn+1)_(kC)/(n+1)_C^k,                            (4)
```

and put

```text
mu=(k-1)!.
```

## 2. Unique tropical face and sharp gap

Every resultant monomial has total moving-coordinate weight

```text
product_(r=2)^k r=k!.                                 (5)
```

Indeed, `(5)` is the character of scaling the moving variable.  A coefficient
containing `j>=1` copies of that variable has, for fixed other offsets,

```text
positive constant * j^(jC) C^beta (1+O(1/C));         (6)
```

coefficients with `j=0` are fixed.  Formula `(6)` follows directly from
Stirling's Gamma-ratio expansion.

Because `j log j` per unit weight is strictly maximized at `j=k`, the unique
top face uses `mu` copies of the pure degree-`k` coefficient `(4)` and only
moving-free coefficients in every lower form.  Poisson specialization of the
resultant gives that face exactly:

```text
sigma S^k U_C^mu, sigma in {+1,-1}.                   (7)
```

Here `sigma` is one fixed resultant-orientation sign.

There is a universal gap with a partition-sharp relaxed constant.  The
relaxed redistribution bound is attained for every `k>=3`.  If one pure
coefficient is removed, its weight `k` is best redistributed as `(k-1)+1`,
so the ratio of exponential bases is

```text
rho_k=(k-1)^(k-1)/k^k<1.                              (8)
```

Removing more copies does no better.  Ignore the additional constraints of
the actual resultant Newton support and write `L` for the number removed.
For `1<=L<=k-1`, convex integer redistribution gives the upper bound

```text
max ratio=(rho_k L)^L<=rho_k,                         (9)
```

with equality in the relaxed partition only at `L=1`; increasing `L` by
`k-1` multiplies the ratio by `((k-1)/k)^(k(k-1))`.  Thus `(9)` bounds every
actual face for every `L>=1`.  The exact three-slot resultant attains the
constant, where `rho_3=4/27`.

The resultant has finite Newton support.  Fixed-shift Gamma-ratio bounds,
`(6)--(9)`, and `S!=0` therefore give effective constants `D,K>0` such that,
after orienting once by `epsilon=sign(sigma S^k)`,

```text
epsilon R_C=|S|^k U_C^mu(1+eta_C),
|eta_C|<=K C^D rho_k^C                               (10)
```

for every sufficiently large integer `C`.  Along `(1)`, the error is
`O(c^D rho_k^(s c))`.

## 3. The top face is exactly Beta-Gamma

Let THM-3047's character be

```text
A=k!(H_k-1), B=k!(k+1-2H_k), I=A+B.                  (11)
```

For the unit moving clock, define the positive carrier

```text
g_C=|S|^k F_(C+h)^(k)(1/n) U_C^mu.                   (12)
```

Gauss multiplication and THM-3047 give the exact factorization

```text
g_C=g_0 Lambda^C
 (n+h)_C^A (n+h+1)_C^B
 product_(q=1)^(k-1)
 [(n+q/k)_C/(n+1)_C]^mu,                             (13)

Lambda=k^(k!)/n^I.                                   (14)
```

Thus `(g_C)` is the moment sequence of a scaled product of

```text
A copies of Gamma(n+h,1),
B copies of Gamma(n+h+1,1),
mu copies of Beta(n+q/k,1-q/k) for every q=1,...,k-1. (15)
```

It has full positive support and is strictly generalized-Hankel totally
positive by the product-support/Andreief argument used in THM-3053; that
strictness argument does not require its shape parameters to lie on the
integer lattice.

The affine subsequence in `(1)` remains strict Stieltjes without any new
calculation:

```text
g_(s c+b)=integral x^(s c+b)dmu(x)
          =integral y^c dnu(y), y=x^s,               (16)
```

where `nu` is the pushforward of `x^b mu(dx)`.  Full support survives.  An
`s`-fold multiplication formula identifies `(16)` explicitly as a product of
`sI` Gamma factors and the corresponding subdivided Beta ratios; the
companion checks both descriptions exactly.

## 4. Uniform finite-order tail holotopy

Combine `(10)` and `(12)` along the affine clock:

```text
H_c=epsilon F_(C_c+h) R_(C_c)
    =g_(C_c)(1+eta_(C_c)).                            (17)
```

For `0<=theta<=1`, define the straight coefficient homotopy

```text
H_c^(theta)=g_(C_c)(1+theta eta_(C_c)).               (18)
```

Fix an order bound `d` and offset bound `B_0`.  Then there is an effective
tail threshold `c_0=c_0(k,n,d_i,s,b,h,d,B_0)` such that, for every

```text
c>=c_0, 0<=theta<=1, 1<=r<=d,
0<=rho_1<...<rho_r<=B_0,
0<=tau_1<...<tau_r<=B_0,                              (19)
```

one has

```text
det[H_(c+rho_i+tau_j)^(theta)]_(i,j=1)^r>0.           (20)
```

The proof exposes the quantifiers.  For any carrier

```text
q_m=C_0 L^m product Gamma/Beta Pochhammer ratios       (21)
```

with `P>=1` net Gamma factors, fixed `rho,tau`, and
`N_r=r(r-1)/2`, factor rows and columns in its generalized minor.  The
remaining functions satisfy

```text
Q_t(x)=q_(x+t)/(L^t q_x)=x^(Pt)(1+O(1/x)).            (22)
```

Taylor expansion at `x=c` shows that the first nonzero alternant uses
derivative orders `0,...,r-1`.  Its exact leading term is

```text
[P^N_r V(rho)V(tau)/product_(ell=0)^(r-1)ell!]
 c^(P sum_j tau_j-N_r).                               (23)
```

It is positive.  For the affine carrier, `P=sI`.  Multilinearity and `(10)`
bound the relative determinant error, uniformly in `theta`, by

```text
O(c^(D+N_r) rho_k^(s c))=o(1).                       (24)
```

The bank `(19)` is finite, so one threshold works for all of it.  Equations
`(23)--(24)` prove `(20)`.

This is a literal path inside every chosen finite strict-total-positivity
chamber.  Consequently each corresponding finite Hankel transform is
variation diminishing, and the associated Stieltjes orthogonal polynomials
through any fixed degree retain positive simple roots and strict interlacing
throughout the tail homotopy.  Both consequences retain the order-dependent
threshold.

## 5. The frozen sidecar is automatic through four slots

The hypothesis `S!=0` is not an extra debt at small slot number.

- At `k=2`, it is the empty convention.
- At `k=3`, `S` is one positive monomial coefficient.
- At `k=4`, `S` is a binary quadratic/cubic factorial resultant, and the
  following exact lemma makes it strictly positive.

Normalize two lower offsets to `(0,c)` at depth `n`.  Put

```text
d=(n+1)_c,
x=(2n+1)_c/d,       y=(2n+1)_(2c)/d^2,
p=(3n+1)_c/d,       q=(3n+1)_(2c)/d^2,
r=(3n+1)_(3c)/d^3.                                  (25)
```

The forms are

```text
G_2=1+2xz+yz^2,       G_3=1+3pz+3qz^2+rz^3.          (26)
```

Define

```text
C_0=y^2+2xr-3yq,
C_1=4x^2r-6xyq+3y^2p-yr.                             (27)
```

Exact Euclidean division gives

```text
Res(G_2,G_3)
 =[(C_1-xC_0)^2+(y-x^2)C_0^2]/y^2.                  (28)
```

Now `y>x^2`.  Also AM--GM gives

```text
y^2+2xr>=3(xyr)^(2/3)>3yq,                           (29)
```

where the strict second inequality is equivalent to `yq^3<x^2r^2`.
After clearing Pochhammers, that inequality follows termwise from

```text
(2n+i)(3n+2c+i)^2
 -(2n+c+i)(3n+i)(3n+c+i)
 =c(3ci+5cn+2i^2+9in+9n^2)>0                        (30)
```

for `i=1,...,c`.  Hence `C_0>0`, and `(28)` is strictly positive.

Any pair `0<=a<c` reduces to `(0,c-a)` at depth `n+a`, multiplied by positive
form scalars.  Thus the lemma covers every fixed lower pair and proves the
`k=4` assertion.

For arbitrary fixed `k`, `S(t)` is a rational analytic low resultant and
THM-3040 proves `S(0)!=0` by its generalized-Vandermonde unit.  Therefore

```text
S(1/n)!=0 for every sufficiently large integer depth n.                (31)
```

At `k>=5`, a finite set of small depths may still hit `S=0`; there `(7)`
cancels and the next tropical face is genuinely new.  No all-depth claim is
made without checking that sidecar.

## 6. Exact three-slot model and sharp correction

For the three-slot support `(0,c,c+1)`, the variables `(25)` give the generic
thirteen-term resultant

```text
R_c=r^2-6xqr-6ypr+6xyr-8x^3r+12x^2pr
    +12x^2yq-6xy^2p-18xypq+y^3+9y^2p^2-6y^2q+9yq^2. (32)
```

The exact leading correction is

```text
R_c/r_c^2
 =1-C_n c^(2n+1/2)(4/27)^c(1+O(1/c)),                (33)
```

where

```text
C_n=6 Gamma(n+1/3)Gamma(n+2/3)
 [ Gamma(n+1)/(Gamma(2n+1)Gamma((3n+1)/2)Gamma((3n+2)/2))
   +1/(Gamma(n+1/2)Gamma(3n+1)) ]>0.                 (34)
```

In particular

```text
C_1=8 sqrt(3 pi)/9,       C_2=40 sqrt(3 pi)/2187.    (35)
```

The two leading error terms in `(32)` are `-6xqr` and `-6ypr`; `+6xyr` is
smaller by `c^(-n)`, and every other exponential base is below `4/27`.
Formula `(28)` strengthens the size-one conclusion: `R_c>0` for every
`n,c>=1`, even before the eventual-minor theorem begins.

The correction itself is not a hidden Stieltjes multiplier.  At `n=1`, if
`E_c=R_c/r_c^2`, then

```text
(E_1,E_2,E_3)=(1/25,53/140,28621/38115),
E_1E_3-E_2^2=-12089453/106722000<0.                  (36)
```

Nor does inversion repair it: the exact nineteenth forward difference of
`1/E_c` at `c=5` has the wrong positive sign for a Hausdorff moment sequence.
Thus `(17)--(24)` genuinely use exponentially small perturbation of the
whole carrier; they do not smuggle in an exact multiplier theorem.

## 7. Sharp scope hostiles

Both affine spacing and fixed gap are load-bearing even for the dominant
carrier.  At `k=3,n=1`, write

```text
g(C,h)=F_(C+h) U_C^2.                                  (37)
```

The nonaffine clock `C=(1,3,4)` gives

```text
g(1,1)g(4,1)-g(3,1)^2
 =-14623245661460508379818491904000000000<0.         (38)
```

Keeping `C=(1,2,3)` but varying the gaps as `(1,2,1)` gives

```text
g(1,1)g(3,1)-g(2,2)^2
 =-80810821433972705093222400000000<0.               (39)
```

Composition with a general increasing index path does not preserve moment
sequences; the affine clock works because it is the single pushforward
`x->x^s` in `(16)`.

Most importantly, the theorem has quantifier order

```text
for every fixed (d,B_0), there exists c_0(d,B_0).      (40)
```

It does **not** say that one tail is Stieltjes at all orders.  The infinite
Hankel cone has no uniform projective interior supplied by this argument.
THM-3051's early negative minors remain exact positive controls for this
boundary.

Everything here concerns the intrinsic coefficientwise formal corner.  It
does not prove positivity of a physical width, raw selected Macaulay chart,
wall-stripped core, global jet polynomial, or arbitrary nonlinear lower path.
It proves no new NC2, GMC, LRC, Jacobian, tournament, or integer-sequence
statement.

## 8. Exact companion

The dependency-free exact referee performs:

- `243` independent specializations of the generic thirteen-term resultant;
- `320` all-positive lower resultants, `3,360` termwise factor identities,
  and `1,050` arbitrary-pair translation identities;
- `120` unit-clock and `864` affine-subsequence Beta-Gamma carrier identities;
- `45` relaxed tropical-gap optimization cells for `k=2..10`;
- six exact polynomial-alternant leading-coefficient identities `(23)`;
- `648` straight-homotopy finite-window recovery controls;
- the correction, inverse-Hausdorff, nonaffine-clock, and varying-gap hostiles
  `(36),(38),(39)`.

Run

```text
python 04-computation/gmc_affine_moving_lower_tail_holotopy_thm3054.py
python -O 04-computation/gmc_affine_moving_lower_tail_holotopy_thm3054.py
```

Both modes equal the stored eleven-line transcript after LF normalization.

Two independent immutable-file audits rederived the all-`k` tropical
optimization, the exact Beta--Gamma factorization, the rational-ratio
alternant, the `k<=4` sidecar, and the three-slot correction by paths
independent of the companion.  Both replayed normal and optimized execution
against the stored transcript and found no defect.

**QED.**
