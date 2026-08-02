---
id: THM-3051
title: "Stieltjes multiplier, Gamma flow, and moving-lower Hankel boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Hadamard
  multiplication by a Stieltjes sequence preserves every
  Stieltjes sequence, and this characterizes universal Stieltjes multipliers.
  For THM-3047's fixed-lower product-Gamma flag, adjacent Gamma-shape flow gives
  a larger explicit preservation cone.  Neither adjacent log-convexity nor the
  universal flag controls literal moving lower offsets: a three-slot family
  fails the adjacent Hankel inequality at every integer depth n>=2, and an
  exact four-slot factorial-resultant corner fails already at n=2.  These are
  formal-corner transport statements, not physical-width or raw-chart claims.
source: kind-pasteur-2026-08-01-moving-lower-moment-boundary
audit: >
  Two independent immutable-file audits ACCEPTED the universal multiplier
  iff, strictness boundary, Gamma inventory and both escape examples, the
  order-three curvature identity, all-depth binary-resultant hostile, and the
  four-slot Macaulay hostile.  Each independently reconstructed the moving
  resultants by a path different from the companion, replayed normal and
  optimized execution against the stored transcript, matched both LF hashes,
  and passed the documentation checker.  Their one shared scope correction --
  examples disprove universal preservation but not every moving path -- is
  incorporated below.
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-3040-formal-corner-resultant-width-quotient-and-all-order-bernoulli-law
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
related:
  - THM-3021-the-hadamard-multiplier-question-is-false-and-sfc2-is-appell-squarefreeness
  - THM-3050-rational-product-gamma-radial-nullcone-and-critical-borel-order
script: 04-computation/gmc_stieltjes_gamma_flow_moving_lower_thm3051.py
output: 05-knowledge/results/gmc_stieltjes_gamma_flow_moving_lower_thm3051.out
script_sha256: 8a37ce6412af4292d400854604f7faf7927fc4db1a6e57fc92c5c96fa6e7d704
output_sha256: 8e1d86f9d66644dcec1f5983b7db2e6e03bc0fca7c28449d1d40a029f5f13b57
hash_basis: LF-normalized bytes
---

# THM-3051 -- the exact preservation cone stops before moving lower offsets

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3047 proves that the universal fixed-lower width flag is a strict
Stieltjes moment sequence.  The natural inheritance question is whether a
lower-resultant transport factor preserves that moment law.  There are two
positive mechanisms: multiplicative convolution with another Stieltjes law,
and literal transfer among adjacent Gamma shapes.  They are genuine, but they
do not cover arbitrary moving lower offsets.

The failure is already visible before any global physical specialization.  A
three-slot moving family crosses the first Hankel wall at every depth `n>=2`,
and the actual `(2,3,4)` factorial system does the same at depth `2`.  Thus the
remaining transport problem is not to strengthen adjacent curvature.  It is
to find extra compatibility carried by a physical lower-offset path.

## 1. Universal Stieltjes multipliers

Let

```text
F_M=integral_(0,infinity) x^M dmu(x)                 (1)
```

be THM-3047's product-Gamma flag at some fixed `k>=2,t>0`.  Its representing
measure `mu` charges every nonempty open interval in `(0,infinity)`.  Let

```text
L_M=integral_[0,infinity) y^M dlambda(y)              (2)
```

for a finite positive measure with all moments finite.  Then

```text
F_M L_M=integral (xy)^M d(mu tensor lambda)(x,y).     (3)
```

Hence `(F_M L_M)` is Stieltjes: its measure is the pushforward of
`mu tensor lambda` by multiplication.  It is strictly generalized-Hankel
totally positive exactly when

```text
lambda((0,infinity))>0, equivalently L_1>0.           (4)
```

Indeed, under `(4)` the product measure charges every open interval in
`(0,infinity)`, so the generalized-Vandermonde/Andreief proof of THM-3047
applies verbatim.  If `lambda` is supported at zero, every positive moment
vanishes and strictness fails.  Notice that `lambda` need not have total mass
one; `(3)` is multiplicative convolution of finite measures, not an implicit
probability normalization.

More generally, for a real sequence `L=(L_M)_(M>=0)`,

```text
L Hadamard-multiplies every Stieltjes sequence to a Stieltjes sequence
  iff L itself is a Stieltjes sequence.                (5)
```

Sufficiency is `(3)`.  For necessity, apply the presumed multiplier to the
constant sequence `1`, represented by the point mass at `1`.

## 2. Adjacent Gamma-shape flow

Write `a=1/t`, and retain THM-3047's inventory

```text
e_0=A, e_1=B, e_j=0 (j>=2), I=A+B.                   (6)
```

For `j>=0`, define the normalized adjacent transfer

```text
T_j(M)=(a+j+M)/(a+j)=(a+j+1)_M/(a+j)_M.              (7)
```

Let `d_j` be finitely supported integers, put `d_(-1)=0`, and let `c>0`.
Set

```text
L_M=c^M product_(j>=0) T_j(M)^(d_j),
n_j=e_j+d_(j-1)-d_j.                                  (8)
```

Direct cancellation gives the exact identity

```text
F_M L_M=(c t^I)^M product_(j>=0)(a+j)_M^(n_j).        (9)
```

Consequently, if every `n_j>=0`, `(F_M L_M)` is the moment sequence of a
scaled product of `n_j` independent `Gamma(a+j,1)` variables.  Since
`sum_j n_j=I>0`, the product law has full support and is strictly
generalized-Hankel totally positive.

Condition `n_j>=0` is an iff only for the displayed **cancellation-free
adjacent Gamma inventory**.  It is not necessary for arbitrary Stieltjes
representability.  The sharp escape is already at `k=2`.  With `d_0=-1`,

```text
F_M/T_0(M)=t^M (a)_M^2/(a+1)_M
           =t^M (a)_M a/(a+M),                       (10)
```

which has `n_1=-1` but is the moment sequence of `t X B`, where
`X~Gamma(a,1)` and `B~Beta(a,1)` are independent.

Nor must a multiplier preserving this one fixed `F` itself be Stieltjes.
For `d_0=1`,

```text
L_M=T_0(M)=1+Mt,
det[[L_0,L_1],[L_1,L_2]]=-t^2<0,                     (11)
```

but `F_M L_M` merely transfers one of the `A>=1` factors from shape `a` to
shape `a+1`, so it remains strict product-Gamma.  Thus `(5)` is a universal
preserver theorem; the multiplier cone of one fixed interior point is larger.

## 3. The first higher-Hankel gate

For any positive sequence `g`, put

```text
q_M=g_(M+1)/g_M,
kappa_j=g_(j-1)g_(j+1)/g_j^2.                        (12)
```

Writing

```text
u=kappa_(M+1), v=kappa_(M+2), w=kappa_(M+3),
Psi(u,v,w)=u v^2 w-u v^2-v^2 w+2v-1,                (13)
```

direct expansion gives

```text
det[g_(M+i+j)]_(0<=i,j<=2)
 =g_M^3 q_M^6 u^3 Psi(u,v,w).                        (14)
```

Hadamard products multiply adjacent curvatures:

```text
kappa_M(F L)=kappa_M(F) kappa_M(L).                  (15)
```

But the three-variable compatibility `Psi>=0` does not follow from the three
separate inequalities `u,v,w>1`.  The exact hostile

```text
g_M=M!(M+1)^2                                         (16)
```

is strictly log-convex, because for `M>=1`

```text
kappa_M=M(M+2)^2/(M+1)^3>1,
M(M+2)^2-(M+1)^3=M^2+M-1>0.                         (17)
```

Nevertheless

```text
det[g_(i+j)]_(0<=i,j<=2)=-24.                        (18)
```

Thus even strict adjacent curvature is only the first wall of the Stieltjes
cone.  Formula `(14)` is the cheapest exact order-three compatibility test.

## 4. Three slots: failure at every depth `n>=2`

Take the literal moving family of three-slot formal corners with support

```text
(0,c,c+1), c=1,2,3, and t=1/n.                       (19)
```

After the terminal face is removed, the normalized binary lower forms are

```text
G_(r,c)(X,Y)=sum_(j=0)^r binom(r,j)
  (rn+1)_(jc)/(n+1)_c^j X^(r-j)Y^j, r=2,3.          (20)
```

Let `R_c=Res(G_(2,c),G_(3,c))`, with the standard descending-coefficient
Sylvester convention.  Exact elimination gives

```text
R_1=A(n)/(n+1)^4,
R_2=B(n)/((n+1)^4(n+2)^4),
R_3=C(n)/((n+1)^4(n+2)^5(n+3)^4),                   (21)
```

where

```text
A=n^4+22n^3+69n^2+44n+8,

B=15625n^8+427500n^7+2850950n^6+8097488n^5
  +11893269n^4+9683740n^3+4378980n^2+1022112n+95616,

C=47045881n^13+1797908516n^12+23919687803n^11
  +167453540986n^10+722042412723n^9+2070797182104n^8
  +4110657369269n^7+5751069261826n^6+5683612606852n^5
  +3918712637960n^4+1829514912288n^3+546467611968n^2
  +93389808384n+6890365440.                          (22)
```

The `t=0` control is

```text
R_c -> (3^c-2^c)^6>0.                                (23)
```

For `k=3`, `(A_3,B_3,I_3)=(5,2,7)`.  Put

```text
H_c=F_(c+1)^(3)(1/n) R_c.                            (24)
```

Then the adjacent curvature is exactly

```text
K(n)=H_1H_3/H_2^2
 =(n+4)^2 A(n)C(n)/((n+2)^2(n+3)B(n)^2).             (25)
```

Let the denominator of `(25)` minus its numerator be `P(n)`.  The exact
certificate is

```text
P(m+2)=sum_(j=0)^19 p_j m^j,

(p_0,...,p_19)=
(3268881760112640000,37643475128286720000,
 146291966980162617344,311943061093012269312,
 433931551733794554880,429003646190355919616,
 316601230236031112592,179745211720414744752,
 79988104378771836568,28208693320266071200,
 7923650163754635817,1772672109839816584,
 314264085331109200,43683300848187864,4678174175000898,
 375655828575392,21687531118480,840873878616,19348674701,
 197094744).                                         (26)
```

Every coefficient in `(26)` is positive.  Hence `K(n)<1` for every integer
`n>=2`.  The boundary is sharp:

```text
K(1)=3577625/2123604>1.                               (27)
```

At the first failing depth,

```text
(R_1,R_2,R_3)=(188/27,1025675/27,44055606),
H_1H_3-H_2^2=-9090480053413828125/512.               (28)
```

Thus literal lower motion destroys even the adjacent Hankel inequality at
every depth `n>=2`, despite strict total positivity of the universal factor.

## 5. Four slots: an exact GMC-shaped hostile

The same failure occurs in the actual `(2,3,4)` factorial system, not only in
the binary toy.  Fix `n=2` and take the moving supports

```text
(0,1,c,c+1), c=2,3,4.                                (29)
```

For `r=2,3,4`, the normalized ternary low form is

```text
G_(r,c)(x_0,x_1,x_2)
 =sum_(alpha_0+alpha_1+alpha_2=r) multinomial(r;alpha)
   (2r+1)_(alpha_1+c alpha_2)
   /(3^alpha_1 (3)_c^alpha_2)
   x_0^alpha_0 x_1^alpha_1 x_2^alpha_2.              (30)
```

Write `R_c=Res(G_(2,c),G_(3,c),G_(4,c))`.  A direct degree-seven Macaulay
determinant, with its nonzero extraneous flag divided exactly, gives

```text
R_2=671696427641384054000000000/1162261467,

R_3=4720147255045226732121521309597418358046720/19683,

R_4=23877237441067576157124399040421454340333007132221361662211000369152
    /2017815046875.                                   (31)
```

For reproducibility, order degree-`d` monomials by increasing first and then
second exponent.  Use all `46` rows `G_2*Mon_5,G_3*Mon_4,G_4*Mon_3` and select
global row indices

```text
0..19, 21..29, 35, 36..41.                           (32)
```

The selected `36x36` determinant `Delta` satisfies

```text
R_c=Delta/(q_200^6 c_300 K_flag),
K_flag=c_120 q_200^2-c_210 q_110 q_200
       -c_300 q_020 q_200+c_300 q_110^2=28/9.        (33)
```

Here `q_200=c_300=1`, so every division gate is explicit and nonzero.  At
`t=0` the response determinant for `c=2` is `2`, and the resultant is the
positive control `2^24`.

Now `(A_4,B_4,I_4)=(26,20,46)` and

```text
H_c=F_(c+1)^(4)(1/2)R_c.                             (34)
```

The universal flag contributes `(6/5)^26(7/6)^20` at the central width
`M=4`.  Exact reduction yields

```text
H_2H_4/H_3^2
=44087099546187941338084870318252428808920867535835042464970681136337702301141312674129
 /148395225288083619662436115015878402556351144600145942926978139965896606445312500000000
<1.                                                    (35)
```

So the four-slot moving-lower formal corner already lies outside the
Stieltjes cone at depth `2`.

## 6. Exact frontier and loss ledger

The source-to-target map is now precise:

```text
lower transport factor L
  --Hadamard multiplication--> full formal-corner width flag F*L.          (36)
```

It preserves Stieltjes positivity when `L` has a positive multiplicative
measure, and it also preserves the explicit non-Stieltjes Gamma-flow cone of
Section 2.  Arbitrary literal lower motion need not preserve this predicate,
as Sections 4--5 show.  The quotient to adjacent curvatures loses the
compatibility polynomial `Psi`; even retaining `Psi` at one index would lose
all larger Hankel minors.

Consequently, a positive moving-lower theorem needs a sidecar such as a
single multiplicative-convolution coupling, a totally-positive Gamma/Beta
kernel, or another path-level compatibility tying successive low
resultants together.  Scalar positivity of each `R_c`, pointwise positivity
of every adjacent factor, and strict log-convexity are insufficient.

This theorem does not weaken THM-3050 or NC2.  Their moment index is the power
of one fixed radial law.  Here the index changes the support itself.  Nor do
`(19)--(35)` assert a negative physical width, a nonzero raw selected chart,
or a failure of GMC: they are coefficientwise formal-corner transport
hostiles.

## 7. Exact companion

The dependency-free companion performs:

- `108` finite-atomic multiplicative-convolution and Hankel-PSD checks;
- `81` adjacent Gamma-flow and `27` Beta-escape identities;
- `12` independent order-three curvature determinant identities and the
  exact `-24` strict-log-convex hostile;
- `96` direct binary Sylvester-resultant checks, the complete positive
  shifted polynomial certificate `(26)`, its frozen checksum, and `(28)`;
- all three `36x36` fraction-free Macaulay determinants in `(31)--(33)`, the
  `2^24` control, and the exact curvature `(35)`.

Run

```text
python 04-computation/gmc_stieltjes_gamma_flow_moving_lower_thm3051.py
python -O 04-computation/gmc_stieltjes_gamma_flow_moving_lower_thm3051.py
```

Both modes equal the stored nine-line transcript after LF normalization.

**QED.**
