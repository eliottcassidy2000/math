---
id: THM-2564
title: "Six-tooth doubly centered tomography and the section-holonomy boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  On the rational doubly centered C_7-by-F_13 interaction space, the
  toothpush maps in any prescribed six nonzero slopes form an isomorphism
  onto six zero-sum root profiles. Every nonzero interaction survives in at
  least seven of twelve slopes, sharply: every prescribed five-slope zero
  set is realized by an integral literal field-trace hostile. Each surviving
  rational profile fires all twelve nontrivial C_13 characters. Integrally,
  every six-chart bank has determinant 13^21 and Smith form
  1^51 direct-sum 13^21. The charts depend on an ordered C_7 section; its
  wrap carry is nonconstant and prevents a canonical physical root action.
  The outputs remain signed and can realize the exact six-replica pattern
  demanded by THM-2562, so no Boolean ancestry coupling, scalar-row
  exclusion, or LRC(14) conclusion is claimed.
source: root-holotopy-2026-07-27-six-tooth-tomography
depends_on:
  - THM-2507-truncated-radon-toothpick-tomography-and-nonaffine-root-boundary
  - THM-2557-doubly-centered-transfer-moment-index-and-hall-cone-boundary
related:
  - THM-2559-target-informed-chord-and-universal-old-repair-packet
  - THM-2561-positive-interval-physical-six-comb-blind-necklace-hostile
  - THM-2562-canonical-duty-commutator-line-rank-and-anchor-rigidity
  - THM-2563-paired-dipole-deep-target-corner-and-partial-bare-boundary
script: 04-computation/lrc14_six_tooth_v00_tomography_thm2564.py
output: 05-knowledge/results/lrc14_six_tooth_v00_tomography_thm2564.out
script_sha256: ab74dfe96998ab2f5efef0b96562c5b87a504a1b6e2d5ab44a2f7b6bcb0bd4b8
output_sha256: b17d69e9c9183b0e2a1394f5bc4fc6d4e2a5cde532fc80e83a55336d3da6fc29
hash_basis: working-tree bytes (LF)
---

# THM-2564 -- six teeth reconstruct the signed interaction

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT
PENDING.**

THM-2507 proves truncated Radon tomography on a row-zero strip, modulo a
vertical kernel.  THM-2557 places the live typed interaction in the smaller
doubly centered lattice.  Transposing the strip and imposing the second
margin kills the entire vertical kernel.  The resulting six-chart atlas is
therefore exact, and its integral defect can also be computed.

This is a reconstruction theorem for a **signed interaction**.  It is not a
physicalization theorem.  The section needed to draw each tooth has a real
wrap cocycle, and six-chart surjectivity is compatible with the precise
six-replica cancellation pattern isolated by THM-2562.

## 1. General doubly centered tooth theorem

Let `p` be prime and

```text
2 <= q < p.
```

Use the ordered set `I_q={0,...,q-1}` and define

```text
V_00(p,q)
 ={f:I_q x F_p -> Q:
     sum_v f(k,v)=0 for every k,
     sum_k f(k,v)=0 for every v}.                              (1)
```

For `tau in F_p`, put

```text
(T_tau f)(v)=sum_(k=0)^(q-1) f(k,v-tau k).                     (2)
```

Every output lies in the zero-sum space `V_0(F_p)`.  For any set of `m`
distinct nonzero slopes,

```text
S subset F_p^*,

rank_Q T_S=(p-1) min(m,q-1),                                  (3)

T_S f=(T_tau f)_(tau in S).
```

In particular, every prescribed `q-1` slopes give an isomorphism

```text
T_S:V_00(p,q)  ~=  direct_sum_(tau in S) V_0(F_p).             (4)
```

Consequently every nonzero `f` survives in at least

```text
p-q+1                                                        (5)
```

of the `p-1` nonzero slopes.  This is sharp over the integral lattice.

### Fourier proof

Fix a primitive `p`-th root `zeta`.  For `alpha in F_p`, write

```text
F_alpha(k)=sum_v f(k,v) zeta^(-alpha v),

P_alpha(X)=sum_(k=0)^(q-1) F_alpha(k) X^k.                    (6)
```

The first margin in (1) says `F_0(k)=0`; the second says

```text
P_alpha(1)=0.                                                 (7)
```

Changing variables `v=r+tau k` in (2) gives

```text
widehat(T_tau f)(alpha)=P_alpha(zeta^(-alpha tau)).            (8)
```

For each `alpha!=0`, the allowed slice is exactly

```text
P_alpha(X)=(X-1)Q_alpha(X),       deg Q_alpha <=q-2,           (9)
```

a `(q-1)`-dimensional space.  The points

```text
zeta^(-alpha tau),                 tau in S,
```

are distinct and avoid one.  After dividing by `X-1`, (8) is ordinary
Vandermonde evaluation, of rank `min(m,q-1)`.  There are `p-1` nontrivial
Fourier slices, proving (3).  The rank is unchanged after extending scalars
from `Q` to `Q(zeta)`.

If a nonzero interaction vanished in `q-1` slopes, (4) would kill it.
Thus at most `q-2` slopes are bad, proving (5).

## 2. The integral atlas has one uniform prime invoice

Let `V_00(p,q;Z)` and `V_0(F_p;Z)` denote the corresponding integral
lattices.  For every `(q-1)`-slope bank, (4) restricts to an injective
integral map with

```text
e=q(q-1)/2,

|det T_S|=p^e,

Smith(T_S)=1^((q-1)(p-1)-e) direct-sum p^e.                   (10)
```

Here the exponent on `p` in the Smith display is its multiplicity: every
nontrivial invariant factor is exactly `p`, not a higher power.

### Determinant

Choose the coefficient basis of polynomials of degree at most `q-1`
vanishing at one.  On a fixed nonzero character, the evaluation determinant
is, up to a unit,

```text
product_(tau in S)(x_tau-1)
 product_(tau<sigma)(x_sigma-x_tau),

x_tau=zeta^(-alpha tau).                                     (11)
```

The Fourier changes of basis occur `q-1` times on both source and target
and cancel.  Taking the product of (11) over `alpha=1,...,p-1`, every anchor
difference and every pair difference has cyclotomic norm `p`.  There are

```text
(q-1)+binom(q-1,2)=q(q-1)/2=e                              (12)
```

such factors, proving the determinant formula.

### Reduction modulo `p`

Over `F_p`, put `u=t-1`.  Since

```text
t^p-1=(t-1)^p,
```

a zero-sum root row is the ideal `u F_p[u]/(u^p)`.  Divide source and target
by `u` and work over

```text
R=F_p[u]/(u^(p-1)).                                           (13)
```

Eliminate the last of the `q` clock rows using the second margin.  Changing
the polynomial columns to

```text
(X-1),(X-1)^2,...,(X-1)^(q-1)                               (14)
```

is unimodular because `q<p`.  At the slope `tau_i`, put

```text
y_i=(1+u)^(tau_i)-1=u(tau_i+O(u)).                            (15)
```

The evaluation matrix is

```text
(y_i^j)_(i,j)
 =C(u) diag(u,u^2,...,u^(q-1)),                              (16)
```

where

```text
C(0)=(tau_i^j)_(i,j=1)^(q-1)
```

is a nonzero Vandermonde matrix because the slopes are distinct and
nonzero.  Hence `C(u)` is invertible over `R`.  Multiplication by `u^j` on
`R` has nullity `j`, so the mod-`p` nullity is exactly

```text
1+2+...+(q-1)=e.                                             (17)
```

Thus exactly `e` integral invariant factors are divisible by `p`.
Their total `p`-valuation is `e` by (10)--(12); all equal `p`.  This proves
the Smith form.

## 3. Sharp prescribed bad sets

Given any prescribed

```text
S subset F_p^*,                       |S|=q-2,
```

define in `Q(zeta)[X]`

```text
P_S(X)=(X-1) product_(s in S)(X-zeta^(-s))
      =sum_(k=0)^(q-1) D_k X^k,                              (18)

f_S(k,v)=Tr_(Q(zeta)/Q)(D_k zeta^v).                         (19)
```

The coefficients `D_k` are cyclotomic integers, so `f_S` is integral.
Both margins vanish:

```text
sum_v f_S(k,v)
 =Tr(D_k sum_v zeta^v)=0,

sum_k f_S(k,v)
 =Tr(P_S(1) zeta^v)=0.                                      (20)
```

Moreover, for `alpha!=0`,

```text
widehat(T_tau f_S)(alpha)
 =p sigma_alpha(P_S)(zeta^(-alpha tau)).                     (21)
```

Therefore the nonzero bad slopes are **exactly** `S`.  Since `D_(q-1)=1`,
the witness is nonzero.  It has precisely `p-q+1` good slopes, proving
sharpness in (5).

The literal trace in (19) is load-bearing.  THM-2507's convenient
single-character coefficient lift is only row-zero and need not have the
second margin.  It cannot be copied into this theorem.

## 4. Rational Galois saturation

For rational `f`, every nonzero tooth profile fires all nontrivial root
characters:

```text
T_tau f !=0
 iff
widehat(T_tau f)(alpha)!=0 for every alpha in F_p^*.          (22)
```

Indeed the output is a rational zero-sum vector.  If one nontrivial Fourier
coefficient vanished, its coefficient polynomial would be divisible by
`Phi_p`.  Its degree is at most `p-1`, so it would be a scalar multiple of
`Phi_p`; the zero sum makes that scalar zero.

Rationality is essential.  Over `C`, a single character is zero-sum and has
only one Fourier colour.

## 5. Specialization to `C_7 x F_13`

For `(p,q)=(13,7)`,

```text
dim V_00=72,

any 6 nonzero teeth:          rational isomorphism,

every nonzero interaction:   at least 7 of 12 teeth,

integral determinant:         13^21,

integral Smith form:          1^51 direct-sum 13^21.          (23)
```

On each nontrivial root character, the twelve slope values form a
generalized Reed--Solomon `[12,6,7]` code: six symbols reconstruct and seven
is the sharp minimum support.  The full rational interaction is the Galois
compatible collection of the twelve character slices.

The exact companion exhausts all

```text
binom(12,5)=792 five-slope banks,
binom(12,6)=924 six-slope banks.                              (24)
```

Every five-bank has rank `60`.  Every six-bank has rational rank `72` and
rank `51` modulo thirteen.  All `792` literal trace witnesses have their
prescribed five-slope zero set; their `L^1` norms range from `396` to `664`.

Conditional on THM-2550(B)'s current proved-candidate exact input, each of
its nonzero rational doubly centered interactions `d_M,d_C` therefore
survives at least seven tooth slopes, and any prescribed six profiles
reconstruct it.  This uses only Part (B); THM-2550(A)'s positive drift is a
different clock and packet.

## 6. Ordered-section shear and holonomy

The construction uses the representatives `k=0,...,6`, not an affine map
from `C_7` to `F_13`.  On this ordered set define

```text
(S_sigma f)(k,v)=f(k,v-sigma k).                              (25)
```

Then exactly

```text
T_tau S_sigma=T_(tau+sigma).                                 (26)
```

This identity lives on the ambient **row-zero** strip.  It is not an
automorphism of `V_00`: the shear generally destroys the column margin.  For
the rectangle

```text
e_(0,0)-e_(0,12)-e_(6,0)+e_(6,12),
```

`S_1` has column-margin vector

```text
(1,0,0,0,0,1,-1,0,0,0,0,0,-1).                              (26a)
```

The coordinate changes between two six-banks on `V_00` are the abstract
rational isomorphisms

```text
T_(S') T_S^(-1),                                             (26b)
```

not the shears (25).  Cyclically moving the owner/clock row does preserve
both margins, but crosses the section seam.  For a row shift `c`,

```text
iota(k+c mod 7)=iota(k)+c-7 q_c(k),                           (27)
```

where the wrap `q_c(k)` depends on the row.  The tooth coordinate therefore
acquires the term

```text
-7 tau q_c(k)                     mod 13.                     (28)
```

It is not a common root translation.  For `c=tau=1`, the wrap vector is

```text
(0,0,0,0,0,0,1),                                           (29)
```

and the companion gives an exact doubly centered witness whose rotated
tooth is not any cyclic translate of the original tooth.

Equations (26)--(29) are the precise section-holonomy boundary.  Six-bank
coordinates glue abstractly on the signed space, but the tempting geometric
shear already leaves that space, and cyclic owner motion introduces the
nonconstant wrap.  A physical gluing must retain the ordered cut and carry.

## 7. Two sharp composition no-goes

### Tooth slope is not target chord

The slope `tau` connects owner/source phase `ell` to intervention co-shift
`s` in the signed THM-2550(B) interaction.  THM-2559's chord displacement
connects an occupied old-predecessor root to the already named literal
target-role failure on a Boolean base.  The two `F_13` labels have different
domains and actions.  THM-2561's physical blind necklace remains a valid
warning: thresholding a signed tooth profile does not create the target
failure mask.

### Tomography does not defeat six replicas

Fix a six-slope bank and a nonzero rational scalar `c`.  The centered profile

```text
b(0)=-6c,
b(1)=...=b(6)=c,
b(7)=...=b(12)=0                                             (30)
```

lies in `V_0(F_13)`.  By surjectivity in (4), there is a unique rational
signed `f` whose output is this same `b` in all six chosen slopes.  After a
root rescaling, (30) is exactly the six-equal-replica shape demanded of a
residual cancelling THM-2562's anchored duty term.

Therefore more signed tooth data cannot refute the THM-2562 cancellation
invoice.  One needs a physical positivity/common-carrier constraint on the
covariance residual, or the carry-labelled semantic 2-cell that identifies
the future root.  The atlas supplies neither.

## 8. Exact referee and stopping boundary

Run

```bash
python3 04-computation/lrc14_six_tooth_v00_tomography_thm2564.py
python3 -O 04-computation/lrc14_six_tooth_v00_tomography_thm2564.py
```

Both executions must reproduce

```text
05-knowledge/results/lrc14_six_tooth_v00_tomography_thm2564.out
```

byte-for-byte.  The dependency-free companion constructs the `72`
rectangle basis, verifies every bank in (24), checks rank `51` modulo
thirteen for all `924` six-banks, computes three independent `72 x 72`
Bareiss determinant controls equal to `13^21`, exhausts the `792` literal
trace hostiles and all their good root colours, checks all `169` ambient
ordered-set shear identities, exhibits their second-margin hostile, and
verifies the nonconstant seam carry.

The theorem reconstructs the whole rational interaction, determines its
integral index, and identifies the exact chart holonomy.  It does not make a
signed current nonnegative, canonically identify a target chord, preserve a
THM-2305 word through the section seam, produce an old-head/later-root map on
one ancestry base, exclude a scalar row, or prove LRC(14). **QED.**
