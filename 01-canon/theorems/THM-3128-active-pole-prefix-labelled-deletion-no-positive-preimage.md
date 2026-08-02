---
id: THM-3128
title: "Active pole-prefix labelled-deletion no-positive-preimage boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The first
  active THM-3120 pole-prefix Hasse hostile, at I2 support (1,3) and degree
  five, has no nonnegative Hasse-boundary preimage with the same THM-3119
  labelled-deletion image.  This remains true after adding an arbitrary
  vector from the full two-dimensional labelled-deletion kernel.  The fixed
  negative upset is exactly one raw block-deletion coordinate.  Selectors
  that change the deletion image or leave the Hasse cone remain open.
source: root/multiscale-newton-flag/low-child-flag-extension-2026-08-02
audit: >
  An independent hostile audit reconstructed the labelled-deletion matrix,
  determinant-360 rank certificate, full two-dimensional kernel, factorial
  gauge, raw conjugated deletion, invariant top-two upset, and the resulting
  no-positive-preimage implication.  Fresh normal and optimized runs
  byte-match the stored output; LF hashes and documentation checks pass.
depends_on:
  - THM-3119-factorial-normalized-labelled-deletion-and-young-carrier-order
  - THM-3120-row-pole-prefix-newton-flag-positivity
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
related:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-3122-labelled-deletion-positive-kernel-ghost-and-no-upward-induction
script: 04-computation/gmc_active_prefix_labelled_deletion_no_preimage_thm3128.py
output: 05-knowledge/results/gmc_active_prefix_labelled_deletion_no_preimage_thm3128.out
script_sha256: 2af848da181af8a674a8e9c1d6600ef4b3b063ada9954abb9c3cbf39c12481cd
output_sha256: fdb4ea35183445d7f48b605e828abb04747b7f05265400495ac634725a9bb037
hash_basis: LF-normalized bytes
---

# THM-3128 -- active pole-prefix labelled-deletion no-positive-preimage boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3120 proves scalar positivity of every active pole-prefix Newton
coordinate in its exact support bank.  The scalar flag does not always lift
to the stronger THM-3115 cone of nonnegative fine-to-coarse Hasse boundaries.
THM-3122 explains abstractly why labelled deletion can hide a positive or
negative kernel component.  The natural next repair is therefore to add an
arbitrary labelled-deletion kernel correction.

For the first active hostile, even that repair is impossible.  The full
kernel is two-dimensional, and a violated coarsening upset is completely
invisible to both kernel directions.  Equivalently, its negative mass is one
coordinate of the raw deletion image and cannot change within the fibre.

## 1. The active pole-prefix current

Use the `I_2` residual bank at

```text
(a,b)=(1,3).                                                   (1)
```

Its reduced row denominator has poles

```text
(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1),                          (2)
```

and its reduced numerator is `t^5P(t)`, where

```text
P(t)=1440-16440t+51264t^2-37176t^3+84240t^4-325248t^5.        (3)
```

Thus `d=5`.  Remove the active top prefix

```text
R=(8,7,6,5,5).                                                (4)
```

For the signed bank functional `Phi` and distinguished residual alphabet
`Q`, put

```text
Phi^R(f)=sum_S epsilon_S f[S-R],             Q^R=Q-R.          (5)
```

At degree five define the raw Young-gap coefficient current

```text
G_mu=Phi^R(h_5)m_mu[Q^R]-Phi^R(m_mu)h_5[Q^R].                 (6)
```

Here `G_mu` is the coefficient of `Lbar_mu`, not yet a factorial coordinate.
In the order

```text
(5),(4,1),(3,2),(3,1,1),(2,2,1),(2,1,1,1),(1^5),             (7)
```

exact evaluation gives

```text
G=(-2324160, 544320, 2237760, -656640,
   -915840, 972000, 142560),                   sum_mu G_mu=0.  (8)
```

The scalar flag itself is positive:

```text
Phi^R(h_5)=1440>0.                                            (9)
```

The coarsest sign is lost through the virtual reference alphabet, not the
row scalar:

```text
Phi^R(p_5)=0,              p_5[Q^R]=-1614,
G_(5)=1440(-1614)=-2324160.                                  (10)
```

The companion reconstructs `(2)--(10)` from the THM-3110/3120 banks; the
current is not inserted as an unexplained fixture.

## 2. The two gauges

THM-3119 labelled deletion acts on factorial coordinates.  Put

```text
w_mu=product_i mu_i!,
W_5=diag(120,24,12,6,4,2,1),
c=W_5^(-1)G.                                                  (11)
```

Thus

```text
c=(-19368,22680,186480,-109440,
   -228960,486000,142560).                                   (12)
```

The associated carrier is unchanged:

```text
sum_mu c_mu Kbar_mu=sum_mu c_mu w_mu Lbar_mu
                    =sum_mu G_mu Lbar_mu.                     (13)
```

This gauge distinction is load-bearing.  Adding a kernel vector directly to
the raw vector `G` would apply labelled deletion in the wrong coordinates.

## 3. The full labelled-deletion kernel

Order the degree-four target types as

```text
(4),(3,1),(2,2),(2,1,1),(1^4).                               (14)
```

In the source and target orders `(7),(14)`, THM-3119 labelled deletion is

```text
    [5 1 0 0 0 0 0]
    [0 4 2 2 0 0 0]
A = [0 0 3 0 1 0 0] .                                        (15)
    [0 0 0 3 4 3 0]
    [0 0 0 0 0 2 5]
```

The minor on columns `0,1,2,3,5` has determinant `360`, so `rank A=5`.
Two independent kernel vectors are

```text
K_1=(-1,5,-2,-8,6,0,0),
K_2=(1,-5,0,10,0,-10,4).                                    (16)
```

Direct multiplication gives `AK_1=AK_2=0`.  Since the nullity is two,
`K_1,K_2` span the full kernel.  The second vector is the degree-five hook
ghost `kappa_5` from THM-3122; the first is the missing independent
direction.

Every factorial-coordinate current with the same labelled-deletion image as
`c` is therefore

```text
c'=c+alpha K_1+beta K_2.                                     (17)
```

Returning to raw `Lbar` coefficients, `(17)` becomes

```text
G'=G+alpha WK_1+beta WK_2,                                   (18)

WK_1=(-120,120,-24,-48,24,0,0),
WK_2=(120,-120,0,60,0,-20,4).                                (19)
```

## 4. Universal top-two deletion cut

The invariant used below exists in every degree.  Let `N>=3` and let `B_N`
be raw unweighted block deletion from partitions of `N` to partitions of
`N-1`.  A source partition can lower to the one-part target `(N-1)` only in
two ways:

```text
(N)       -> (N-1)       by lowering its unique part,
(N-1,1)   -> (N-1)       by deleting its unique singleton.    (19a)
```

Both raw coefficients are one.  Therefore the `(N-1)` row of `B_N` is

```text
[e_(N-1)]B_N e_lambda
 =1_{U_N}(lambda),
U_N={(N),(N-1,1)}.                                           (19b)
```

The set `U_N` is the top-two coarsening upset.  Hence every raw Hasse
boundary `Y` satisfies

```text
(B_NY)_(N-1)=Y(U_N)>=0.                                      (19c)
```

It follows universally that

```text
(B_NG)_(N-1)<0
 ==> no raw Hasse boundary Y has B_NY=B_NG.                   (19d)
```

Via THM-3119 factorial conjugacy, `(19d)` is equally a no-preimage theorem
for the corresponding labelled-deletion fibre.  This is a reusable
isotone-row obstruction; the rest of the proof identifies its first active
pole-prefix instance.

## 5. The invariant negative upset

Consider the coarsening upset

```text
U={(5),(4,1)}.                                                (20)
```

The hostile and both full-kernel directions have masses

```text
G(U)=-1779840,
(WK_1)(U)=0,                 (WK_2)(U)=0.                     (21)
```

Consequently every current in `(18)` satisfies

```text
G'(U)=-1779840<0.                                             (22)
```

By THM-3127 finite Strassen--Hasse duality, no such `G'` is a nonnegative
fine-to-coarse Hasse boundary.

There is an even shorter structural reading.  Let

```text
W_4=diag(24,6,4,2,1),
B=W_4 A W_5^(-1).                                             (23)
```

Then

```text
    [1 1 0 0 0 0 0]
    [0 1 1 2 0 0 0]
B = [0 0 1 0 1 0 0] ,                                        (24)
    [0 0 0 1 2 3 0]
    [0 0 0 0 0 1 5]
```

which is raw unweighted block deletion.  Its `(4)` row is exactly the
indicator of `U`.  Hence

```text
(BG)_(4)=G(U)=-1779840.                                      (25)
```

Same labelled-deletion image in factorial coordinates is equivalent, after
the invertible gauges, to same `B` image in raw coordinates.  Equation `(25)`
therefore fixes the negative upset mass throughout the entire affine fibre.
The no-preimage result is an invariant-cut theorem, not a failure of one
chosen kernel basis or one attempted max flow.

## 6. Exact conclusion and scope

For the active prefix `(1)--(4)`, there is no raw current `G'` such that

```text
A W_5^(-1)G' = A W_5^(-1)G
```

and `G'` is a nonnegative fine-to-coarse Hasse boundary.  This quantifies over
the **full** labelled-deletion kernel.

The theorem does **not** prove that every active pole prefix has this
obstruction.  It does not rule out a selector which changes the labelled-
deletion image, changes the pole/reference alphabet, couples different
supports or degrees, or proves positivity in the larger Young-gap operator
PSD cone without entering the Hasse cone.  Failure of a sufficient Hasse
certificate is not a proof that the represented operator is non-PSD.

In particular, the result does not weaken THM-3120 scalar row positivity and
does not close or refute arbitrary anchored product-Gamma width three, SFC,
the Gaussian Moment Conjecture, NC2, LRC(14), JC(2), or DC(2).  It identifies
the exact missing datum more sharply: a successful transverse selector must
alter the deletion-visible `(4)` coordinate, not merely choose another point
of the existing labelled-deletion fibre.

## 7. Exact companion

Reproduce with

```text
python 04-computation/gmc_active_prefix_labelled_deletion_no_preimage_thm3128.py
python -O 04-computation/gmc_active_prefix_labelled_deletion_no_preimage_thm3128.py
```

The companion reconstructs the virtual prefix current, both gauge vectors,
the exact `A` and `B` matrices, the rank-five minor, the full two-dimensional
kernel in both gauges, the invariant upset masses, and the conjugacy of the
two deletion images.  Both executions must byte-match the stored output.

**QED.**
