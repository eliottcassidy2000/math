---
id: THM-2944
title: "Width-nine/ten two-chart Macaulay resultant closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  n>=0, first-window SFC(4) holds on all 64
  translated four-slot supports of exact width nine or ten.  The
  original and stable-mutated Macaulay charts have a common factor
  associated to q200^5*c300*g*Res, where g=gcd(K,P_alt); it has
  nonnegative coefficients and positive constant term, while the
  residual chart cofactors are coprime.  Independently, for the
  primitive positive-leading associate H~g*Res, the radical of
  gcd(H,H') is an explicit product of negative-root
  linear factors depending only on the first gap.  The coefficientwise
  common-factor gate is the primary closure; the repeated-divisor
  identity is a negative-root sidecar, not a bare-resultant formula.
  No width-eleven or arbitrary-width claim is made.
source: codex-gmc-uniform-width-extension-2026-07-29
audit: >
  An independent hostile audit rederived the two separate common-factor
  gates, checked that only positivity of the full
  q200^5*c300*g*Res content excludes the genuine resultant wall,
  verified the g*Res rather than bare-Res typing of the repeated-divisor
  sidecar, recomputed the source and output hashes, and replayed normal,
  optimized, and stored output byte-for-byte with empty stderr.  It
  found no remaining defect.
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
  - THM-2945-nonnegative-complete-intersection-norm-and-repeated-divisor-gate
related:
  - THM-2849-four-slot-first-window-macaulay-box
  - THM-2924-diameter-six-macaulay-newton-atlas
  - THM-2927-general-width-flagged-macaulay-leading-coefficient
script: 04-computation/gmc_width_nine_ten_two_chart_resultant_closure_thm2944.py
output: 05-knowledge/results/gmc_width_nine_ten_two_chart_resultant_closure_thm2944.out
script_sha256: 4332ca487e96dff3a63fc3b83089e517a6bc6eed1762e3c3e57a54f8aaeeb0d2
output_sha256: 3830429671604f2cfd4b930f7e94cb7911beb687ae3c0f377c359962fa1aa964
two_chart_dependency_sha256: d2f8afeba7dd6c7950405a4845d7bf112b6c9872dd8161146446be8bbdaae0ba
hash_basis: LF-normalized bytes
---

# THM-2944 -- width-nine/ten two-chart Macaulay resultant closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let

```text
L:C[s] -> C,                         L(s^j)=j!.          (1)
```

Fix

```text
M in {9,10},             0<a<b<M,               n>=0. (2)
```

For every nonzero

```text
H_0=c_0 s^n+c_1 s^(n+a)+c_2 s^(n+b)+c_3 s^(n+M),    (3)
```

at least one of

```text
L(H_0),              L(H_0^2),              L(H_0^3), L(H_0^4) (4)
```

is nonzero.  Thus first-window SFC(4) holds on all

```text
C(8,2)+C(9,2)=28+36=64                               (5)
```

translated four-slot supports of exact width nine or ten.

## 2. The inherited two charts

Normalize by `f_j=s^j/j!`, eliminate the mean against the top slot,
and use PROVED THM-2925 to clear the order-`m` form by

```text
D_(m,M)
 =[prod_(j=1)^(M-1)(n+j)^(m-1)](n+M)^(m-2).          (6)
```

Call the resulting homogeneous forms of degrees `2,3,4`

```text
Q,                              C,                   F. (7)
```

The factors in `(6)` are positive for `n>=0`, so the scaling preserves
projective common zeros.

Use the degree-seven Macaulay map

```text
S_5 direct_sum S_4 direct_sum S_3 -> S_7             (8)
```

and the two row charts of PROVED THM-2943:

```text
Q rows: 0,...,19,
C rows: 21,...,29,35,

J0 local F rows: {0,1,2,3,4,5},
J1 local F rows: {0,3,4,5,6,7}.                      (9)
```

Their global quartic rows are respectively

```text
{36,37,38,39,40,41},
{36,39,40,41,42,43}.                                 (10)
```

Write the determinants as `Delta_0,Delta_1`.  THM-2925 gives

```text
deg Delta_i<=58M-36.                                  (11)
```

## 3. The two-gate common-factor proof

PROVED THM-2942 gives the global identities

```text
Delta_0=q200^6*c300*K*R,
Delta_1=q200^5*c300*P_alt*R,                          (12)
```

where

```text
R=Res(Q,C,F)                                          (13)
```

is the ternary resultant and `K,P_alt` are the explicit flag
coordinates in THM-2942 `(7),(26)`.

For one support put

```text
g=gcd(K,P_alt),              K=gK_0,       P_alt=gP_0. (14)
```

The exact **flag-separation gate** is

```text
gcd(K_0,P_0)=1,
gcd(q200*K_0,P_0)=1.                                  (15)
```

The second equality retains the unmatched extra `q200` in the
original chart.  Equations `(12)--(15)` identify the complete common
factor over `Q[n]`, up to a positive rational scalar, as

```text
G ~ q200^5*c300*g*R.                                  (16)
```

Choose the positive-leading associate `G`.  The exact
**common-resultant-content gate** is

```text
G(0)>0,
every coefficient of G is nonnegative.                (17)
```

Thus

```text
G(n)>0                              for every n>=0.    (18)
```

If `R(n_0)=0` at a nonnegative integer, `(16)` would give
`G(n_0)=0`, contradicting `(18)`.  Therefore

```text
Res(Q,C,F)(n)!=0                    for every n>=0.    (19)
```

The resultant criterion now says that `Q,C,F` have no common
projective zero.  Restoring the eliminated mean proves `(4)`.

This is the primary proof.  It uses both gates: cofactor coprimality
identifies the whole common factor, while coefficientwise positivity
removes that factor, including the genuine resultant, from the
nonnegative depth ray.

## 4. Independent negative-root sidecar

Remove scalar content and choose the positive-leading associate

```text
H ~ g*R.                                               (20)
```

The exact replay independently proves

```text
rad gcd(H,H')
 = [prod_(j=1)^a (2n+2j-1)]
   [prod_(r=a)^M (n+r)].                              (21)
```

The right side is independent of the second interior offset `b`.
Every one of its roots is strictly negative.

Equation `(21)` is typed deliberately: it is a formula for the
associate `H~g*Res`, not for the bare resultant.  The positive extraneous factor
`g` can contribute negative repeated roots.  PROVED THM-2945 explains
why `(21)`, together with the positive endpoint value, is an
independently useful repeated-divisor certificate for the nonnegative
complete-intersection norm.  No division of `g` and no bare-resultant
discriminant formula is asserted here.

The sidecar is not needed for the closure in Section 3.  Its value is
structural: it isolates all repeated divisors on the negative ray and
exposes the striking first-gap-only law behind the finite atlas.

## 5. Exact 64-family verification

For each normalized support in `(5)`, the companion:

1. interpolates `Delta_0` from the complete prefix required by `(11)`;
2. divides exactly by `q200^6*c300*K` to recover `R`;
3. reconstructs `Delta_1` using the proved second identity in `(12)`;
4. verifies `(15)--(18)` and checks the reconstructed raw reduced
   chart cofactors are coprime;
5. computes `gcd(H,H')`, removes repeated multiplicities, and verifies
   `(21)` exactly;
6. compares both charts with direct `36`-by-`36` determinants at the
   three outside-interpolation depths

   ```text
   58M-35, 58M-34, 58M-33,                            (22)
   ```

   for `64*3=192` paired controls; and
7. checks the same three-term Grassmann--Pluecker exchange as
   THM-2943 once per family, for `64` exact exchanges.

The original reduced chart has coefficients of both signs in

```text
width 9:  22 of 28,
width 10: 30 of 36,
total:    52 of 64.                                  (23)
```

This is why a one-chart coefficientwise argument is not the right
invariant.  The mutation removes every common flag wall without
asserting that each mixed cofactor has an integral root.

The family digests are

```text
width 9:
a75989df9b63e477239d912caf5acdfe643ef4424d9bd8a327a23fe3b7d20ca4,

width 10:
598ce2948dae63a93f70ab42f8a57913f8c48fc732ee1e95cd94b92dd702c8d2,

combined:
af2950ea595583411a602b7bf0258e97cc20687de52d852007ade88322b3744a.
                                                               (24)
```

## 6. Scope and stopping boundary

The proved-candidate scope is exactly:

```text
slots:        four distinct exponents;
widths:       exactly nine or ten;
window:       moments one through four;
depth:        every integer n>=0;
certificate:  two Macaulay charts plus full common-factor positivity.
                                                               (25)
```

No exact width eleven, arbitrary-width SFC(4), shifted moment window,
SFC(5), arbitrary-radial GMC(2), or Jacobian-conjecture conclusion is
claimed.  The finite first-gap-only pattern `(21)` is not extrapolated
past the audited bank.

## 7. Exact companion

The companion hash-pins the independently audited THM-2943 two-chart
artifact, which in turn pins THM-2942 and THM-2925.  It uses explicit
`require` checks under both ordinary and optimized Python.  Run

```text
python 04-computation/gmc_width_nine_ten_two_chart_resultant_closure_thm2944.py
python -O 04-computation/gmc_width_nine_ten_two_chart_resultant_closure_thm2944.py
```

Both modes must byte-match the stored transcript and the declared
LF-normalized hashes.

**QED.**
