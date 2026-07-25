---
id: THM-2166
title: "Hybrid whole-core smoothing and low-carry crossing"
status: >
  PROVED from THM-2145's Jackson polynomial, a 1,716-core exact
  mass/component sweep, and THM-2163's radix identity. For a defect-six
  split with seven retained speeds in {1,...,13}, zero safe measure forces a
  crossing relation whose far coefficients have height 298 and are 7-units,
  whose core has support at most two and height 57, and whose nonzero cut
  carry has magnitude at most 708. Core height 57 is certificate-minimal for
  representing every possible carry in this range. Thus either one far speed
  is at most 708
  or at least two far speeds have a height-298 near-cancellation of magnitude
  at most 708. After four dyadic digits the full relation carry has magnitude
  at most 1787 and only six far owners remain. This does not finish defect
  six or LRC(14).
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-2145
  - THM-2163
related:
  - THM-2054
  - THM-2162
script: 04-computation/lrc14_hybrid_core_low_carry_referee_codex_20260724.py
output: 05-knowledge/results/lrc14_hybrid_core_low_carry_referee_codex_20260724.out
script_sha256: 025c9d2a516587c12bba120fce06bb6579869bc0bf37577b2d507f336703c4a8
output_sha256: 044c687196e37851f81ccbee6695cc01b780d7fef022e348326649726e9c1e30
hash_basis: working-tree bytes (LF)
---

# THM-2166 -- hybrid whole-core smoothing and low-carry crossing

At radius `1/14`, put

```text
G={x in R/Z:||x||>=1/14},
G_A={t:at in G for every a in A}.                      (1)
```

Let

```text
E subset {1,...,13},       |E|=7,                      (2)
F=(f_1,...,f_6)                                         (3)
```

where the full thirteen speeds are distinct positive integers. Assume

```text
mu(G_E intersect G_F)=0.                               (4)
```

THM-2145 smooths both blocks coordinate by coordinate. Here the far block is
still treated that way, but the highly structured core is first assembled
as one interval union and only then smoothed. This preserves its signed
boundary cancellation and makes the crossing frequency itself small.

## 1. The six-factor far polynomial

Use THM-2145's Jackson notation

```text
J_N=F_N^2/integral F_N^2,
q_N=J_N*1_G,
eta_N=2 integral_T ||x||J_N(x)dx.                     (5)
```

Then

```text
0<=q_N<=1,
degree(q_N)<=2N-2,
||q_N-1_G||_1<=eta_N.                                 (6)
```

Take `N=150` and define

```text
P_F(t)=product_(i=1)^6 q_150(f_i t).                  (7)
```

The product telescope and Haar invariance under integer dilation give

```text
||P_F-1_(G_F)||_1<=6eta_150.                           (8)
```

THM-1234, as routed through THM-2145, gives

```text
mu(G_F)>=alpha:=61/273.                                (9)
```

Consequently

```text
integral P_F>=alpha-6eta_150.                          (10)
```

Every factor has coefficient-index support `|k|<=298`, producing scalar
modes `k f_i`. Moreover, by THM-2145,

```text
q_150 hat(k)!=0 for 0<|k|<=298
iff 7 does not divide k.                               (11)
```

## 2. Smooth the whole core

Let

```text
beta_E=mu(G_E),
K_E=number of positive-length circle-interval
    components of G_E.                                (12)
```

Isolated weak-safe points are deliberately omitted from `K_E`: their
indicators vanish in `L^1`. Define

```text
R_E=J_355*1_(G_E).                                     (13)
```

Then

```text
0<=R_E<=1,
integral R_E=beta_E,
degree(R_E)<=708.                                      (14)
```

For a union of `K_E` circle intervals,

```text
mu(G_E triangle (G_E+x))<=2K_E||x||.                  (15)
```

Indeed, translating by `x` moves each of the `2K_E` interval endpoints
through length at most `||x||`; overlaps can only reduce the symmetric
difference. Minkowski's inequality and Fubini therefore give

```text
||R_E-1_(G_E)||_1
 <=integral J_355(x)
      mu(G_E triangle (G_E+x))dx
 <=K_E eta_355.                                       (16)
```

This is the exact analytic value of retaining the core as one object.

## 3. Exact all-core margin

The boundary arrangement of the danger combs for speeds `1,...,13` has
`178` rational points and `177` open cells. Exact midpoint evaluation on
those cells gives `(beta_E,K_E)` for all

```text
binomial(13,7)=1716                                    (17)
```

cores. THM-2145's rational-`pi` Jackson formula gives the compact bounds

```text
eta_150<439/156250,
eta_355<371/312500.                                    (18)
```

For every core, the exact sweep verifies

```text
(alpha-6eta_150)beta_E
  >6eta_150+K_E eta_355.                               (19)
```

It is enough to insert the upper bounds in (18). The unique minimum of that
stronger rational ledger occurs at

```text
E_*=(1,5,7,8,9,11,13),
beta_(E_*)=45107/229320,
K_(E_*)=20,                                            (20)
```

and its exact remaining margin is

```text
41050267/1222741406250>0.                              (21)
```

The same rational-`pi` ledger is negative with core kernel `J_354`; hence
`355` is the first of these two adjacent degrees to close this certificate.
No optimality among other kernels is claimed.

## 4. Extract the low carry

By (8), (16), and (4),

```text
integral P_F R_E
 <=6eta_150+K_E eta_355.                               (22)
```

Equations (10), (14), and (19) imply

```text
(integral P_F)(integral R_E)>integral P_F R_E.         (23)
```

Finite Fourier orthogonality now forces a nonzero common scalar frequency:

```text
P_F hat(nu)!=0,       R_E hat(-nu)!=0,
0<|nu|<=708.                                           (24)
```

Expanding the first coefficient as a finite convolution supplies

```text
u=(u_1,...,u_6) in [-298,298]^6,
sum_i u_i f_i=nu,                                      (25)
```

where every nonzero `u_i` is prime to seven by (11).

It remains to turn `-nu` into a sparse core coefficient vector. For an
integer height `B`, define the exact two-coordinate sumset

```text
C_B(E)=union_({e,f} subset E)
       {a e+b f: |a|,|b|<=B}.                          (26)
```

An exhaustive integer bitset computation over all `1716` cores proves

```text
[-708,708] subset C_57(E)          for every E.        (27)
```

This universal height is sharp for the unrestricted frequency range.
At height `56` exactly three cores fail:

```text
(1,2,3,4,5,6,7)       misses +/-699,+/-705,+/-706,
(1,2,3,4,5,6,8)       misses +/-701,
(1,2,4,5,6,8,10)      misses +/-701.                  (28)
```

Thus (24) and (27) supply a vector `w` such that

```text
sum_(e in E)w_e e=-nu,
|supp w|<=2,       ||w||_infinity<=57.                 (29)
```

Combining (25) and (29) proves the crossing relation

```text
sum_i u_i f_i+sum_(e in E)w_e e=0,                    (30)
sum_i u_i f_i=nu!=0.                                  (31)
```

The support-two conclusion uses only the integer `nu` from (24), not a
claim about the Fourier support of `1_(G_E)` or `R_E`. The core coefficients
need not be `7`-units; only the far coefficients inherit (11).

## 5. Defect-six and radix consequences

If exactly one far coefficient in (25) is nonzero, then

```text
|u_i|f_i=|nu|<=708,
```

so that far speed itself satisfies `f_i<=708`. Otherwise at least two far
speeds have a height-298, `7`-unit integer combination whose nonzero
magnitude is at most `708`. Thus

```text
one bounded far speed
or
two-or-more-term bounded near-cancellation.            (32)
```

This replaces THM-2145's raw core carry bound `20860` by the intrinsic cut
carry `|nu|<=708`.

The full coefficient vector `m=(u,w)` satisfies

```text
||m||_1<=6*298+2*57=1902.                              (33)
```

THM-2163 therefore gives the generic radix bound `|kappa_j|<1902`.
There is a sharper dyadic statement. Divide only the far speeds:

```text
F=2^j Z_j+R_j,       D_j=Z_j mod 2,
lambda_j=(u.R_j-nu)/2^j=-u.Z_j.                       (34)
```

Then `lambda_j` is integral,

```text
lambda_0=-nu,
lambda_(j+1)=(lambda_j+u.D_j)/2,                       (35)
```

and it terminates at zero. Once `j>=4`, every core speed in (2) has
terminated. The full relation carry equals `lambda_j`, and

```text
|lambda_j|
 <=((2^j-1)||u||_1+|nu|)/2^j
 <=1788-1080/2^j
 <1788.                                                (36)
```

Hence

```text
|kappa_j|<=1787               for j>=4,                (37)
```

with only the six far-coordinate owner mask still active. This is `3575`
possible carry values, but the owner mask and unbounded digit depth remain;
THM-2163's termination obstruction is not removed.

The theorem therefore sharpens the defect-six state dramatically without
making it finite. In the multi-far branch, (32) still permits arbitrarily
large speeds and does not by itself prove LRC(14).

## 6. Exact referee

The companion reconstructs the full `178`-boundary arrangement, counts
positive-length components rather than closed-list entries, evaluates both
Jackson moments using exact `Fraction` arithmetic, proves (18)--(21), and
checks the adjacent `N=354` failure. An independent pair-sum bitset path
checks (27)--(28), including all `1417` possible signed carries per core.
The script uses explicit raising checks, so optimized Python preserves the
validity gate. Normal and optimized transcripts agree with the stored output.
