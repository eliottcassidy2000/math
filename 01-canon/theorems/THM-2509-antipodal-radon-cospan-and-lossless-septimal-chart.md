---
id: THM-2509
title: "Antipodal Radon cospan and the lossless septimal chart"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED; FINITE-EXACT ATLAS
  SHARPENING. Every nonvertical row-zero array on F_13 x {0,...,6}
  has a complete antipodal pair of nonzero truncated-Radon slopes. The
  joint two-leg chart is injective, recovers the whole array, has those
  two Radon transforms as its marginals, and turns slope reversal into
  leg swap. CRT translation is conjugated faithfully to a permutation of
  the 91-point strip, retaining rather than erasing the septimal carry.
  Some fixed pair survives on at least one sixth of any essential locus;
  the complete THM-2436 atlas improves this to five sixths and gives
  parent-locus floors 5/21 and 5/14. These are pointwise nonvanishing
  bounds, not integrated-current bounds. The packet belongs to an already
  empty high-septimal branch and does not produce a lawful owner, target,
  terminal-word, or deep current on any live row. LRC(14) remains open.
source: codex-2026-07-27-antipodal-radon-cospan
depends_on:
  - THM-2435-top-blocker-essential-parent-and-punctured-stalk-carrier
  - THM-2507-truncated-radon-toothpick-tomography-and-nonaffine-root-boundary
related:
  - THM-2368-owner-pivot-root-fibre-radon-invertibility
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
  - THM-2506-punctured-stalk-primitive-module-saturation-and-thirteen-primary-pushforward-no-go
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
script: 04-computation/lrc14_antipodal_radon_cospan_thm2509.py
output: 05-knowledge/results/lrc14_antipodal_radon_cospan_thm2509.out
script_sha256: 58cb8f241349924c48e2f4c6f0cd05d3e2be31e803b85455dd91933d0800eeb0
output_sha256: 7f31f6cd7cdea2982f834a243243f5f9c77290f1d6fadcedb3be0cb4d3dc89ed
hash_basis: working-tree bytes (LF)
---

# THM-2509 -- the first intrinsic Radon object is a cospan

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED; FINITE-EXACT ATLAS
SHARPENING.**

THM-2507 shows that a nonvertical `13 x 7` row-zero defect is visible in at
least seven of the twelve nonzero toothpick slopes.  Seven objects among six
antipodal pairs force more than a single detector: one pair survives in both
directions.  Keeping its two outputs jointly reveals the natural object that
one scalar Radon fold had hidden:

```text
one signed 13 x 7 strip
  -> one lossless two-leg chart in F_13 x F_13
  -> two nonzero antipodal Radon marginals.                         (1)
```

This cospan retains the ordered septimal section and its carry exactly.  It
does not turn either leg into a lawful physical target current.

## 1. A complete antipodal pair is forced

Let

```text
d:F_13 x {0,...,6} -> C
```

satisfy the pointwise row-zero law

```text
sum_(r=0)^6 d(h,r)=0                                             (2)
```

and suppose that `d` is not independent of `h`.  For `tau in F_13`, put

```text
R_tau d(v)=sum_(r=0)^6 d(v-tau r,r).                            (3)
```

Then there is a `tau!=0` for which

```text
R_tau d!=0                     and R_(-tau)d!=0.                 (4)
```

Indeed, choose a nonzero horizontal character `alpha` on which `d` is
active and set

```text
P_alpha(X)=sum_(r=0)^6 d_hat_r(alpha)X^r.                       (5)
```

The polynomial is nonzero, has degree at most six, and (2) gives
`P_alpha(1)=0`.  THM-2507's slice law is

```text
(R_tau d)_hat(alpha)=P_alpha(zeta_13^(-alpha tau)).             (6)
```

The twelve roots in (6) are distinct.  Since one of the at most six roots
of `P_alpha` is already `1`, at most five nonzero slopes can be bad.  Thus
at least seven slopes are good.  The twelve slopes form the six pairs
`{tau,-tau}`, so the pigeonhole principle proves (4).

This count is sharp in the abstract integral class.  THM-2507's stored
polynomial witness with bad set `{1,2,3,4,5}` has good set
`{6,7,8,9,10,11,12}` and exactly one complete antipodal pair, `{6,7}`.
That witness has `L1=56`; it is not asserted to lie in the much smaller
THM-2436 cover atlas.

## 2. The paired chart is lossless

Let `rho:F_7->{0,...,6} subset F_13` be the ordered representative section.
For a good antipodal pair define

```text
Pi_tau(h,r)=(u,v)
           =(h+tau rho(r),h-tau rho(r)).                         (7)
```

This is an injection onto the `91`-point strip

```text
I_tau={(u,v):(u-v)/(2 tau) lies in rho(F_7)}.                    (8)
```

Its inverse is explicit:

```text
h=(u+v)/2,
r=rho^(-1)((u-v)/(2 tau)).                                      (9)
```

Define the joint signed pushforward by

```text
J_tau(u,v)=sum_(Pi_tau(h,r)=(u,v)) d(h,r).                       (10)
```

Because `Pi_tau` is injective, (10) is merely a relabelling of all `91`
entries of `d`; no value, norm, or ancestry coordinate is identified.  Its
two marginals are exactly

```text
sum_v J_tau(u,v)=R_tau d(u),
sum_u J_tau(u,v)=R_(-tau)d(v).                                  (11)
```

Both are nonzero by (4).  Reversing the slope swaps the two legs:

```text
Pi_(-tau)=swap o Pi_tau,
J_(-tau)(u,v)=J_tau(v,u).                                       (12)
```

Thus the **unordered** cospan is independent of the orientation gauge.
Naming one leg source and the other arrival would consume additional data
and is not part of the theorem.

For rational or integral `d`, every nonzero marginal in (11) has all twelve
nontrivial `C_13` Fourier colours by the prime-cyclotomic argument of
THM-2507.  This remains a pointwise statement; its phases may cancel after
integration over a moving parent.

## 3. The carry is retained as a strip permutation

Consider the forward CRT translation

```text
(h,r)->(h+A,r+C),              A in F_13, C in F_7.             (13)
```

Choose `rho(C) in {0,...,6}` and write

```text
q_C(r)=floor((rho(r)+rho(C))/7),
rho(r+C)=rho(r)+rho(C)-7q_C(r).                                 (14)
```

Conjugating (13) through (7) gives

```text
u -> u+A+tau rho(C)-7tau q_C(r),
v -> v+A-tau rho(C)+7tau q_C(r).                               (15)
```

The right side is a permutation of `I_tau`, and (9) recovers both `r` and
the carry.  In midpoint/section coordinates it is simply

```text
(h,rho(r))->(h+A,rho(r+C)).                                    (16)
```

Equations (15)--(16) are the exact gain: the cospan does not forget the
normalization seam.  They are also the exact boundary.  Unless `C=0`, the
quantity `q_C(r)` depends on the septimal row, so neither leg undergoes one
common `F_13` translation.  The cospan is intrinsic in the category of the
`91`-point strip together with its conjugated action, not in the category of
two independent physical target roots.

## 4. A fixed pair before seeing the parent

Let `Y` range over a finite-measure essential locus `E`, and let `d_Y` be a
measurable family satisfying the hypotheses above almost everywhere.  For
the six unordered antipodal pairs `p`, let

```text
E_p={Y:both slopes in p are good for d_Y}.                       (17)
```

Pointwise, Section 1 gives

```text
sum_p 1_(E_p)(Y)>=1.                                            (18)
```

After integration, one fixed pair chosen before `Y` satisfies

```text
mu(E_p)>=(1/6)mu(E).                                             (19)
```

For the THM-2436 essential-parent floors `mu(E)>=(1+k)/7`, this gives the
generic absolute locus bounds

```text
k=1: 1/21,                    k=2: 1/14.                         (20)
```

These are measures of the **pointwise two-marginal nonvanishing locus**.
They are not lower bounds for an integrated signed current.

## 5. Complete-atlas sharpening

The complete THM-2436 atlas imported and audited in THM-2507 contains
`14,952` distinct nonflat defect matrices.  Every one has eleven or twelve
good nonzero slopes.  Therefore every atlas defect has at least five of its
six antipodal pairs completely good.  Equation (18) strengthens to

```text
sum_p 1_(E_p)(Y)>=5,                                            (21)
```

so one fixed pair survives on at least

```text
(5/6)mu(E).                                                     (22)
```

Combining (22) with the two THM-2435 essential-mass invoices gives

```text
k=1: 5/21,                    k=2: 5/14.                         (23)
```

The count `5/6` is sharp from the zero-mask information alone.  The atlas
boundary

```text
d(0)=-e_0+e_1,              d(12)=e_1-e_2                      (24)
```

has bad nonzero slope `1` and exactly five complete pairs.  No assertion is
made that one physical parent distribution attains equality in (22).

## 6. Why a tournament does not supply the missing orientation

There is no canonical Paley orientation on these slopes: `-1` is a square
modulo `13`, so quadratic-residue difference and ratio tests are symmetric
on each antipodal pair.  The usual additive cyclic tournament does not force
a directed triangle either.  The sharp good set

```text
{6,7,8,9,10,11,12}                                             (25)
```

is a transitive interval in that gauge while having only one complete
antipodal pair.  The useful order-two structure is therefore (12), not a
forced tournament orientation.

The same distinction explains depth parity.  Although `13^L=(-1)^L mod 7`,
odd parity acts on ordered representatives with a wrap term.  The strip and
cut record that term; a quadratic-residue slope colour does not.

## 7. Physical stopping boundary

THM-2509 proves a static, lossless coordinate theorem on a signed defect in
the already-empty high-septimal branch.  It does **not** prove:

- positivity or noncancellation after parent integration;
- a source/arrival orientation of the two legs;
- a lawful THM-2365 factor-wise target action;
- a THM-2471 owner-loop or first-collision current;
- retention of the terminal word, deep ancestry sheet, or old address; or
- exclusion of any one of the `165` live scalar rows.

The missing sidecar is now precise: realize both legs on one live parent, at
one owner, time, and terminal word, while retaining the full cut index,
normalization phase, source-versus-arrival type, and deep ancestry.  Without
that typed intertwiner, two nonzero derived marginals are not two physical
currents.

## 8. Exact companion and independent audit

Run

```bash
python3 04-computation/lrc14_antipodal_radon_cospan_thm2509.py
python3 -O 04-computation/lrc14_antipodal_radon_cospan_thm2509.py
```

Both executions reproduce the stored output byte-for-byte.  The companion
checks all `792` seven-slope subsets, obtaining complete-pair histogram
`{1:192,2:480,3:120}`; all `1,092` chart rows; `99,372` conjugated CRT
translation rows; the six fixed-pair atlas counts; and the mass arithmetic
in (22)--(23).

The independent audit rederived the polynomial root count, injection and
inverse, both marginal signs, carry formula, abstract and atlas invoices,
and all scope boundaries.  It also reran the generic THM-2507 referee under
normal and optimized Python and the full `41,379`-assignment atlas at `-O2`;
every transcript matched its stored output. QED.
