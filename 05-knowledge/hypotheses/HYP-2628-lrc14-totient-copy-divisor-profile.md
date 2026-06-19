---
id: HYP-2628
title: LRC(14) totient-copy divisor profile - exact-period packets are the squarefree mask atoms
status: OPEN
source: codex-2026-06-19-S21
depends_on:
  - HYP-2625
  - HYP-2626
  - HYP-2627
  - THM-538
related:
  - HYP-2617
  - HYP-2619
  - THM-523
  - THM-522
  - HYP-2561
  - OPEN-Q-108
---

# HYP-2628 - LRC(14) Totient-Copy Divisor Profile

## Claim

The user's copy rule

```text
sum_{d|n} c(d) = n,  c(d) >= 1
```

has the unique solution

```text
c(n) = phi(n).
```

This is not merely a number-theory aside.  It gives the canonical measure on
the LRC14 squarefree divisor-profile address:

```text
q-grid residues -> reduced denominator d|q -> exact-period packet of size phi(d)
```

Thus prime masks should be read as exact-period packets, not as raw divisor
counts.  The useful convolution is

```text
1 * phi = id,       phi = mu * id.
```

This upgrades the HYP-2625/HYP-2626 language:

```text
prime-mask inclusion-exclusion over {2,3,5,7}
  = Mobius inversion on exact-period q-grid packets
  = the finite transfer before the signed mod-7 coimage tail.
```

## Computation

Script:

- `04-computation/lrc14_totient_copy_operator_codex_s21.py`
- output: `05-knowledge/results/lrc14_totient_copy_operator_codex_s21.out`

The script verifies the recursive copy rule against `phi(n)` for `n<=80`,
checks the exact-period residue census on `q=14,30,210`, computes compressed
prime-mask profiles for the LRC14 denominators, and records Tournament
Analysis for the proof quotients and mask-mass profiles.

## Exact Findings

The exact-period interpretation is literal.  For every `q`, residues `a/q`
split by reduced denominator `d|q`; each packet has size `phi(d)`.

For the LRC14 seam:

```text
q=14=2*7:
d=1:1, d=2:1, d=7:6, d=14:6.
```

The exact-period unit packet has size

```text
phi(14)=6,
```

which is exactly the unit seam `(Z/14Z)^* -> F_7^*` used in HYP-2626.

For the squarefree carrier:

```text
q=210=2*3*5*7, phi(210)=48.
```

The mask profile over `{2,3,5,7}` is

```text
{}:1
{2}:1, {3}:2, {5}:4, {7}:6
{2,3}:2, {2,5}:4, {2,7}:6, {3,5}:8, {3,7}:12, {5,7}:24
{2,3,5}:8, {2,3,7}:12, {2,5,7}:24, {3,5,7}:48
{2,3,5,7}:48.
```

These are exactly the products `prod_{p in mask}(p-1)`.

For the raw Hill `K_14` product:

```text
1260 = 2^2 * 3^2 * 5 * 7.
```

Compressing exact-period packets by the same squarefree masks gives selected
prime weights `p^a-1`, not merely `p-1`.  So the repeated `6,6` blocks in
HYP-2627 amplify the mask profile by

```text
2: (2^2-1)/(2-1) = 3
3: (3^2-1)/(3-1) = 4.
```

The full-mask copy mass in the raw profile is therefore

```text
(2^2-1)(3^2-1)(5-1)(7-1) = 3*8*4*6 = 576.
```

This is not `phi(1260)`, which is `288`; it is the mass of all exact-period
divisors whose reduced denominator touches every prime in the `{2,3,5,7}`
address.

## 1260/2520 Resonance

The THM-523 half-clash is

```text
15/36 - 2/5 - 1/70 - 1/504 = 1/2520.
```

The symmetry-doubled value is `1/1260`.  The totient-copy readout adds a new
exact-period explanation:

```text
phi(1260) = 288
phi(2520) = 576
full-mask mass of the raw 1260 profile = 576.
```

Thus the half-clash denominator `2520` sees exactly the all-four-prime packet
mass inside the raw Hill product `1260`.  After `tau <-> 1-tau` symmetry, the
measure denominator collapses back to the raw Hill product.

This is a stronger bridge than "1260 appears twice."  The four-prime packet
count before symmetry is the exact-period unit count at the half-gap
denominator.

## Interpretation

HYP-2628 separates three divisor notions that were easy to blur:

1. Raw divisor count, `tau(q)`.  Too coarse; it treats every divisor equally.
2. Exact top-period units, `phi(q)`.  Useful for unit seams such as
   `(Z/14Z)^*`.
3. Mask-compressed exact-period mass, `sum_{rad(d)=m} phi(d)`.  This is the
   correct carrier for raw denominators like `1260`, where repeated primes
   matter before projecting to the squarefree `210` coimage address.

The proposed proof transfer is therefore:

```text
raw denominator ledger
-> exact-period phi packets
-> squarefree mask compression
-> mod-7 unit/coimage quotient
-> signed repeated-root tail.
```

This says why HYP-2627's raw denominator discipline matters.  Dividing the
crossing value `1260/4=315` too early loses the dyadic exact-period layer; using
only `rad(1260)=210` too early loses the repeated `2,3` amplification.

## Tournament Analysis

The computation again uses proof quotients as vertices, not runners.

Hamiltonian path:

```text
totient_copy_operator
> exact_period_q_grid
> raw_1260_profile
> squarefree_210_profile
> unit_seam_phi14
> coimage_character_tail
> raw_divisor_counts
> raw_runner_vertices
```

Pairwise observable: preservation of the divisor-copy identity, the exact-period
LRC witness packet, squarefree transfer value, and compatibility with the
HYP-2626 coimage seam.

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles = 0
```

The script also builds mask-mass tournaments for `q=210` and `q=1260`.  The
`q=1260` top mask is the all-four-prime mask with mass `576`; for `q=210`, ties
put `{3,5,7}` just above `{2,3,5,7}` because the prime `2` contributes only
`p-1=1` in the squarefree profile.  This is a useful warning: squarefree
projection can make the dyadic coordinate invisible unless the raw denominator
ledger is retained first.

## Status

Open / structural.  The copy rule itself is theorem-level, but the LRC14
transfer is still a proof route.

Next tests:

1. Re-index the HYP-2626 repeated packets `(1,1,1,1,a,a)` and
   `(1,1,1,1,a,b)` by exact-period denominator packets rather than only residue
   classes.
2. Check whether the `576` all-mask packet controls the signed mass of the
   `12 -> 36` local champion after finite wall deletion.
3. Rewrite the HYP-2625 prime-mask recurrence as an explicit zeta/Mobius
   transfer matrix on exact-period packets, then attach the mod-7 coimage
   quotient as a functorial projection.
