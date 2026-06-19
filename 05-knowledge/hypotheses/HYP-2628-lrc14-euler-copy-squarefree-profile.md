---
id: HYP-2628
title: LRC(14) Euler-copy squarefree profile - totient multiplicities refine the mod-210 crossing carrier
status: OPEN
source: codex-2026-06-19-S21
depends_on:
  - HYP-2627
  - HYP-2626
  - HYP-2625
  - THM-538
related:
  - HYP-2624
  - HYP-2617
  - HYP-2619
  - THM-523
  - OPEN-Q-108
---

# HYP-2628 - LRC(14) Euler-Copy Squarefree Profile

## Core Identity

HYP-2628 / T876 records the user's divisor-copy reframe:

```text
Find copy counts c(n) >= 1 such that sum_{d|n} c(d) = n.
```

By Mobius inversion, the unique solution is

```text
c(n) = phi(n),
```

because `sum_{d|n} phi(d)=n`.  The intended LRC14 use is not to rediscover the
Euler totient, but to use the identity as a prime-mask copy ledger:

```text
N = sum_{d|N} phi(d)
  = sum_{squarefree masks M} copy_mass_N(M).
```

## Working Claim

The Euler-copy profile may refine HYP-2627's raw `K_14` crossing carrier:

```text
P_14 = 5*6*6*7 = 1260,   rad(P_14)=210.
```

Instead of recording only that the squarefree core is `{2,3,5,7}`, record how
many totient copies live on each squarefree mask.  This should distinguish:

1. the raw Hill product `1260`, which retains all four primes;
2. the divided crossing value `315`, which loses the dyadic copy gate;
3. the mod-6/mod-30/mod-210 address recurrence from HYP-2625;
4. the prime-7 coimage seam from HYP-2626.

## Computation

Script:

- `04-computation/lrc14_totient_copy_operator_codex_s21.py`
- output: `05-knowledge/results/lrc14_totient_copy_operator_codex_s21.out`

The script verifies the recursive copy rule against `phi(n)` for `n<=80`,
checks the exact-period residue census on `q=14,30,210`, computes compressed
prime-mask profiles for the LRC14 denominators, and records Tournament Analysis
for the proof quotients and mask-mass profiles.

## Exact-Period Interpretation

The copy rule is stronger than a static divisor identity.  On the denominator
`q` grid, residues `a/q` split by reduced denominator `d|q`; each packet has
size `phi(d)`.  Therefore the LRC squarefree profile should be read as

```text
raw denominator -> exact-period phi packets -> squarefree mask compression.
```

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

## Mask Profiles

For the squarefree carrier

```text
210 = 2*3*5*7,
```

mask weights are `prod_{p in mask}(p-1)`.  The full mask has mass

```text
phi(210)=48.
```

For the raw Hill `K_14` product

```text
1260 = 2^2 * 3^2 * 5 * 7,
```

compressing exact-period packets by the same squarefree masks gives selected
prime weights `p^a-1`, not merely `p-1`.  Thus the repeated `6,6` blocks in
HYP-2627 amplify the squarefree profile by

```text
2: (2^2-1)/(2-1) = 3
3: (3^2-1)/(3-1) = 4.
```

The full-mask copy mass in the raw profile is

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

## Exact-Period Safe-Center Transfer

The computation also applies the exact-period packet lens to the HYP-2625
modular-center predicate.  For each `Q`, it counts residues `a/Q` that are
strict-safe for every speed in `P subset {1,...,13}`.

The unweighted survivor rows record whether at least one exact-period packet is
present; the total rows retain the number of such packets.  In particular:

```text
Q=210 all residues, survivor counts by |P|:
1,13,78,286,715,1287,1716,1716,1287,715,286,78,11,0

Q=1260 all residues, survivor counts by |P|:
1,13,78,286,715,1287,1716,1716,1287,715,286,78,13,0
```

So `Q=1260` clears every `12`-of-`13` AP core, while still leaving the full AP
`{1,...,13}` strict-tight on this grid.

The drop-by-drop table is the local transversality diagnostic:

```text
drop  Q=210 all  Q=1260 all
1         14          88
2          6          38
3          8          48
4          4          30
5          8          40
6          0          10
7         18         104
8          8          56
9         12          72
10         4          30
11         8          74
12         0          16
13         4          38
```

Thus the radical/mod-210 grid misses exactly the AP drops `6` and `12`, while
the raw `1260` denominator recovers them.  This is a concrete reason to keep
the raw denominator ledger: the two troublesome AP one-drop cores are invisible
at `210` but visible at `1260`.

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

The computation uses proof quotients as vertices, not runners.

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
put `{3,5,7}` just above `{2,3,5,7}` because prime `2` contributes only
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
2. Explain why the `Q=210` misses are exactly AP drops `6` and `12`, and why
   the raw `Q=1260` packet repairs both without opening the full AP13 row.
3. Check whether the `576` all-mask packet controls the signed mass of the
   `12 -> 36` local champion after finite wall deletion.
4. Rewrite the HYP-2625 prime-mask recurrence as an explicit zeta/Mobius
   transfer matrix on exact-period packets, then attach the mod-7 coimage
   quotient as a functorial projection.
