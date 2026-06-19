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

## Namespace Claim

This stub reserves HYP-2628 / T876 for the user's divisor-copy reframe:

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

## Still Missing

This stub has not yet run the computation.  The next concrete task is to compute
the mask-weight profile

```text
copy_mass_N(M) = sum_{d|N, rad(d) has mask M} phi(d)
```

for `N` equal to `6`, `30`, `210`, `315`, `1260`, and the Hill products near
`K_14`; then compare that profile against the HYP-2625 prime-mask rows and the
HYP-2626 live `{7}` coimage transfer.

## Tournament Analysis Placeholder

Vertices to challenge before settling:

```text
totient_copy_masks, raw_squarefree_radicals, divisor_lattice_layers,
Hill_four_blocks, crossing_values, mod7_coimage_classes, raw_runner_vertices.
```

The expected preserved predicate is the LRC14 proof address after mod-210 and
coimage projection.  The expected destroyed information is exact witness time
geometry and individual divisor representatives.
