---
id: HYP-2671
title: LRC14 shell-full constant stratification around the one open p1-tax constant
status: OPEN; namespace claimed, exact extension in progress
source: codex-2026-06-20-S45
depends_on:
  - HYP-2670
  - HYP-2669
  - HYP-2668
  - HYP-2666
  - HYP-2661
related:
  - THM-545
  - HYP-2665
  - HYP-2664
  - HYP-2662
  - HYP-2643
  - OPEN-Q-108
---

# HYP-2671 - Shell-Full Constant Stratification

## Claim Being Tested

HYP-2670 split the shell-1-full p1-tax obligation into three constants:

```text
finite high pocket: max(E') <= 14,
new-speed decay:    max(E') > 14 should satisfy Delta_w^+ <= p1(E')/3,
far-tail decay:     max(E') > 24 should satisfy Delta_w^+ <= p1(E')/4.
```

This session tests whether the apparent "one open constant" is really a small
finite apex-packet ledger plus a monotone decay regime.  The intended exact
scout will extend the HYP-2670 shell-full quotient and stratify rows by
`max(E')`, fold-target reciprocal mass, Glaisher odd-carry profile, and phase
packet concentration.

## Early Observation

A local exact re-read of the HYP-2670 `B=30` rows suggests a sharper picture
than the original three-band statement:

```text
rows above 3/10: only 3
rows above 1/3: only the B13 leader
rows with max(E') > 24 and above 1/4: 0
```

The only genuine new-speed row above `3/10` in that ledger is the dyadic apex
avatar

```text
E'=(0,1,2,4,8,12,16,20), w=24,
Delta_w^+/p1 = 1371/4319.
```

The working target is therefore not just "prove `2p1/5` after the shell gate".
It is:

1. certify the finite `max(E')<=14` high packet ledger;
2. isolate the dyadic apex avatar as the only new-speed row above `3/10`;
3. prove that the later tail has a true `p1/4` decay mechanism.

## Status

This is a claimed namespace stub.  The next commit should attach the exact
script, stored output, and a proof-obligation refinement.  LRC(14) is not
proved.
