---
id: HYP-2672
title: LRC14 shell-full tail stratification after the new-speed dyadic block
status: OPEN; renumbered continuation after HYP-2671
source: codex-2026-06-20-S45
depends_on:
  - HYP-2671
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

# HYP-2672 - Shell-Full Tail Stratification

## Claim Being Tested

HYP-2671 localizes the post-shell-gate new-speed `1/3` constant to the dyadic
block row

```text
E'=(0,1,2,4,8,12,16,20), w=24,
Delta_w^+/p1 = 1371/4319.
```

The remaining shell-full constant work should now split around that canonical
exception:

```text
finite high pocket: max(E') <= 14,
dyadic block:       HYP-2671's m=4 row controls the new-speed 1/3 barrier,
intermediate band:  21 <= max(E') <= 24 has only finite >1/4 exceptions,
far-tail decay:     max(E') > 24 should satisfy Delta_w^+ <= p1(E')/4.
```

This continuation tests whether, after removing the HYP-2671 dyadic block, the
apparent "one open constant" becomes a finite apex-packet ledger plus a true
tail decay regime.  The intended exact scout will extend the HYP-2670/HYP-2671
shell-full quotient and stratify rows by `max(E')`, fold-target reciprocal
mass, Glaisher odd-carry profile, and phase packet concentration.

## Early Observation

A local exact re-read of the HYP-2670 `B=30` rows suggests a sharper picture
than the original three-band statement:

```text
rows above 3/10: only 3
rows above 1/3: only the B13 leader
rows with max(E') > 24 and above 1/4: 0
```

The only genuine new-speed row above `3/10` in that ledger is exactly the
HYP-2671 dyadic block

```text
E'=(0,1,2,4,8,12,16,20), w=24,
Delta_w^+/p1 = 1371/4319.
```

The working target is therefore:

1. certify the finite `max(E')<=14` high packet ledger;
2. use HYP-2671 for the dyadic block/new-speed `1/3` exception;
3. certify the finite `21..24` intermediate rows above `1/4`;
4. prove that the later tail has a true `p1/4` decay mechanism.

## Status

This is a renumbered continuation stub after the incoming HYP-2671/T910 result.
The next commit should attach the exact extension script, stored output, and a
proof-obligation refinement.  LRC(14) is not proved.
