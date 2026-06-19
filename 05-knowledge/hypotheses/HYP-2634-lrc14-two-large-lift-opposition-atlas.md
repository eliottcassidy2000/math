---
id: HYP-2634
title: LRC(14) two-large lift-opposition atlas - bounded reciprocal signs are controlled by integer lift shape, not finite packet weight alone
status: OPEN
source: codex-2026-06-19-S25
depends_on:
  - HYP-2633
  - HYP-2632
  - HYP-2614
  - THM-538
related:
  - HYP-2619
  - OPEN-Q-108
---

# HYP-2634 - LRC(14) Two-Large Lift-Opposition Atlas

## Claim Stub

HYP-2633 found the first concrete lift-opposition warning: two `4+2` QR
finite packets with the same HYP-2632 weight `-25U` can have opposite bounded
reciprocal signs after integer lifting.  The example is:

```text
(1,2,8,9,15,22)   packet (1,1,1,1,2,2)   lift at H=16 positive
(1,4,8,11,15,22)  packet (1,1,1,1,4,4)   lift at H=16 negative
```

The next proof obligation is to decide whether this is an isolated accident or
a stable lift-shape phenomenon.  The proposed atlas should enumerate families
of integer representatives for the same finite two-large packet, compare their
bounded reciprocal signs, and record the low-height relation features that
drive sign opposition.

The working target is:

```text
finite packet class
+ integer lift profile
+ low-height relation-shell signature
=> predicted bounded reciprocal sign / cancellation pressure.
```

This is still a guardrail, not a proof of LRC(14).  Its purpose is to identify
the right discrepancy statistic for the HYP-2633 residue-lift
equidistribution / Abel-summation lemma.

## Status

Claim reserved for a stored lift-opposition scan over two-large representative
families.  The scan should explicitly consider alternate Tournament Analysis
vertices: runners, residues, lift offsets, relation-shell events, boundary
faces, additive-frequency shells, and proof obligations.
