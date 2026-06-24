---
id: HYP-2953
title: LRC14 source-spectrum pullback
status: CLAIMED / synthesis proof target; evidence and proof obligations still missing
source: codex-2026-06-24-S149
related:
  - HYP-2952
  - HYP-2951
  - HYP-2950
  - HYP-2947
  - HYP-2936
  - HYP-2929
  - HYP-2928
  - HYP-2927
  - HYP-2920
  - HYP-2917
  - HYP-2486
  - THM-523
  - THM-566
  - THM-572
  - THM-420
  - THM-430
  - OPEN-Q-108
---

# HYP-2953: LRC14 Source-Spectrum Pullback

This placeholder reserves the S149 synthesis target.

The proposed missing object is a **source-spectrum pullback**:

```text
Farey/Stern-Brocot binding node
  + threshold observer-source lift
  + Haar/Baire boundary-vs-interior code
  + packet labels retained until discharge.
```

The claim is not a proof.  The working idea is that the old proof attempts
have been describing different projections of one object:

```text
tournament spectrum      remembers magnitude by moving through phase
source-cone lift         remembers the LRC witness predicate
q/Farey branch           remembers the binding scale
Haar/Baire carrier       remembers boundary-only versus open witness
C27/unital/K33 packets   remember which boundary debt is dischargeable
state lift               gives the forbidden-H=7 endpoint
```

S149 will test whether this pullback gives the right theorem shape for both
irreducible branches:

```text
non-covering tight branch: deepest source-spectrum node stays at apex only for
AP/GW;

covering branch: no fixed denominator source works, but a positive Haar
source interval or a forbidden state-lift packet must appear.
```

Required before promotion: record the history links, state the preserved and
destroyed predicates for each quotient, run Tournament Analysis over proof
carriers, and post the synthesis attempt to `poke-forum/`.
