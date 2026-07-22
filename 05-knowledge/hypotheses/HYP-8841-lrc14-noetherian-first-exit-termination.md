---
id: HYP-8841
title: "LRC14 Noetherian first-exit termination"
status: >
  OPEN / exact local no-go and bounded one-lift base cases available. THM-2050
  proves that complete period-14 top germs are globally blind. THM-2043 proves
  every audited one-lift alias has a labelled strict exit by q<=42, including
  an infinite Hasse-indistinguishable family. The target is a height-decreasing
  termination theorem for arbitrary multi-lift/covering rows, ending either at
  an AP/GW boundary atom or at a strict rational lonely phase.
source: codex-2026-07-21-DC2-LRC14-termination
related:
  - THM-523
  - THM-597
  - THM-2043
  - THM-2047
  - THM-2049
  - THM-2050
---

# HYP-8841 -- termination, not another local invariant

For a finite speed set `S`, define the first strict exit

```text
q_exit(S)=min{q>=2: exists a, gcd(a,q)=1,
                       min_(v in S)||av/q||>1/14},   (1)
```

with value infinity if the set is empty.

THM-2049 shows why this sidecar is necessary: `AP13` and `12->26` have the
same complete local germs at every unit point `a/14`, but

```text
q_exit(AP13)=infinity,       q_exit(12->26)<=12.     (2)
```

THM-2043 supplies a stronger finite base case: all 156 one-lift aliases in its
scope either are the AP/Goddyn--Wong boundary rows or have an exact labelled
exit by denominator `42`; its infinite `12->96+3444n` family retains the same
exit `(q,a)=(41,17)` at every Hasse depth.

The proposed LRC14 termination theorem is:

> Every primitive 13-speed row admits a finite sequence of labelled
> phase-height layer deletions/restrictions that strictly decreases a
> Noetherian height tuple and ends either at an AP/GW boundary certificate or
> at a strict exit (1).

A candidate height tuple is

```text
(unpaid q-witness count,
 first-exit denominator or infinity,
 maximum speed in the active magnitude fiber,
 active relation-layer rank,
 deleted-essential-point count).                   (3)
```

The order and even the correct coordinates in (3) are open.  What is fixed by
the DC(2) comparison is the proof architecture:

```text
local layer complex          usually solvable / acyclic
valuation or magnitude      records what localization forgot
strict height descent       prevents an infinite repair series
terminal atom               AP/GW boundary or explicit lonely phase
```

THM-2048 proves that local associated-graded solvability does not imply finite
polynomial closure.  For LRC14, THM-597 and THM-2049 give the corresponding
warning: local safe-component opening and complete top germs do not imply a
global lonely phase.  The missing theorem in both cases is termination.

The first concrete experiment is to compute the deletion/restriction height
tuple on the AP/GW, `12->26`, `12->36`, `12->96`, covering-min, and K33 banks,
then require every move to preserve endpoint owners and the exact LRC
predicate.  A scalar residue, collapse coefficient, or unlabelled toric layer
is an invalid vertex because it merges THM-2049's hostile pair.
