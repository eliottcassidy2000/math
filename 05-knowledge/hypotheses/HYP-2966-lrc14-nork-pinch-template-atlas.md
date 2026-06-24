---
id: HYP-2966
title: LRC14 NORK pinch-template atlas
status: FINITE ATLAS / NORK obstruction audit; not a proof
source: codex-2026-06-24-nork-pinch
related:
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2962
  - HYP-2961
  - HYP-2960
  - HYP-2956
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2951
  - HYP-2950
  - HYP-2949
  - HYP-2948
  - HYP-2947
  - HYP-2937
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_nork_pinch_template_audit_codex_20260624.py
  - 05-knowledge/results/lrc14_nork_pinch_template_audit_codex_20260624.out
external:
  - https://arxiv.org/abs/2606.22636
---

# HYP-2966: LRC14 NORK Pinch-Template Atlas

## Claim

The remaining LRC14 counterexample-shaped bucket should be treated as a NORK
bucket:

```text
NORK = No Open Residual Kernel
     = qdiv >= 14, no strict safe open interval,
       and not the AP/Goddyn-Wong equality boundary atom.
```

The proposed proof lemma is:

```text
Every qdiv>=14 AP-mutation packet either
  (1) is AP/GW boundary-only,
  (2) has a positive labelled pinch front,
  (3) carries a named petal/K33/state-lift label, or
  (4) is a genuine F6/F7 obstruction.
```

The computation in this session finds no case (4) in an enlarged audit bank.

## Computation

Script:

```text
04-computation/lrc14_nork_pinch_template_audit_codex_20260624.py
```

Stored output:

```text
05-knowledge/results/lrc14_nork_pinch_template_audit_codex_20260624.out
```

Default audit:

```text
one-swap add<=420
two-swap add<=60
three-swap add<=34
four-swap add<=24
workers=8
```

This extends the S150 packet-migration gauntlet in two ways:

1. three-swaps increase from `add<=30` to `add<=34`;
2. a first four-swap AP-neighborhood bank is added through `add<=24`.

Main counts:

```text
generated rows              705940
exact qdiv>=14 rows         141351
non-AP/GW F6 NORK rows      0

F1 AP/GW boundary           2
F2 positive unit-petal      28762
F3 positive K33             340
F4 positive q14 front       78651
F5 positive covering        33596
```

Per bank:

```text
0-swap: exact=1, boundary_only=1
1-swap add<=420: exact=2740, boundary_only=1, positive_open=2739
2-swap add<=60: exact=25884, positive_open=25884
3-swap add<=34: exact=73931, positive_open=73931
4-swap add<=24: exact=38795, positive_open=38795
```

So the only zero-open threshold-support rows in the audit are AP and GW.

## Pinch Templates

The new object is the **pinch template** of a positive-open row: the shortest
strict safe interval, with exact endpoints, active endpoint owners, width,
slack, q-class, C27-normalized owner labels, and packet atom labels.

This shifts the proof target from:

```text
safe mass is positive
```

to:

```text
which labelled endpoint-owner pinch forces positivity?
```

Recurring templates include:

```text
13L -> 12R    q14/off-apex front
14L -> 13R    covering front
11L -> 16R    q14/four-swap front
7L  -> 20R    unit-petal front
5L  -> 36R    K33 front
```

The tightest fronts in the large one- and two-swap scans can become very narrow
when a large added speed supplies one endpoint, but they remain positive and
labelled.  This is why scalar mass alone is a bad proof carrier: the owner
labels distinguish a harmless narrow front from an F6 zero-open kernel.

## Proof Route

HYP-2956 says LRC14 reduces to F0-F5, unless F6 or F7 appears.  HYP-2966 turns
that into a sharper route:

```text
F6 cannot be attacked by scalar lower bounds first.
F6 should be attacked by endpoint-owner pinch templates.
```

A plausible proof skeleton is:

1. qdiv gate removes `qdiv<14`.
2. AP/GW boundary-skeleton rigidity handles the only zero-open equality atoms.
3. Every local mutation away from AP/GW creates at least one labelled pinch
   template.
4. Petal/K33 templates route to the existing C27/unital/state-lift machinery.
5. Covering templates route to boundary-moment/gK8 positivity.
6. Any packet escaping those templates is, by definition, the next F6/F7
   obstruction to isolate.

This is not yet a proof of LRC.  It is progress because the negation has been
made highly local and labelled: a counterexample must now avoid both the
qdiv-gate and the observed endpoint-owner pinch atlas.

## Tournament Analysis

Vertices are proof carriers and pinch templates, not runners.

Pairwise observable:

```text
which carrier preserves NORK status, boundary owners, packet labels,
state-lift visibility, and anti-scalarization before discharge.
```

Switch:

```text
majority over the six retention coordinates; ties follow the listed order.
```

Fingerprint from the run:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles=0
Hamiltonian paths=1
SCCs all singleton
```

The tournament is deliberately transitive: first preserve F6/NORK status and
boundary skeletons, then pass through pinch templates and fixed-margin packets,
and only then scalarize.

## Missing Step

The next theorem target is:

```text
NORK pinch theorem:
For every primitive reduced LRC14 packet in the AP/GW source core, if the
packet is not AP/GW boundary-only, then some retained endpoint-owner pinch
template has positive width, unless the packet constructs HYP-2908/THM-572.
```

That theorem would close the F6 bucket for the local source core and would make
HYP-2956 substantially closer to an actual LRC14 proof.
