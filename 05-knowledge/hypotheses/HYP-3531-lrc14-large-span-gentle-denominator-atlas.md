---
id: HYP-3531
title: LRC14 large-span gentle / denominator atlas
status: RESERVED / extends HYP-3530; evidence pending exact larger-span and structured-family audit
source: codex-2026-06-29 continuation of HYP-3530 multi-denominator rho_D and k<12 attacker-floor scout
tangent: T1531
technique: LTI-531
tournament_technique: LTT-431
script: 04-computation/lrc14_large_span_gentle_denominator_atlas_codex_20260629.py
result: 05-knowledge/results/lrc14_large_span_gentle_denominator_atlas_codex_20260629.out
reflection: 07-reflections/lrc14-large-span-gentle-denominator-atlas-codex-20260629.md
related:
  - HYP-3530
  - THM-530
  - THM-531
  - THM-578
  - HYP-3529
  - OPEN-Q-108
---

# HYP-3531: LRC14 Large-Span Gentle / Denominator Atlas

## Claim Reservation

HYP-3530 makes the next gap precise: below `k=12`, the bounded-core bank
through `span <= k+5` has no attackers, no below-consecutive shapes, and
consecutive minimizers.  HYP-3531 reserves the continuation lane:

```text
1. extend exact primitive attacker scans to larger spans;
2. stress structured families that could evade the bounded bank;
3. test whether the same gentle rows carry small or patterned rho_D witnesses.
```

The intended preserved predicate is the THM-530 union-bound floor:

```text
mu_1/7(E) >= thr_k := 1 - min_{|P|=13-k} meas(G_P)
```

for `k=8,9,10,11`.  A row is an attacker if it violates this inequality.
The denominator side remains a certificate bank: `max_D rho_D` may route
finite rational witnesses, but it is not by itself a positive-measure theorem.

## Assumption Challenge

Tournament vertices should be proof carriers and row families: exact bounded
span banks, one-tail rows, two-tail rows, perforated APs, even-lattice bridge
rows, sparse random probes, and denominator families.  Runners, raw spans, and
single denominator counts are shadows unless they preserve the attacker
predicate plus the witness sidecar.
