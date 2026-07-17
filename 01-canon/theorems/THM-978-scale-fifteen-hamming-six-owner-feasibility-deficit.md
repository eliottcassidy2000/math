---
id: THM-978
title: Scale-fifteen Hamming-six owner-feasibility deficit
status: CLAIMED — an exact Python census reduces 10,320,710,400 primitive proper AP-centred common-scale-fifteen contexts to 2,184 scalar rows, and no row is locally feasible at more than four of its six owners; independent C++ hardening is in progress
source: codex-2026-07-17-S66 scale-fifteen frontier probe
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-974, THM-976, THM-977, HYP-6820]
---

# THM-978 — scale fifteen has at least two impossible owners

This number is reserved for the scale-fifteen closure now under independent
certificate hardening.

For `c=15`, the effective orders are `1,3,5,15`, with fifteen `(D,e)` states.
Exact leave-one-out lcm enumeration gives 3,249 hereditary order words and
11,169,600 state words per support, hence

```text
924*11,169,600 = 10,320,710,400
```

labelled raw contexts.  Unit-independent scalar owner capacity leaves 2,184
contexts on 462 supports across sixteen order-multiplicity patterns.  Exact
owner-local union-mask reachability gives the following number of feasible
owners per scalar context:

```text
0 owners: 750 contexts,
1 owner : 456 contexts,
2 owners: 912 contexts,
4 owners:  66 contexts.
```

In particular no context is feasible at all six owners.  Across all 13,104
owner rows, the largest reachable union has histogram

```text
11 sheets:  804,   12 sheets: 7,512,   13 sheets: 1,812,
14 sheets:  432,   15 sheets: 2,544.
```

The faithful terminal datum is the six-owner feasibility subset, strengthened
by the exact maximum-union vector.  For tournament telemetry, orient two owners
by decreasing maximum union and break ties by label order.  The completion is
always transitive, hence it forgets both the absolute threshold fifteen and
the fact that at least two owners fail.  Runner, provider, divisor, residue,
and unlabelled sheet vertices lose still more of the shared-unit incidence.

The independent Python certificate is
`lrc13_scale_fifteen_hamming_six_referee_codex_c15.py`.  A frozen C++
reconstruction, optimized/unoptimized/sanitizer replay, frozen Python output,
and hashes remain required before promotion to `PROVED FINITE-EXACT`.
