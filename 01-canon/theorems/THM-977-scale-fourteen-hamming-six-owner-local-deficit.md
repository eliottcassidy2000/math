---
id: THM-977
title: Scale-fourteen Hamming-six owner-local deficit
status: CLAIMED — an exact exploratory census reduces 6,194,388,816 primitive proper AP-centred common-scale-fourteen contexts to 576 scalar rows and then eliminates every row at every owner; independent certificate hardening is in progress
source: codex-2026-07-17-S66 scale-fourteen frontier probe
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-970, THM-974, THM-976, HYP-6820]
---

# THM-977 — scale fourteen dies before the global fibre

This number is reserved for the scale-fourteen closure now under independent
certificate hardening.  Scale thirteen is not a live face: THM-860 proves that
`13` cannot divide the common scale of a primitive proper AP-centred packet,
because all retained and replacement speeds would then be divisible by
thirteen.  Thus `c=14` is the first legal common scale after THM-976.

For `c=14`, the effective orders are `1,2,7,14`, with fourteen `(D,e)`
states.  Exact leave-one-out lcm enumeration gives 3,249 hereditary order
words and 6,703,884 state words per support, hence

```text
924*6,703,884 = 6,194,388,816
```

labelled raw contexts.  Unit-independent scalar owner capacity leaves 576
contexts on 36 supports.  They form three multiplication orbits of twelve
supports; every support carries a unique pair of order-two providers and has
exactly sixteen scalar rows, obtained by assigning each of the other four
providers order seven or fourteen.  The order-pattern histogram is

```text
(#D1,#D2,#D7,#D14):
(0,2,0,4):  36,   (0,2,1,3): 144,   (0,2,2,2): 216,
(0,2,3,1): 144,   (0,2,4,0):  36.
```

Literal owner-local union-mask reachability eliminates all 576 rows at every
one of their six owners.  More sharply, the largest attainable union has size

```text
owner of order 7 or 14: at most 12 of 14 sheets,
owner of order 2:       at most 11 of 14 sheets
                         (and at most 10 in 96 owner rows).
```

Thus no global unit fibre remains to replay.  Here runner supports are only an
algebraic indexing device; the faithful obstruction vertices are the fourteen
local sheets at a fixed owner.  A tournament on runners or owners destroys the
cardinality deficit and is recorded only as deliberately lossy telemetry.
The exploratory certificate is
`lrc13_scale_fourteen_hamming_six_frontier_scout_codex_c14.py`.  A frozen
independent implementation, optimization/sanitizer replay, hashes, and an
exact finite formalization remain required before promotion to
`PROVED FINITE-EXACT`.
