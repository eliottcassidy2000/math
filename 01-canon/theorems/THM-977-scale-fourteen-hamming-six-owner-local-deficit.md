---
id: THM-977
title: Scale-fourteen Hamming-six owner-local deficit
status: PROVED FINITE-EXACT — the complete 6,194,388,816-context primitive proper AP-centred common-scale-fourteen Hamming-six bank reduces to 576 scalar rows, and exact owner-local reachability eliminates all 3,456 owner rows with a two-sheet deficit or larger; independently reconstructed in C++ and Python
source: codex-2026-07-17-S66 scale-fourteen exact C++ certificate and independent Python scout
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-970, THM-974, THM-976, HYP-6820]
verification:
  - 04-computation/lrc13_scale_fourteen_hamming_six_owner_local_deficit_codex_c14.cpp
  - 05-knowledge/results/lrc13_scale_fourteen_hamming_six_owner_local_deficit_codex_c14.out
  - 04-computation/lrc13_scale_fourteen_hamming_six_frontier_scout_codex_c14.py
  - 05-knowledge/results/lrc13_scale_fourteen_hamming_six_frontier_scout_codex_c14.out
---

# THM-977 — scale fourteen dies before the global fibre

Scale thirteen is not a live face: THM-860 proves that
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
## Completeness and independent replay

The C++ certificate constructs the fourteen effective states by literal CRT
search.  It enumerates all 3,249 actual ordered divisor words and literally
traverses all 6,703,884 state words per support once, obtaining the independent
FNV64 audit `3f93d84053a6bdd3`.  It checks unit-independence of every scalar
mask cardinality before scanning all `924*3,249=3,002,076` labelled divisor
contexts.  For each of the 576 surviving contexts it runs exact set-union
reachability separately at all six owners.  No symmetry quotient is used in
the scalar or local verdict.  Optimized, unoptimized, and
AddressSanitizer/UndefinedBehaviorSanitizer builds reproduce the frozen output
byte-for-byte; strict warnings are clean.

The standard-library Python scout is algorithmically separate.  It derives
the CRT base algebraically, constructs masks as Python integers, enumerates
actual divisor tuples, and uses Python sets for owner-local reachability.  Its
normal and `python -O` runs reproduce a second frozen output byte-for-byte.
The implementations agree on every census, all three multiplication-orbit
representatives, and the complete maximum-union histogram.

Frozen SHA-256 values:

```text
C++ source     ad24fcf1e0d39f533859b98dff73ea92565b07fc867e456ec0a1abb2132feb83
C++ output     6bc0196c5eadf95814966a6027326f6d829d17198e7eb2385796a39a1c56eaf1
Python source  8b68f1615f330808d3dba2a96a0039942ca5d35899776e56424e7961797bef93
Python output  d4cbe8290f2a0e446145ea70f6ccbca15e44f0f6f9c5d8ecff4c4d38b2af7d1d
```

## Tournament and vertex audit

Tournament Analysis is honestly inapplicable after the terminal gate: no
context is feasible at even one owner, hence there is no global unit fibre,
no owner-obligation set, and no pair-intersection observable to orient.  The
faithful vertices are the fourteen sheet positions at one fixed owner, with
maximum reachable union size as the invariant.  Runner, provider, divisor,
residue, or owner vertices destroy the shared-unit incidence that produces
the deficit.  The three support multiplication orbits are useful algebraic
telemetry, but no orbit quotient enters the proof.

This closes only the primitive proper AP-centred common-scale-fourteen
Hamming-six face under THM-860's hereditary-lcm and common-sheet conclusions.
It does not close `c>=15`, H5 ramification, non-AP-centred/deep-sheet packets,
or global sporadic emptiness.  The next common-scale frontier is `c=15`.  ∎
