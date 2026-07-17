---
id: THM-978
title: Scale-fifteen Hamming-six owner-feasibility deficit
status: PROVED FINITE-EXACT — the complete 10,320,710,400-context primitive proper AP-centred common-scale-fifteen Hamming-six bank reduces to 2,184 scalar rows, and exact owner-local reachability proves that no row is feasible at more than four of its six owners; independently reconstructed in C++ and Python
source: codex-2026-07-17-S66 scale-fifteen exact C++ certificate and independent Python referee
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-974, THM-976, THM-977, HYP-6820]
verification:
  - 04-computation/lrc13_scale_fifteen_hamming_six_owner_local_deficit_codex_c15.cpp
  - 05-knowledge/results/lrc13_scale_fifteen_hamming_six_owner_local_deficit_codex_c15.out
  - 04-computation/lrc13_scale_fifteen_hamming_six_referee_codex_c15.py
  - 05-knowledge/results/lrc13_scale_fifteen_hamming_six_referee_codex_c15.out
---

# THM-978 — scale fifteen has at least two impossible owners

The primitive proper AP-centred common-scale-fifteen Hamming-six sheet bank is
empty.

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

The 66 closest rows are feasible at four owners.  They occur on eighteen
supports and have order-pattern census

```text
(#D1,#D3,#D5,#D15):
(0,0,2,4):12,   (0,4,0,2):18,   (0,4,1,1):24,   (0,4,2,0):12.
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

## Completeness and independent replay

The C++ certificate reconstructs all fifteen states by literal CRT search,
enumerates all 3,249 ordered divisor words, and literally traverses the
11,169,600 state words per support once, obtaining the independent FNV64 audit
`555244fbea49f335`.  It proves unit-independence of scalar cardinality before
checking all `924*3,249=3,002,076` labelled divisor contexts.  Each of the
2,184 scalar survivors then receives six exact set-union reachability DPs.
No multiplication quotient is used.  The support projection has 38
multiplication orbits of size twelve and one orbit of size six; these are
telemetry only.  Optimized, unoptimized, and
AddressSanitizer/UndefinedBehaviorSanitizer builds reproduce the frozen output
byte-for-byte under strict warnings.

The standard-library Python referee derives CRT bases algebraically,
reconstructs actual divisor words, verifies the provider/owner ratio reduction,
and uses Python sets for all owner-local DPs.  It freezes the complete scalar
bank, maximum-union vectors, and four-owner core with separate digests.  Normal
and `python -O` runs reproduce its output byte-for-byte.

Frozen SHA-256 values:

```text
C++ source     1ec49fd3bf2efab6ec68a4a3984b49dc68672e5434470e1ea618e9c0526112d9
C++ output     941d6be379be4d5f19a6e6d08099b259c250f94d40742870ee2891d955f88203
Python source  915c712443ff63c9cf18ac1467024863e35a0335c13299dfe71e2c0b6e3afb3b
Python output  97112efa22ca9a1e35ace13faaf38b666b93a420163db4acbfaea48dc717a12e
```

## Tournament and vertex audit

Take owners as vertices and compare the lexicographic pair
`(locally feasible, maximum union size)`; break ties along support-coordinate
order.  All 2,184 completed tournaments are transitive, with scores
`0,1,2,3,4,5`, no directed triangles, six singleton SCCs, and one Hamiltonian
path.  Their edge-flip counts range from zero to twelve.  This ranking cannot
prove the theorem: it discards both the absolute threshold fifteen and the
fact that at least two owners fail.  The faithful terminal carrier is the
six-owner feasibility subset, strengthened by the exact maximum-union vector;
provider vertices lose shared unit choices, while individual sheet vertices
are faithful only inside one labelled owner DP.

This closes only the primitive proper AP-centred common-scale-fifteen
Hamming-six face under THM-860's hereditary-lcm and common-sheet conclusions.
It does not close `c>=16`, H5 ramification, non-AP-centred/deep-sheet packets,
or global sporadic emptiness.  The next common-scale frontier is `c=16`.  ∎
