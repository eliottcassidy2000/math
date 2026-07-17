---
id: THM-980
title: Scale-sixteen Hamming-six owner-feasibility deficit
status: PROVED FINITE-EXACT — the complete 13,806,600,192-context primitive proper AP-centred common-scale-sixteen Hamming-six bank reduces to 2,540 scalar rows, and exact owner-local reachability proves that no row is feasible at more than two of its six owners; independently reconstructed in C++ and Python
source: codex-2026-07-17-S66 scale-sixteen exact C++ certificate and independent Python referee
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-976, THM-977, THM-978, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_sixteen_hamming_six_frontier_scout_codex_c16.cpp
  - 05-knowledge/results/lrc13_scale_sixteen_hamming_six_frontier_scout_codex_c16.out
  - 04-computation/lrc13_scale_sixteen_hamming_six_referee_codex_c16.py
  - 05-knowledge/results/lrc13_scale_sixteen_hamming_six_referee_codex_c16.out
---

# THM-980 — scale sixteen has at least four impossible owners

The primitive proper AP-centred common-scale-sixteen Hamming-six sheet bank is
empty.

For `c=16`, the effective orders are `1,2,4,8,16`, with sixteen literal
`(D,e)` states.  Leave-one-out lcm enumeration gives 5,385
hereditary order words and 14,942,208 literal state words per support, hence

```text
924*14,942,208 = 13,806,600,192
```

labelled raw contexts.  Unit-independent scalar owner capacity leaves 2,540
contexts on 404 supports across twenty order-multiplicity patterns.  Exact
owner-local union-mask reachability
gives

```text
0 feasible owners: 2,006 contexts,
1 feasible owner :   432 contexts,
2 feasible owners:   102 contexts.
```

Thus no context is feasible at all six owners;
indeed every row has at least four impossible owners.  Across all 15,240 owner
rows its maximum reachable-union histogram is

```text
 9 sheets:   144,  10 sheets:   468,  11 sheets:   876,
12 sheets: 2,316,  13 sheets: 4,068,  14 sheets: 4,740,
15 sheets: 1,992,  16 sheets:   636.
```

The faithful terminal certificate is the ordered six-owner vector of local
maximum unions (or its thresholded feasibility subset), before any shared-unit
global replay.  This is a **pre-nerve** obstruction: because at least one owner
obligation set is already empty, the nerve of globally covariant unit words
never needs to be formed.  A tournament that ranks owners by
`(feasible, maximum union)` is useful telemetry but loses the absolute
sixteen-sheet threshold unless the full pair observable is retained.

## Completeness and independent replay

The C++ certificate reconstructs every CRT base by literal search, verifies
the complete divisor/unit grammar, traverses all 14,942,208 literal state words
per support, and checks all `924*5,385=4,975,740` labelled divisor contexts.
Its scalar projection has multiplication-orbit histogram

```text
size 2: 1 orbit,   size 6: 1 orbit,   size 12: 33 orbits,
```

but no multiplication quotient is used in the proof.  Each scalar survivor
receives six exact set-union reachability DPs.  Since the all-six local bank is
empty, the subsequent globally covariant replay is vacuously empty.
Optimized, unoptimized, and AddressSanitizer/UndefinedBehaviorSanitizer builds
reproduce the frozen output byte-for-byte under strict warnings.

The standard-library Python referee is independently structured.  It derives
CRT bases algebraically, rebuilds the literal mask table, enumerates the full
state grammar with a separate SHA-256 digest, packs scalar capacities into six
independent integer fields, and uses Python sets for the owner-local DPs.  It
reproduces the twenty multiplicity patterns, support census, orbit census,
owner-feasibility histogram, maximum-union histogram, and all 1,210 distinct
ordered owner-maximum vectors.  Normal and `python -O` runs reproduce its
frozen output byte-for-byte.

Frozen SHA-256 values:

```text
C++ source     97a841bd2841a8dadaf932c813755fb50efb4235c08186f61322a4c8ef409ec9
C++ output     68d3a09379dd7dd37e74134c34a1f8ae4954d4e4af1fd2c075aa8e3b898e4552
Python source  9598ac0060e01bfe6893f21345630abb6465944014ccac1542dfc587985ec47b
Python output  c6ef771a288aac091176f43d70a43906a47f5132de6dd53054b266c10835d54c
```

## Tournament and vertex audit

Take owners as vertices and retain the exact ordered pair of summaries
`(locally feasible, maximum union)` as the pair observable.  Orient by its
lexicographic sign and use support-coordinate order as the tie Hamiltonian
path.  All 2,540 completions are transitive: their score multiset is
`{0,1,2,3,4,5}`, they have no directed triangle, six singleton SCCs, and one
Hamiltonian path.  The frozen outputs record the full tie-edge and edge-flip
histograms.  The oriented tournament alone is not a proof: it forgets the
absolute threshold sixteen and the deficit magnitudes.  Provider, residue,
divisor-word, and isolated sheet vertices lose still more of the labelled
six-owner projection data.

This theorem concerns only the primitive proper AP-centred common-scale-sixteen
Hamming-six face under THM-860.  It does not close `c>=17`, the H5 bank,
non-AP-centred/deep-sheet continuations, or global sporadic emptiness.  The
scale-seventeen continuation is recorded separately in THM-981.  ∎
