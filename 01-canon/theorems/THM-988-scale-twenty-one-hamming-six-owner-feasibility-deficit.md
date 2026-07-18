---
id: THM-988
title: Scale-twenty-one Hamming-six owner-feasibility deficit
status: PROVED STRUCTURAL + FINITE-EXACT — independent literal-search C++ and algebraic-CRT Python certificates exhaust the 77,810,217,408-context primitive proper AP-centred common-scale-twenty-one Hamming-six bank and prove that no scalar row is owner-locally feasible at more than four owners
source: codex-2026-07-17-S66 scale-twenty-one continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-982, THM-983, THM-986]
related: [THM-978, THM-980, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_one_hamming_six_frontier_scout_codex_c21.cpp
  - 05-knowledge/results/lrc13_scale_twenty_one_hamming_six_frontier_scout_codex_c21.out
  - 04-computation/lrc13_scale_twenty_one_hamming_six_referee_codex_c21.py
  - 05-knowledge/results/lrc13_scale_twenty_one_hamming_six_referee_codex_c21.out
---

# THM-988 — scale twenty-one has at least two impossible owners

The primitive proper AP-centred common-scale-twenty-one Hamming-six sheet bank
is empty before the global covariant unit-fibre or metric-continuation stages.

For `c=21`, the effective orders are `1,3,7,21`, with twenty-one literal
`(D,e)` states.  Hereditary leave-one-out lcm enumeration gives 3,249 divisor
words and 84,210,192 literal state words per support, hence

```text
924*84,210,192 = 77,810,217,408
```

unquotiented labelled state contexts.  Both certificates leave 350 scalar
contexts on 80 supports across nine multiplicity profiles.  Exact owner-local
reachability gives the census

```text
0 feasible owners: 182 contexts,
1 feasible owner :  96 contexts,
4 feasible owners:  72 contexts.
```

Thus every scalar row has at least two empty owner projections.  The 2,100
owner rows have maximum-union histogram

```text
17:1056, 18:648, 20:12, 21:384.
```

The two all-order-twenty-one scalar rows are especially revealing.  All six
owners are exactly scalar-capacity tight, yet every owner-local union has
maximum twenty.  This is a pure overlap/projection obstruction, invisible to
scalar slack: every owner owes exactly one sheet, but the owed sheet is visible
only before cardinality compresses the mask family.

## Independent replay

The primary certificate constructs CRT representatives by literal search,
checks each sheet cardinality against an independent period formula, and uses
a packed exact union-mask table.  The Python referee was derived before reading
the primary.  It instead solves CRT algebraically, enumerates Cartesian grammar
fibres, stores each reachable owner bank as an immutable set, and hashes every
mask in increasing order.  After freezing its answer, it replayed the C++ FNV
serialization and matched the primary's mask, order, multiplicity, scalar-bank,
capacity-vector, row-capacity, and owner-profile digests exactly.

The implementations agree on all decisive counts: the nine multiplicity
profiles, support-context histogram `1:14,2:48,8:12,24:6`, multiplication
orbits `2:1,6:1,12:6`, 131 capacity vectors, 34 owner maximum vectors, all
owner histograms above, and zero all-six contexts.

| artifact | SHA-256 |
|---|---|
| `04-computation/lrc13_scale_twenty_one_hamming_six_frontier_scout_codex_c21.cpp` | `08772ed1a3e5acf59eb1f33ac0420de6394fa6a7648b8c6efcdb79a3d36268bc` |
| `05-knowledge/results/lrc13_scale_twenty_one_hamming_six_frontier_scout_codex_c21.out` | `65a55c93aaea9d53ee98ff2d3b36fa61afcc4511eecf90904d621a415e4fb284` |
| `04-computation/lrc13_scale_twenty_one_hamming_six_referee_codex_c21.py` | `14a5d181a93befc5516711a3340b1ceaa4899d3e478d649cb4504a7e3b027ae7` |
| `05-knowledge/results/lrc13_scale_twenty_one_hamming_six_referee_codex_c21.out` | `fbf1680761bf9d42e2b95e4ca4f37e239a7e9840b40fad46816870c641b1987b` |

Optimized, unoptimized, and ASan/UBSan C++ outputs agree byte-for-byte.  Normal
and `python -O` referee outputs also agree byte-for-byte.

## Faithful carrier and terminal implication

The faithful carrier is the labelled owner capacity/feasibility/max-union
vector.  Its completed tournament is transitive in every scalar survivor and
forgets both the absolute twenty-one-sheet threshold and the number of empty
owners.  Provider, divisor, residue, isolated-sheet, and wall-event vertices
lose shared-unit incidence.

For every context a global unit word would induce a feasible choice in all six
owner projections, so the exact maximum of four is terminal.  The converse is
not used: separately feasible owner witnesses need not glue to one global
word.  This theorem concerns only the primitive proper AP-centred
common-scale-twenty-one Hamming-six face; it does not close scale 22 or higher,
the H5 bank, non-AP/deep-sheet continuations, or global sporadic emptiness.
