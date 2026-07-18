---
id: THM-988
title: Scale-twenty-one Hamming-six owner-feasibility deficit
status: PROVED STRUCTURAL + REFEREED FINITE-EXACT — independent literal-search C++ and algebraic-CRT Python certificates exhaust the 77,810,217,408-context primitive proper AP-centred common-scale-twenty-one Hamming-six bank and prove that no scalar row is owner-locally feasible at more than four owners
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

This namespace is reserved for the exact scale-twenty-one continuation of the
primitive proper AP-centred common-scale Hamming-six classification.  THM-987
is left available for the already named canonical exact-deep-count target.

For `c=21`, the effective orders are `1,3,7,21`, with twenty-one literal
`(D,e)` states.  Hereditary leave-one-out lcm enumeration gives 3,249 divisor
words and 84,210,192 literal state words per support, hence

```text
924*84,210,192 = 77,810,217,408
```

unquotiented labelled state contexts.  The primary certificate leaves 350
scalar contexts on 80 supports across nine multiplicity profiles.  Exact
owner-local reachability gives the exact census

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

The two all-order-twenty-one scalar rows are especially revealing: all six
owners are exactly scalar-capacity tight, yet every owner-local union has
maximum twenty.  This is a pure overlap/projection obstruction, invisible to
scalar slack.

## Independent replay

The primary constructs every CRT representative by literal search and carries
reachable unions in an epoch-marked packed table.  The independent Python
referee instead solves the CRT algebraically, verifies every mask cardinality
against a separate period-count formula, constructs the hereditary grammar by
Cartesian products plus the prime-provider criterion, and carries owner-local
banks as immutable sets.  It reproduces the primary mask, order, multiplicity,
scalar-bank, capacity, and owner-profile FNV checkpoints while also hashing
every reachable mask in sorted order.

Across the 2,100 owner rows the referee hashes 1,393,896 reachable masks.  Its
full reachable-bank SHA-256 is
`c3e7ce9c46dffd5846bb663f9f490ac9bf3fae76885a7085d5c3172a2907a8aa`.
The two implementations also agree on the nine multiplicity profiles, all 350
labelled scalar survivors, 131 capacity vectors, support-context histogram,
multiplication-orbit census, 34 owner maximum vectors, slack histograms, and
the complete owner-feasibility census.

The C++ source/output SHA-256 values are respectively
`08772ed1a3e5acf59eb1f33ac0420de6394fa6a7648b8c6efcdb79a3d36268bc`
and `65a55c93aaea9d53ee98ff2d3b36fa61afcc4511eecf90904d621a415e4fb284`;
optimized, unoptimized, and ASan/UBSan outputs agree byte-for-byte.

The Python source/output SHA-256 values are respectively
`9bc5e99523f7601d98d23e9090f4351791033ae58f945ed92a95f2d04be605f7`
and `4ee08c57cf70d248c46f60cc4db10c30d342878a91d03ea72e21c7f22b6a76c6`;
normal and `python -O` outputs agree byte-for-byte with the frozen result.

The faithful carrier is the labelled six-tuple of exact owner-local mask banks;
it preserves sheet identities and each full-mask predicate.  The completed
owner-summary tournament uses the exact pair observable
`(feasible, maximum union, scalar capacity)` and is transitive in every scalar
survivor, with one Hamiltonian path.  It forgets the absolute twenty-one-sheet
threshold, exact masks, unit witnesses, and simultaneous cross-owner
compatibility.  Runner/provider, gap, fixed-section, section-boundary,
wall-event, residue, cover-arc, Fourier-mode, and matroid-circuit vertices all
lose the labelled shared-unit incidence needed by the terminal predicate.

This claim concerns only the primitive proper AP-centred common-scale-twenty-one
Hamming-six face; it does not close scale 22 or higher, the H5 bank,
non-AP/deep-sheet continuations, or global sporadic emptiness.
