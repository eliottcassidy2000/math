---
id: THM-990
title: Scale-twenty-four Hamming-six owner-feasibility deficit
status: CLAIMED + FROZEN PRIMARY — all 66,984 scalar survivors have at most four owner-local feasible projections; the exact primary is byte-stable in normal and optimized modes, and an independently structured referee is still required
source: codex-2026-07-17-S66 scale-twenty-four continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988, THM-989]
related: [THM-978, THM-980, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_four_hamming_six_frontier_scout_codex_c24.py
  - 05-knowledge/results/lrc13_scale_twenty_four_hamming_six_frontier_scout_codex_c24.out
---

# THM-990 — scale twenty-four has at least two impossible owners

This namespace and companion computation filename reserve the next legal
primitive proper AP-centred common-scale Hamming-six face after scale 22;
scale 23 is uniformly prime-excluded by THM-983.

For `c=24`, the effective orders are
`1,2,3,4,6,8,12,24`, with twenty-four literal `(D,e)` states.  Hereditary
leave-one-out lcm enumeration gives 108,813 divisor words and 167,165,952
literal state words per support, hence

```text
924*167,165,952 = 154,461,339,648
```

unquotiented labelled state contexts.  A complete algebraic-CRT scratch
reconstruction leaves 66,984 scalar contexts on 854 supports across 202 order
multiplicity profiles.  Exact owner-local set-union reachability gives

```text
0 feasible owners: 64,962 contexts,
1 feasible owner :  1,800 contexts,
2 feasible owners:    192 contexts,
4 feasible owners:     30 contexts.
```

Thus every scalar row has at least two empty owner projections.  The 401,904
owner rows have maximum-union histogram

```text
12:72, 14:2136, 15:1644, 16:15876, 17:24420, 18:76296,
19:94872, 20:104592, 21:53040, 22:24948, 23:1704, 24:2304.
```

There are 20,302 distinct owner maximum vectors, and no reachable-mask bank
exceeded 7,728 states.  A global unit word would project to a feasible word at
all six owners, so the local deficit is terminal for this common-scale face.

## Frozen primary and carrier audit

The exact primary checks every algebraic CRT representative against literal
search, checks every local mask cardinality against an independent period
formula, audits the prime-power hereditary grammar against all six
leave-one-out lcms, and enumerates all support/order capacity rows without a
symmetry quotient.  Python normal and `-O` runs are byte-identical.  The
frozen artifact hashes are

```text
primary source  c3d20203dea9c36396db4a9975759b3d65b101c0aaf5a2c053d1370669861db1
primary output  3bb20b04576bcdb293ef474c4a083b43e8ab2e228ce87107b7bae36c02311ee4
```

Owner obligations are the terminal-faithful tournament vertices.  Their exact
ordered observable is `(feasible,max-union,capacity)`; lexicographic switching
with coordinate ties produces a transitive tournament in all 66,984 rows, with
score word `0,...,5`, no directed triangle, six singleton SCCs, and one
Hamiltonian path.  This is diagnostic telemetry: it forgets the absolute
24-sheet threshold, magnitudes, and exact unit-mask incidence.  Provider,
divisor, residue, isolated-sheet, and wall-event vertices lose shared-unit
incidence even earlier.

Promotion requires an independently structured replay.  This claim does not
cover scale 25 or higher, H5 ramification, non-AP/deep sheets, or global
sporadic emptiness.
