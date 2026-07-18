---
id: THM-982
title: Scale-eighteen Hamming-six owner-feasibility deficit
status: PROVED STRUCTURAL + FINITE-EXACT — independent literal-search C++ and algebraic-CRT Python certificates exhaust the complete 27,490,799,952-context primitive proper AP-centred common-scale-eighteen Hamming-six bank, reduce it to 13,098 scalar rows on 684 supports, and prove that no row is owner-locally feasible at more than four of its six owners
source: codex-2026-07-17-S66 scale-eighteen exact probe
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-977, THM-978, THM-980, THM-981, HYP-6820]
---

# THM-982 — scale eighteen has at least two impossible owners

The primitive proper AP-centred common-scale-eighteen Hamming-six sheet bank is
empty before the global unit-fibre or metric-continuation stages.

For `c=18`, the effective orders are `1,2,3,6,9,18`, with eighteen literal
`(D,e)` states.  Hereditary leave-one-out lcm enumeration gives 26,961 divisor
words and 29,751,948 literal state words per support, hence

```text
924*29,751,948 = 27,490,799,952
```

unquotiented labelled state contexts.  Exact scalar
owner capacity leaves 13,098 rows on 684 supports across 63 order-multiplicity
patterns.  Exact owner-local set-union reachability gives

```text
0 feasible owners: 8,922 contexts,
1 feasible owner : 3,108 contexts,
2 feasible owners: 1,056 contexts,
4 feasible owners:    12 contexts.
```

Thus no row is locally feasible at all six owners; every row has at least two
empty owner projections.  Across all 78,588 owner rows, the exact
maximum-union histogram is

```text
10:48, 11:96, 12:2472, 13:6732, 14:21696,
15:19584, 16:16068, 17:6624, 18:5268.
```

The faithful certificate is the labelled six-owner
feasibility subset, strengthened by the maximum-union vector.  Ranking these
summaries produces tournament telemetry but loses the absolute eighteen-sheet
threshold.  Divisor or provider vertices are especially lossy at this scale:
the richer divisor lattice has 63 surviving multiplicity patterns, while the
terminal contradiction is an owner-labelled projection statement.

## Completeness and independent replay

The primary C++ certificate reconstructs every CRT mask by literal search,
enumerates all 26,961 hereditary divisor words, filters the 24,911,964 labelled
order contexts by unit-independent owner capacity, and computes exact reachable
unions for each of the six owner projections.  The Python referee instead
solves the CRT algebraically, generates the grammar with Cartesian products,
uses immutable-set reachability, and hashes every layer with SHA-256.  They
agree on all censuses, including the context-per-support histogram,
multiplication orbits `6:2,12:56`, 4,575 distinct maximum-union vectors, and
zero six-owner survivors.  The C++ `-O3`, `-O0`, and ASan/UBSan runs are
byte-identical; normal and optimized Python runs are byte-identical.

| artifact | SHA-256 |
|---|---|
| `04-computation/lrc13_scale_eighteen_hamming_six_owner_local_deficit_codex_c18.cpp` | `694562dfb5c79ef0dd7085e779712dbb032c96449648825888556455707a1b35` |
| `05-knowledge/results/lrc13_scale_eighteen_hamming_six_owner_local_deficit_codex_c18.out` | `9a1d26df03d4ce1d0730d28aeeead6f87e79410f6fa0ebcaac02d80f03114d82` |
| `04-computation/lrc13_scale_eighteen_hamming_six_referee_codex_c18.py` | `e9a52837ef35af0ae9f7afec1a208cca80dfffde9157587a4b012193ba278f2d` |
| `05-knowledge/results/lrc13_scale_eighteen_hamming_six_referee_codex_c18.out` | `9e2c2ea433c081ce0a56d14f6d436202f56da6823d5703309e1bd81d4c745268` |

For every scalar context, a global covariant unit word would project to a
locally feasible choice at every owner.  The exact maximum of four feasible
owners therefore supplies the terminal contradiction.  Local choices are
allowed to vary with the owner, so this implication uses the correct
direction: local failure is decisive, while six local successes would only be
necessary, not sufficient.

This theorem concerns only the primitive proper AP-centred common-scale-eighteen
Hamming-six face under THM-860.  It does not close `c>=19`, the H5 bank,
non-AP-centred/deep-sheet continuations, or global sporadic emptiness.
