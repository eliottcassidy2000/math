---
id: THM-989
title: Scale-twenty-two Hamming-six owner-feasibility deficit
status: CLAIMED + FROZEN PRIMARY — all 984 scalar survivors have at most one owner-local feasible projection; the standard-library primary is byte-stable in normal and optimized modes, and an independently structured referee is still required
source: codex-2026-07-17-S66 scale-twenty-two continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988]
related: [THM-978, THM-980, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_two_hamming_six_frontier_scout_codex_c22.py
  - 05-knowledge/results/lrc13_scale_twenty_two_hamming_six_frontier_scout_codex_c22.out
---

# THM-989 — scale twenty-two has at least five impossible owners

This namespace and the companion computation filename are reserved for the
next exact primitive proper AP-centred common-scale Hamming-six face.

For `c=22`, the effective orders are `1,2,11,22`, with twenty-two literal
`(D,e)` states.  Hereditary leave-one-out lcm enumeration gives 3,249 divisor
words and 100,975,500 literal state words per support, hence

```text
924*100,975,500 = 93,301,362,000
```

unquotiented labelled state contexts.  A from-scratch algebraic-CRT probe
leaves 984 scalar contexts on 180 supports across eight multiplicity profiles.
Exact owner-local set-union reachability gives

```text
0 feasible owners: 792 contexts,
1 feasible owner : 192 contexts.
```

The 5,904 owner rows have maximum-union histogram

```text
16:864, 17:1584, 18:2784, 19:480, 22:192.
```

Thus every candidate has at least five empty owner projections.  This is
terminal for the common-scale Hamming-six face: a global unit word would
project to a feasible word at every owner.

## Frozen primary and carrier audit

The standard-library Python primary solves CRT representatives algebraically,
cross-checks every representative by literal search, and checks every mask
cardinality against an independent one-period formula.  It traverses all
3,002,076 support/order rows without quotienting, hashes every scalar survivor
and every sorted owner-local reachable-mask bank, and reproduces the following
stronger telemetry:

```text
scalar support orbits under multiplication: 15, all of size 12;
distinct capacity vectors: 484;
distinct owner max-union vectors: 127;
reachable-union-bank SHA-256:
baf8aa9ee67d7686b25e8665bec8f94514d7abb6e7be780bebc2f98039675f1b.
```

Normal and `python -O` runs reproduce the frozen output byte-for-byte.  Frozen
artifact SHA-256 values are

```text
primary source  bd521cea444a19cae759e8ac2b4251bc400f481cbf957877b53e2a622444b592
primary output  40cda9290d98fcbb819cec3e63d7cfe16014c266397bb38b8f851fd4caa81a50
```

Owner obligations are the terminal-faithful vertices: their exact summary is
`(feasible,max-union,capacity)`.  Lexicographic comparison with coordinate tie
breaks gives a transitive tournament in all 984 rows (score word `0,...,5`,
zero directed triangles, six SCCs, one Hamiltonian path).  That tournament is
telemetry only; it forgets the absolute 22-sheet threshold and the number of
empty owners.  Provider, divisor, residue, isolated-sheet, and wall-event
vertices lose shared-unit incidence earlier still.

Promotion requires an independently structured replay.  The claim does not
cover scale 24 or higher (scale 23 is prime and already excluded by THM-983),
H5 ramification, non-AP/deep sheets, or global sporadic emptiness.
