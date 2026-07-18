---
id: THM-989
title: Scale-twenty-two Hamming-six owner-feasibility deficit
status: CLAIMED FROM AN INDEPENDENT SCRATCH RECONSTRUCTION — all 984 scalar survivors have at most one owner-local feasible projection; the frozen primary certificate and an independent referee are in progress
source: codex-2026-07-17-S66 scale-twenty-two continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988]
related: [THM-978, THM-980, THM-981, HYP-6820]
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

The frozen primary is being written at
`04-computation/lrc13_scale_twenty_two_hamming_six_frontier_scout_codex_c22.py`.
Promotion requires a byte-stable primary, an independently structured replay,
and cross-mode checks.  The claim does not cover scale 24 or higher (scale 23
is prime and already excluded by THM-983), H5 ramification, non-AP/deep sheets,
or global sporadic emptiness.

