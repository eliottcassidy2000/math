---
id: THM-988
title: Scale-twenty-one Hamming-six owner-feasibility deficit
status: CLAIMED — a frozen exact C++ certificate exhausts the 77,810,217,408-context primitive proper AP-centred common-scale-twenty-one Hamming-six bank and finds no scalar row feasible at all six owners; an independently derived algebraic-CRT referee is in progress
source: codex-2026-07-17-S66 scale-twenty-one continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-982, THM-983, THM-986]
related: [THM-978, THM-980, THM-981, HYP-6820]
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
owner-local reachability gives the provisional census

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

The C++ source/output SHA-256 values are respectively
`08772ed1a3e5acf59eb1f33ac0420de6394fa6a7648b8c6efcdb79a3d36268bc`
and `65a55c93aaea9d53ee98ff2d3b36fa61afcc4511eecf90904d621a415e4fb284`;
optimized, unoptimized, and ASan/UBSan outputs agree byte-for-byte.

The faithful carrier is the labelled owner capacity/feasibility/max-union
vector.  Its completed tournament is transitive in every scalar survivor and
forgets both the absolute twenty-one-sheet threshold and the number of empty
owners.  Provider, divisor, residue, isolated-sheet, and wall-event vertices
lose shared-unit incidence.

Promotion requires the in-progress independent algebraic-CRT Python referee.
This claim concerns only the primitive proper AP-centred common-scale-twenty-one
Hamming-six face; it does not close scale 22 or higher, the H5 bank,
non-AP/deep-sheet continuations, or global sporadic emptiness.
