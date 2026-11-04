# Independent review of the order-five quartic classification

**Verdict: PASS. No discrepancy found.**

This review read `response/CLASSIFICATION.md` and `classify_five.py`, replayed the producer in normal and optimized Python, and compared its output byte-for-byte with the supplied `classification.json`.

## Independent construction

`independent.py` does not import the producer or its cover normal form. It enumerates literal Hamilton words on the five named old vertices and two ports. A quartic word must have both ports internal and nonadjacent. There are 1,440 such word skeletons; each tests exactly two old arcs and contributes +1 to its named bidegree-(2,2) coordinate. The independent checker accumulates these literal contributions for all 1,024 tournament codes.

All 102,400 named coordinates, including zeros, are compared as exact integer vectors with the published records. Grouping uses full tuples, not digests or sorted coefficient multisets. Adjacency rows and the common physical-name basis are checked for every code.

## Classification results

The independent count gives exactly 502 classes: 492 two-member whole-converse pairs and ten four-member classes. For each four-member class it explicitly verifies the same named triangle, its two orientations, and both possible source/sink orders on the exterior names. All ten named triangle sets occur exactly once.

A separate Floyd-style reachability computation gives 544 strong tournaments. Each belongs to a two-member converse class. Thus the strong-tournament conclusion is justified at order five only; no all-order reconstruction conclusion follows.

## The exceptional pair is not a full-response collision

For T=(10,0,3,6,15) and U=(6,0,10,3,15), direct enumeration of all permutations of the seven extended names gives 31 and 29 Hamilton paths at the same query a_0=1,c_2=1, with all other attachment bits zero and x->y fixed.

Independent full Boolean Hamilton-word expansion yields exactly six nonzero quadratic and 36 nonzero cubic difference coefficients, with no constant, linear, or quartic difference. Both old Hamilton counts are 3. This verifies the claimed example, its physical names, and the distinction between top-degree and full-response equality.

## Controls and artifacts

Normal and optimized producer outputs agree. Normal and optimized independent outputs also agree. The independent receipt retains every difference coefficient for the example and all ten exceptional classes. All checks remain active under optimization.

Artifacts: `classification.normal.json`, `classification.optimized.json`, `independent.py`, `independent.json`, and `independent.optimized.json`.
