---
id: THM-825
title: Quadratic positional moments are a literal mirror-layer codec through size fifteen
status: RESERVED / PROOF AND EXACT REPLAY IN PROGRESS
source: codex-2026-07-15-S13-postjoin
depends_on: [THM-553, THM-809, THM-818]
related: [THM-796, THM-814, HYP-6880]
---

# THM-825 — quadratic positional moments close the layer codec through `n=15`

Namespace reservation after a live-main check.  THM-824 was claimed
concurrently by the fixed-pair symmetric-radius project, so this result moved
to THM-825 before publication.

The proved core lemma will say that a subset of `{0,...,5}` is determined by
its cardinality, position sum, and squared-position sum.  Applied separately
to the four mirror-pair states, this makes `S2+M1+M2` a literal layer codec
whenever the layer has at most six positions.  The staircase layer bounds
then give exactness through `n=15`.  The first possible failure is the
length-seven collision

```text
{0,4,5} and {1,2,6}: cardinality 3, sum 9, square sum 41.
```

Still missing at reservation time: the written disjoint-difference proof,
the exhaustive verifier/output, the exact translation to nonfixed and
apex-normalized fixed B2 layers, and the preservation boundary.
