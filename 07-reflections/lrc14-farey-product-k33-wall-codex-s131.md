---
source: codex-2026-06-23-S131
status: new LRC14 proof split from the Farey-product complete-bipartite graph
tags: [lrc14, farey, complete-bipartite, K33, product-ledger, tournaments]
---

# LRC14 Farey Product: the K3,3 Wall

The product payload finally has a clean job.

For `M(S)=p/q`, the theorem still belongs to

```text
M(S)-1/14 = (14p-q)/(14q).
```

So the denominator `q` is not negotiable.  But `p*q` is not just a noisy area
coordinate either.  It is literally the edge count of `K_{p,q}`.

That makes the user's Farey observation click into the LRC14 chain.  Ordinary
Farey order first sees complete-bipartite nonplanarity at `F_4`, through
`3/4 -> K_{3,4}`.  The LRC14 unit-excess chain is stretched by the apex:

```text
p/q = p/(14p-1).
```

The first two escapes are planar:

```text
1/13 -> K_{1,13}
2/27 -> K_{2,27}
```

Then the known near miss appears:

```text
3/41 -> K_{3,41}.
```

That is the first `K_{3,3}` product wall, and it is exactly the `12->36` row
from the S128 Farey-mediant atlas.

This is useful because it separates the unit-excess branch into two different
species.  The `2/27` rows are two-block planar strips; they want the existing
petal/Jacobsthal/two-block rigidity tools.  The `3/41` and higher rows are
three-owner incidence objects; they should be attacked as finite obstruction
packets, potentially feeding the HYP-2908 state-lift endpoint.

The S131 row bank was pleasantly sharp: among `749` AP/GW/petal/single-swap
rows, the only unit-excess `K33-wall` row is `near-miss 12->36`.  The only
unit-excess two-block rows are `10->20` and `13->26`.  Everything tight remains
in the star ledger `K_{1,14}`.

So the split I would carry forward is:

```text
e=0        floor candidate
e=1,p=1    coarse q-threshold parent
e=1,p=2    planar two-block strip
e=1,p>=3   K3,3 wall
```

Nonplanarity is not a proof by itself.  But it is the first honest
three-owner signal in the Farey product ledger, and it gives the next finite
packet a name.
