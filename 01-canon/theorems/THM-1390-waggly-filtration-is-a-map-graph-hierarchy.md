---
id: THM-1390
title: The waggly filtration is a map-graph hierarchy — the d=1 (wiggly) metagraph is the "edge-contact" dual and the deeper layers are "point contacts", and adding them drives the metagraph to the COMPLETE graph, the clique explosion that separates map graphs from planar duals. This defines a new invariant, the SATURATION DEPTH d_sat(n) = the S_n-quotient's tile-flip diameter, computed exactly: 2, 3, 4, 7 for n = 4,5,6,7. The clean pattern d_sat = n−2 fits n ≤ 6 and is FALSE at n = 7.
status: VERIFIED-EXACT (full enumeration of the tiling hypercube and iso classes for n = 4,5,6,7 — 8/64/1024/32768 tilings, 4/12/56/456 classes, matching A000568; layered edge counts and densities exact; the n=7 diameter computed by exhaustive Hamming-ball intersection, unreachable-pair counts 21165/2687/250/20/0 at d ≤ 3/4/5/6/7). The map-graph correspondence is a structural ANALOGY, not a proved equivalence — see Honest scope.
source: mac-mini-2026-07-20-S126 (owner: extend ideas abstractly similar to map graphs — vertices = faces, edges = faces meeting at a vertex OR edge)
related:
  - the tiling model and waggly-layer structure (CLAUDE.md)
  - METAGRAPH-ATLAS
---

# THM-1390 — the waggly filtration as a map-graph hierarchy

> ### CORRECTION ACCEPTED (2026-07-20, kind-pasteur-S128c108, THM-1400 §I) — `d_sat` IS NOT NEW
>
> **The saturation depth `d_sat` is the metagraph DIAMETER, and the diameter was already
> known.** `G^(≤d)` is complete exactly when every pair of classes is within `d` flips, so
> `d_sat = diam(G_n)` — trivially, once stated. And `opus-2026-03-24-S306` had already
> identified `diam(G_n) = max_T min-FAS(T) = A003141(n)` (growth `~n²/4`), with
> `OPEN-QUESTIONS.md` listing it **RESOLVED** and `README.md` carrying it as the Waggly
> Completeness Theorem. I cited neither. The correction is **right and I accept it in full**:
>
> - `d_sat` is a **rediscovery**, not a new invariant. Strike "new invariant" wherever it
>   appears above and below.
> - My handoff *"compute n=8 before conjecturing"* needs no computation: `A003141(8) = 8`.
> - My *"no linear formula"* is a **known quadratic**, `n(n+1)/4 − Θ(n^{3/2})`.
> - My refutation of `d_sat = n−2` at `n=7` restates opus-S306's *"diam = n−2 is WRONG for
>   large n."* Independent, but not first.
> - My `2,3,4,7` is the **unmerged** diameter; canon's merged values are `1,3,4` at `n=4,5,6`.
>   These agree except at `n=4`, as they must — complementation can only shorten distances.
>
> **What survives:** the map-graph *framing* (the table below, point-contact vs edge-contact,
> the clique-explosion reading of saturation). THM-1400 explicitly grants this. The lesson for
> me is the standing one — **grep the reflections and OPEN-QUESTIONS for the invariant before
> naming it**; the diameter was sitting in canon marked RESOLVED. Logged to MISTAKES.
>
> The point-adjacency thread continued in **THM-1405**, which asks what the *star* (clique-at-
> a-point) move quotients by, and finds a genuine invariant there: the cycle space.

**One line.** A map graph generalizes the planar dual by letting faces meet at a *point*,
not only along an edge, and the payoff is unbounded cliques. The repo's waggly filtration
does the same thing to the metagraph — and the clique explosion is total: the metagraph
saturates to `K_V`.

## The correspondence

| map graphs | the metagraph |
|---|---|
| faces of a planar map | tournament iso classes |
| faces sharing an **edge** (the planar dual) | `d = 1`, the **wiggly** layer (flip one tile) |
| faces meeting at a **point** | `d ≥ 2` waggly layers (flip `d` tiles) |
| `k`-map graph (≤ `k` faces per point) | the truncation `G^{(≤k)}` |
| planar duals have bounded cliques; **map graphs have unbounded cliques** | `G^{(1)}` is sparse; `G^{(≤d)}` becomes **complete** |

## The computation (exact, full enumeration)

Tiling model: base path `n → n−1 → … → 1`, tiles = arcs `(x,y)` with `x−y ≥ 2`,
`m = C(n−1,2)`, all `2^m` tilings enumerated and canonicalised into iso classes
(class counts 4, 12, 56, 456 for `n = 4,5,6,7` — matching A000568).

**Layer-1 density collapses as `n` grows** (the "planar dual" is getting sparse):

| `n` | `m` | classes | `|E(d=1)|` | density of `G^{(1)}` |
|---|---|---|---|---|
| 4 | 3 | 4 | 5 | 0.833 |
| 5 | 6 | 12 | 30 | 0.455 |
| 6 | 10 | 56 | 290 | 0.188 |

**Adding point-contact layers saturates to the complete graph:**

| `n` | cumulative density at `d ≤ 1,2,3,4` |
|---|---|
| 4 | 0.833 → **1.000** |
| 5 | 0.455 → 0.909 → **1.000** |
| 6 | 0.188 → 0.696 → 0.964 → **1.000** |

## The new invariant, and the mirage

> **Saturation depth** `d_sat(n)` := the least `d` with `G^{(≤d)}` complete
> = the diameter of the `S_n`-quotient of `Q_m` in the tile-flip metric
> (the fewest tile flips relating *any* two iso classes, after relabelling).

```
n        :  4   5   6   7
d_sat(n) :  2   3   4   7
m=C(n-1,2):  3   6  10  15
```

`n = 4,5,6` give exactly `n−2`, which is a clean and tempting formula. **It is false.**
At `n = 7` the value is **7**, not 5: exhaustive Hamming-ball intersection leaves
21165 / 2687 / 250 / 20 unreachable class pairs at `d ≤ 3/4/5/6`, reaching 0 only at
`d = 7`. So `2,3,4,7` has no linear formula, and this is another instance of the
small-`n` trap that has recurred across this project's extremal work — a three-point
pattern that dies at the fourth point.

What survives is weaker but real: `d_sat` stays well below `m`, so the quotient does
compress the hypercube's diameter (`15 → 7` at `n = 7`), just not linearly.

## Honest scope

- **Verified:** every number above, by full enumeration; the class counts independently
  match A000568, which validates the canonicalisation.
- **Analogy, not theorem:** "the waggly filtration *is* a map-graph hierarchy" is a
  structural correspondence, not a proof that `G^{(≤k)}` is a `k`-map graph in the
  Chen–Grigni–Papadimitriou sense (which would require exhibiting a planar bipartite
  half-square representation). That is the natural next question and is **not** settled
  here; note the metagraph's saturation to `K_V` is consistent with it, since `K_V` *is*
  a map graph for every `V`.
- **A caught error:** the first `n = 7` run reported diameter 99 because I had truncated
  its loop (`if mask > 4000: break`) for tractability — an artefact of my own code, not a
  property of the object. The value 7 comes from the corrected exhaustive computation.

*Artifacts:* `04-computation/waggly_map_graph_hierarchy_macmini_S126.py` (+out).
