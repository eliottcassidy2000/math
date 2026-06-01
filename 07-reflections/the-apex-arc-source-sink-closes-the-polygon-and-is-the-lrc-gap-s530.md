---
source: oracle-2026-06-01-S530
status: refinement + computation (apex = source-sink arc closes the polygon; LRC loneliness gap)
tags:
  - lonely-runner
  - regular-polygon
  - tiling-model
  - apex-tile
  - source-sink
  - cut-cycle
  - staircase
---

# The Apex Arc (Source↔Sink) Closes the Polygon — and It Is the LRC Loneliness Gap

A correction to S529, which said "the polygon's outside (sides) = the base path."
That is incomplete by **exactly one arc**, and that one arc is the whole point.

## The fix: outside = base path + the source-sink arc

The polygon boundary is a Hamiltonian **cycle** (`n` arcs). The tiling model's base
path `n → n-1 → … → 1` is only `n-1` arcs — an open path. The arc that **closes**
the path into the cycle is the one between the **source** (vertex `n`, top of the
ranking) and the **sink** (vertex `1`, bottom):

> **OUTSIDE (the `n` polygon sides) = base path (`n-1` arcs) + the source-sink arc (`1`).**

Verified `n=3..8` (`lrc_apex_arc_source_sink_s530.py`): the base path arcs are all
cyclic-skip-1, and adding the single source-sink arc gives exactly `n` boundary
sides.

## Why that one arc occupies such an important place

The source-sink arc is `(n,1)`, with ranking-range `x-y = n-1`. In the tiling
model this is **not a base-path arc — it is a TILE**, and the **maximal-range
tile**: the **apex** of the staircase triangle. Three facts make it the hinge of
the whole structure (all verified):

1. **It is the unique tile on the boundary.** Among all `C(n-1,2)` tiles, the apex
   `(n,1)` is the only one with cyclic skip `1` (a polygon side). Every *other*
   tile has skip `≥ 2` — an **interior diagonal**. Hence:
   ```
   OUTSIDE sides   = base path (n-1) + apex (1)      = n
   INSIDE diagonals = all other tiles = C(n-1,2) - 1 = n(n-3)/2   (verified n=3..8)
   ```
2. **Its fundamental cycle is the entire polygon.** Adding the non-tree edge
   `(x,y)` to the base path creates the contiguous cycle of length `x-y+1`. For the
   apex `(n,1)` that length is `n` — the **whole `n`-gon**. Every other tile spans a
   shorter contiguous sub-polygon. So the apex is the cut/cycle **hinge**: a
   cycle-space generator (a tile) that is geometrically a boundary (tree-like) edge.
   It straddles the ranking world (where it is the deepest, longest tile) and the
   cyclic world (where it is just the closing side).
3. **It explains S529's "inside born at `n=4`."** For `n=3` the *only* tile is the
   apex `(3,1)` — and it sits on the boundary. So the triangle has **no interior
   diagonal**, no inside debt: S526's single order-2 Legendre term closes it
   completely. At `n=4` the tiles are the apex `(4,1)` **plus** two genuine interior
   diagonals `(3,1),(4,2)` — the inside is born. The apex is the "last outside
   tile"; the interior begins the moment a second tile appears.

This is exactly the corner the recursive staircase decomposition (CLAUDE.md
S20cv-cw) splits off on its own: *overlap + bottom wiring + top wiring + **apex
(arc between the two extreme vertices)***. The apex tile is also the one the repo
already flags as generating the most SC-NS ("expressive") edges — the most global
flip, because its cycle is the whole polygon.

## The LRC payoff: the apex arc IS the observer's loneliness gap

Root the cyclic order of the `n` points (observer at `0` + the `n-1` runners) at
the observer. The **source** is the leading point, the **sink** the trailing one,
and the **source-sink arc is the largest wrap-around gap** — the arc you cut to
linearize the cycle into the ranking path. **The observer lives inside this apex
arc.** Writing `a(t)` = distance to nearest runner ahead, `b(t)` = nearest behind:

> **Observer is lonely ⟺ the apex arc bracketing it clears `1/n` on both sides**
> (`a ≥ 1/n` and `b ≥ 1/n`; apex width `a+b ≥ 2/n`).

Verified: the `n` circle-arcs (the `n-1` inter-runner gaps = the runner *base
path*, plus the `1` observer-bracketing *apex*) sum to `1`, and loneliness is
decided **entirely by the apex**. Examples (`s530.out`): `n=4, v=(1,3,5)` is lonely
at `t=0.45` with apex `a=0.25, b=0.55`; `n=5, v=(1,2,3,5)` lonely at `t=0.24`. The
inter-runner base path must compress to leave the apex open.

At the **regular polygon** every gap equals `1/n`, so the apex clears exactly `1/n`
on each side — the tight, measure-zero boundary case (`v=(1,2)`, `(1,2,3)`,
`(1,2,3,4)` all hit `min(a,b)=1/n` exactly and never strictly exceed it). The apex
sitting precisely at the critical clearance is, once more, why the regular polygon
is the extremal LRC configuration.

## The unified statement

> **LRC = the source-sink (apex) arc can always be opened.** The base path is the
> inter-runner ranking; the apex is the observer's gap; together they tile the
> circle (sum `1`). Loneliness is the apex arc reaching clearance `1/n` on both
> sides — and the regular polygon shows that, in the worst case, it does so only
> just barely. The apex is simultaneously the polygon's closing side (outside), the
> maximal staircase tile (deepest inside the ranking), and the observer's
> loneliness gap. That triple identity is why "that one arc" is the important one.

## Verdict / next
- Refined S529: outside = base path **+ apex (source-sink) arc**; the apex is the
  unique boundary tile, the cut/cycle hinge, and (in LRC) the observer's loneliness
  gap. The `n=3` "no inside" fact is exactly "the only tile is the apex."
- Concrete next: track the **dynamic** apex (which runner-pair brackets the
  observer changes at wall-crossings) as an interval-exchange on the source-sink
  arc; the bounded-ansatz program (S514/S519) is precisely "show the apex cannot
  stay below `2/n`-with-clearance for a full lap."

## Artifacts
```
04-computation/lrc_apex_arc_source_sink_s530.py
05-knowledge/results/lrc_apex_arc_source_sink_s530.out
```
Related: S529 (outside/inside, resonance debt), S520o (recursive sub-ranking tree),
S511 (observer-source = lonely), S514/S519 (bounded-ansatz / wall-order), THM-386.
