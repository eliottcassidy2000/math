---
source: oracle-2026-06-01-S530
status: result (3 structural identities verified) + synthesis
tags: [tiling-model, apex, source-sink, hamiltonian-cycle, principal-line, block-extremes, LRC, THM-382, polygon]
---

# The source-sink arc: where the polygon boundary, the staircase apex, and the LRC observer-gap are the same arc

**Prompt (user):** the outside of the polygon = the base path of the tiling model
PLUS the one arc between source and sink; that arc closes the boundary cycle and
occupies an important place in the tiling model.

This is exactly right, and the arc turns out to be the same object in five guises.
Fix the base path `0 -> 1 -> ... -> (n-1)` (source `0`, sink `n-1`); tiles are the
non-path arcs `(i,j), j>i+1`; the **source-sink arc `(0,n-1)`** has maximal range
`n-1`.

## 1. Outside = base path + the source-sink arc = the boundary n-cycle

The base path is a Hamiltonian *path* — `n-1` of the polygon's boundary edges.
Adding the source-sink arc closes it into the full **boundary n-cycle** = the
regular polygon's perimeter = the "outside." So the user's correction is precise:
the outside is not the path alone, it is the path **plus this one closing arc**.

## 2. It is the staircase APEX = the principal-line tile (verified)

In the staircase/tiling model the source-sink arc `(0,n-1)` is the **apex tile**:
maximal range, the hypotenuse corner of the triangle `delta_{n-2}`. Its importance
is quantitative — flipping it on the transitive tournament is the **maximal single-
tile H jump**:

```
H(transitive + apex back-arc) = 1 + 2^{n-2}   for n=3..8: 3,5,9,17,33,65.  (verified)
```

That is the transitive class's "big strongly-connected neighbour" on the principal
line (CLAUDE.md: `Delta H = 2^{n-2}` from transitive to first SC neighbour). No
other single tile moves H so far. The source-sink arc is the hypotenuse =
`H = 1 + 2^d` at `d = n-2` (max).

## 3. It is the master switch between the two block-extremes (verified)

Round tournaments (the LRC-realizable necklace, S523) live only at `#SCC in {1,n}`
(S524/S525). The source-sink arc is what selects which extreme. Building round
tournaments from circle gaps (3000 samples each, n=4,5,6, **0 mismatches**):

```
transitive (#SCC = n)        <=>  the largest boundary gap > 1/2  (a SEMICIRCLE; THM-374)
strongly connected (#SCC = 1) <=>  the largest gap < 1/2 (runners wrap; Moon: Ham cycle = boundary)
```

and the **source-sink arc is the chord across that largest gap**. Aligned (over a
big empty gap) -> transitive; flipped into a back-arc (no gap exceeds 1/2) -> the
boundary becomes a directed Hamiltonian cycle -> one strong block. So flipping the
single apex arc is the move between the transitive vertex and the strong-block
extreme — the two ends of the round necklace.

## 4. In LRC it is the observer's loneliness gap

The observer sits at a point on the circle; cut there and the runners get a linear
order whose first and last elements (source, sink) bracket the **gap containing the
observer** — the largest relevant gap. So the source-sink arc is *the observer's
gap*, and by THM-382 loneliness is exactly:

```
observer lonely  <=>  both observer-adjacent gaps >= 1/n
                 <=>  the two endpoints of the source-sink arc are each >= 1/n from the observer.
```

Hence the two lonely regimes are selected by this one arc:
- **apex aligned (big gap) = transitive = observer OUTSIDE the runner cluster** (the
  semicircle case; the runners hide in a half and the observer is lonely in the
  empty half);
- **apex back-arc (all gaps < 1/2) = one strong block = observer INSIDE a wrap-
  around** — the tight `AP / regular-polygon` case (n=14's unique wall-only set,
  S525). The runners surround the observer and loneliness is barely achieved.

## The five guises of one arc

```
source-sink arc (0,n-1)
  = the closing edge of the boundary n-cycle      (polygon outside)
  = the apex tile, max range, hypotenuse          (staircase / tiling model)
  = the maximal single-tile jump  H = 1 + 2^{n-2} (principal line of the metagraph)
  = the switch between #SCC = n and #SCC = 1       (the two block-extremes / round necklace)
  = the observer's loneliness gap                  (LRC, THM-382)
```

That is the "important place" the user pointed to. The whole inside/outside picture
(S529) hinges on this single arc: it is the one piece of the OUTSIDE whose
orientation already decides the gross INSIDE structure (transitive vs one big
cycle), and in LRC it is precisely the arc loneliness is measured on. The tight
extremal case is `apex = back-arc` (the polygon boundary is a full directed cycle =
the regular polygon), which is why the regular polygon / AP keeps reappearing as
the hard witness.

## Anchor
`04-computation/source_sink_apex_arc_s530.py` (+ `.out`): apex flip H=1+2^{n-2}
(n=3..8); transitive <=> largest gap > 1/2 (0 mismatches, n=4,5,6). Builds on the
tiling model (CLAUDE.md), S523/S524 (round, block-extremes), S525 (n=14 wall-only
AP), S529 (outside/inside), THM-382, Moon's theorem.
