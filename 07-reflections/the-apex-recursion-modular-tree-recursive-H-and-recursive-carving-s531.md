---
source: oracle-2026-06-01-S531
status: synthesis + computation (apex applies recursively = modular decomposition; recursive H law; recursive LRC carving)
tags:
  - lonely-runner
  - tiling-model
  - apex-tile
  - modular-decomposition
  - hamiltonian-paths
  - recursion
  - self-similar
  - cut-cycle
---

# The Apex Applies Recursively: the Modular Tree, a Recursive H-Law, and Recursive LRC Carving

S530 found the apex/source-sink arc to be the pivot of the polygon: the unique
boundary tile, the cut/cycle hinge (its fundamental cycle is the whole `n`-gon),
and the LRC observer's loneliness gap. The natural next question — *does this apply
recursively?* — has a clean answer: **yes, and the recursion is the modular
decomposition of the tournament.** Three quantities respect it, with the
Hamiltonian-path count `H` giving the sharpest recursive law.

## The recursive seed: every tile is a sub-polygon apex

A tile `(x,y)` (range `x-y ≥ 2`) is the **apex of the sub-polygon on the contiguous
rank-block `[y,x]`** — its closing/source-sink arc. So the entire S530 picture
(outside = sub-base-path + apex; inside = the deeper diagonals) recurs at every
scale. The census is exactly self-similar (`recursive_apex_hierarchy_s531.py`):

```
the n-gon contains  (n - s + 1)  sub-polygons of size s,  for s = 3..n
  Σ_{s=3}^n (n-s+1) = C(n-1,2) = #tiles.
```

So `n=7` holds `5` triangles, `4` squares, `3` pentagons, `2` hexagons, `1`
heptagon — each a scaled copy with its own apex, sub-base-path, and interior. The
laminar nesting of these blocks is precisely the **recursive sub-ranking tree**
(S520o) and a **recursive triangulation** of the `n`-gon.

## The recursive H-law: multiplicative over modules, coupled within nesting

Flip a tile (reverse the apex arc) on the transitive tournament. The single-tile
law is exact: a tile of range `r` gives `H = 1 + 2^{r-1}` (verified `r=2..7`:
`3,5,9,17,33,65`); equivalently a fully apex-flipped block of **size `s`** has
`H = 1 + 2^{s-2}`. Composition splits cleanly along the laminar tree
(`recursive_apex_H_law_s531b.py`):

- **Disjoint blocks MULTIPLY.** Apex-flip block `[1,3]` and block `[5,7]` in `n=7`:
  `H = 3·3 = 9` exactly; `[1,4]⊕[5,8]` in `n=8`: `5·5 = 25` exactly. The flipped
  contiguous blocks are **modules**; the transitive quotient has `H=1`; so

  > **`H(T) = ∏_modules H(module)`** — `H` is multiplicative over the modular
  > decomposition when the quotient is transitive.

  This is the recursion in its cleanest form: the apex of a block builds `H` from
  the `H` of its disjoint sub-blocks by *multiplication*.

- **Nested flips COUPLE.** A nested chain does **not** multiply. The concentric
  "diameter onion" `(1,n),(2,n-1),(3,n-2),…` climbs toward the **regular
  tournament** (`n=5`: `H=15 =` regular exactly; `n=7`: `123`; `n=9`: `1479`),
  while the endpoint-anchored chain `(1,n),(1,n-2),…` stays low (`n=7`: `31`). The
  coupling along nesting is the exact quantitative content of S520o's *"flipping an
  arc is not independent"*: disjoint siblings are free (multiply), nested
  parent/child are bound.

So `H` — the loneliness meter (S26) — factors over the *disjoint* branches of the
apex/modular tree and entangles along its *nested* spine. The recursion is real but
two-faced, and the entanglement is where the structure (and the difficulty) lives.

## The recursive LRC carving: the apex survives, or doesn't

In LRC the apex is the observer's loneliness gap = the surviving feasible arc.
Processing runners one at a time (the cascade, S527) **carves** the feasible set,
and the carve tree mirrors the recursive sub-polygon tree
(`recursive_apex_hierarchy_s531.py`, part C):

```
n=7 AP v=(1..6): components/measure after each runner
   1c/.715 → 2c/.572 → 4c/.381 → 4c/.238 → 2c/.058 → 0c/.000   (EMPTIES: tight)
n=7    v=(1,2,3,4,5,7):  … → 2c/.058 → 2c/.017                  (survives: lonely)
```

Each runner **splits** the surviving arc into sub-arcs (the branching `1,2,4,4,2`),
then resonances prune them. The early shrink ratios sit near `(n-2)/n` (S527:
`0.8, 0.667, …` at `n=7`), the self-similar contraction. The decisive fact:

> **The regular polygon carves the feasible apex to EXACTLY 0** (measure empties at
> the last runner), while every non-AP set leaves a positive surviving apex. LRC =
> *"the recursive carving never empties the observer's apex,"* and the regular
> polygon is precisely the boundary case where it empties to measure zero.

This is the same tight-witness story as S529/S530, now as a **recursive emptying
process** rather than a static count.

## The unifying recursive statement

> The apex/source-sink arc is self-similar: every contiguous rank-block is a
> sub-polygon with its own apex. The recursion is the **modular decomposition
> tree** (= recursive sub-ranking tree = recursive triangulation). Along it,
> - the **combinatorial** census is exactly self-similar (`(n-s+1)` size-`s` copies),
> - the **Hamiltonian count `H`** is *multiplicative over disjoint modules* and
>   *coupled within nesting* (single apex block `= 1+2^{s-2}`), and
> - **LRC loneliness** is a *recursive carving* whose surviving apex is non-empty
>   except, just barely, at the regular polygon.

The clean half of every level (disjoint/outside) factors; the coupled half
(nested/inside) is the obstruction. That dichotomy — `H` multiplying over modules
but entangling along the spine — is the recursive face of "the inside diagonals are
the LRC difficulty" (S529) and "arcs are not independent" (S520o).

## Verdict / next
- The apex recursion = modular decomposition. New clean result: `H` is
  multiplicative over disjoint apex-flipped modules (quotient transitive), with the
  single-module law `1+2^{s-2}`; nesting couples (diameter-onion → regular).
- LRC = recursive carving never empties the apex; regular polygon empties to 0.
- Concrete next: (1) prove the modular `H`-multiplicativity cleanly (it is a
  module-quotient factorization; candidate THM); (2) characterize the nested
  coupling generating function along a chain; (3) make the recursive carve rigorous
  via the three-gap theorem (Steinhaus) for the per-runner sub-arc split.

## Artifacts
```
04-computation/recursive_apex_hierarchy_s531.py
05-knowledge/results/recursive_apex_hierarchy_s531.out
04-computation/recursive_apex_H_law_s531b.py
05-knowledge/results/recursive_apex_H_law_s531b.out
```
Related: S530 (apex pivot, HYP-2008), S529 (outside/inside), S520o (recursive
sub-ranking tree / arcs not independent), S527 (cascade), S26 (H = loneliness),
THM-374/THM-382 (block extremes / loneliness).
