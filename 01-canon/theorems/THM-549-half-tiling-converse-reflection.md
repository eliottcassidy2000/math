---
id: THM-549
title: The converse-with-relabel acts on fixed-HP tilings as a PURE coordinate reflection rho:(a,b)->(n+1-b,n+1-a) of the staircase (no GF(2) bit-flip); the half-tiling is its fundamental domain
status: PROVED (elementary; verified exhaustively n=3..7 in 04-computation/half_tiling_framework_kps.py)
source: kind-pasteur-2026-06-20
depends_on:
  - definitions.md   # tiling model: base path P_0, tile (a,b), bit convention
related:
  - THM-550   # half-tiling count and grid-sym fraction
  - THM-551   # half-tiling recursions and full<->half folding
  - MISTAKE-033  # complement-tiling vs T^op confusion (this thm pins the exact map)
external: Redei / self-converse tournaments; Burnside on the reversal anti-automorphism.
---

# THM-549 — The converse is a pure reflection of the staircase

**User seed (kind-pasteur-2026-06-20):** "Taking the mirror image of a tiling over the
line y=x is equivalent to reversing all arcs in the tournament (including those in the
fixed path)."  This theorem makes that exact.

## Setup

Fix `n >= 2` and the base path `P_0 = n -> n-1 -> ... -> 1`.  Tiles are pairs `(a,b)`
with `a >= b+2`; bit 0 = forward arc `a->b`, bit 1 = backward arc `b->a`
(definitions.md).  Let `phi(i) = n+1-i` be the order-reversing relabel.  For a tiling
`t`, let `C(t)` be the tiling of the tournament `phi . (T_t)^op` (reverse ALL arcs,
including base-path arcs, then relabel by `phi` so `P_0` is again the base path).

## Statement

`C` is realized as the **pure coordinate involution**
```
        rho(a,b) = (n+1-b, n+1-a),       C(t)(P) = t(rho(P)) for every tile P,
```
with **no GF(2) bit-flip**.  `rho` is an involution of the staircase `delta_{n-2}`.
Its fixed cells are exactly the tiles on the anti-diagonal `a+b = n+1`, and they number
`d = floor((n-1)/2)`.  Equivalently (CLAUDE.md `isGridSym`): `rho` is the reflection
`(x,y) -> (n+1-y, n+1-x)`.

**Corollary.** A tiling is grid-symmetric (`isGridSym`; equivalently `T_t` is
self-converse via the single anti-automorphism `phi`) **iff `t` is constant on every
`rho`-orbit.**

## Proof

**(i) `rho` is an involution on tiles.**  For `(a,b)` with `a >= b+2`, since `b < a` we
have `n+1-b > n+1-a` and `(n+1-b)-(n+1-a) = a-b >= 2`, so `rho(a,b)` is a valid tile.
`rho^2 = id` is immediate.

**(ii) No bit-flip.**  Take tile `(a,b)`, bit 0, i.e. arc `a->b` with `a` the upper
(larger) vertex.  Arc reversal gives `b->a`; relabel by `phi` gives
`phi(b)->phi(a) = (n+1-b)->(n+1-a)`.  In the image tile `{n+1-b, n+1-a}` the upper
vertex is `n+1-b` (the larger), so the arc points upper->lower = **bit 0**.  Thus bit 0
at `(a,b)` maps to bit 0 at `rho(a,b)`; identically bit 1 maps to bit 1.  The two
order-reversals (arc reversal AND label reversal) cancel on the "which endpoint is
upper" bookkeeping, so `C` carries no translation: `C(t) = t . rho`.  Base-path arcs
`k->k-1` map to `phi(k-1)->phi(k) = (n+2-k)->(n+1-k)`, i.e. `(j+1)->j` with
`j = n+1-k`, recovering `P_0`; hence `C(t)` is again a tiling. **QED (ii).**

**(iii) Fixed cells.**  `rho(a,b)=(a,b)` iff `a+b=n+1`.  On `a+b=n+1` with `a>=b+2`:
`a-b = n+1-2b >= 2` iff `b <= (n-1)/2`, so `b in {1,...,floor((n-1)/2)}` and
`d = floor((n-1)/2)`. **QED.**

## Why this matters

- It pins down EXACTLY the map MISTAKE-033 warned about: the staircase reflection is
  the genuine `T^op` (all arcs reversed) **composed with the canonical relabel `phi`** —
  not the "complement tiling" (flip tile bits, keep labels).  Because `phi` undoes the
  base-path reversal, the net action on the *tiles* is a pure permutation.
- The fixed-point set (grid-symmetric tilings) is a coordinate subcube
  `{t : t = t . rho}` of dimension = #orbits = the **half-tiling** dimension (THM-550).
- It is the labeled-tiling-level realization of the fixed-HP model's `Z_2` symmetry
  (07-reflections/full-vs-fixed-hp-tiling-duality.md), i.e. the complement flip of the
  staircase.

## Verification

`04-computation/half_tiling_framework_kps.py`, block [B]: exhaustive over all `2^m`
tilings for `n=3,4,5,6,7` (2,8,64,1024,32768 tilings) confirms `C(t)(P)=t(rho(P))`
for every tiling and every tile (100% match).  Block [C] confirms the fixed tilings are
exactly the `phi`-self-converse tournaments.
