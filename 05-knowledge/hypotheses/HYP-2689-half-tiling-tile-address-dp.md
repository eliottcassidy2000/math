---
id: HYP-2689
title: Half-tiling tile-address dynamic programs for complement-even tournament invariants
status: OPEN; THM-551 proves the local address layer, invariant-level DP remains to be built
source: codex-2026-06-20-S57
depends_on:
  - THM-551
  - THM-549
  - THM-550
  - THM-513
related:
  - HYP-2685
  - HYP-2686
  - HYP-2687
  - HYP-2688
  - THM-442
---

# HYP-2689 - Half-Tiling Tile-Address Dynamic Programs

## Claim Being Tested

The two-clock tile address

```text
(beta,tau) = (a, a+b-1)
```

should support faster incremental computations for complement-even tournament
invariants as `n` grows.  THM-551 proves the local address layer:

- `beta` is the full-staircase birth strip;
- `tau` is the half-tiling mirror crossing;
- `(beta,tau)` recovers the tile and its one-flip interval-root facts.

The open claim is that many cycle-space computations can be routed through the
new `tau=n` crossing layer before expanding to the full strip.  On the
grid-symmetric/self-converse subcube this is a literal one-bit-per-orbit update.
On the full complement quotient, nonfixed mirror pairs still carry unordered
pair states, so the saving is an address/state compression rather than a naive
one-bit replacement.

## Exact Evidence

When passing from `n-1` to `n`:

```text
full new cells = n-2,
half-address crossing coordinates = floor((n-1)/2).
```

The half update is not merely "half as many cells"; it is the fixed-line
crossing layer.  Every older discarded-side tile already has a canonical
representative under

```text
(beta,tau) -> (n+beta-tau, 2n-tau).
```

For one-flip packets, the address alone gives

```text
gap = 2beta-tau-2,
c3 = gap,
H_1flip = 1+2^gap,
score defect = e_(tau-beta+1)-e_beta.
```

So at least the FKN radius-1 ledger can be updated exactly from address
layers.

## Program

1. Build an address-indexed DP for the grid-symmetric/self-converse subcube,
   where each mirror orbit is genuinely one bit.
2. Extend to complement-even full-quotient computations by replacing nonfixed
   mirror pairs with unordered pair states.
3. Treat fixed-line `tau=n` cells as the new crossing coordinates when `n`
   increments.
4. Store discarded-side contributions through their reflected representative
   address, not by duplicating full cells.
5. Add THM-513 interval-root packet data to each address before scalarizing.
6. Compare against KPS HYP-2688: if H-maximizers really lie in the
   phi-self-converse half cube, this address DP gives the search coordinates.

## Tournament Analysis

Vertices are address/update layers:

```text
birth_strip
crossing_spine
mirror_orbit
gap_root
endpoint_score
cycle_packet
complement_even_dp
```

Pairwise observable: shared preserved predicates, predicate count, declaration
order.  Switch/gauge: larger observable orients the edge.  The S57 script uses
this only as a workflow ranking: establish birth/crossing/mirror address first,
then import root and cycle packets.

## Assumption Challenge

Do not assume the DP vertices are tournament vertices or arcs.  The natural
vertices for this route are address clocks, mirror orbits, crossing layers,
interval roots, and invariant packets.  The quotient preserves complement
pairs and fixed-line events; it destroys the independent label of the discarded
side.  That is acceptable for complement-even invariants and dangerous for
score-odd or orientation-sensitive quantities.

## Status

OPEN.  THM-551 proves the exact local address calculus.  The next test is a
small `n<=8` DP that computes `c3` and `H` on the grid-symmetric half cube using
crossing layers, then separately tests the unordered-pair-state quotient for
all complement-even tilings.
