---
id: HYP-2690
title: Half-tiling tile-address dynamic programs for complement-even tournament invariants
status: OPEN; THM-553 proves the local address layer, first grid-symmetric DP verified n<=8
source: codex-2026-06-20-S57
depends_on:
  - THM-553
  - THM-549
  - THM-550
  - THM-513
related:
  - HYP-2685
  - HYP-2686
  - HYP-2687
  - HYP-2688
  - HYP-2689
  - THM-442
---

# HYP-2690 - Half-Tiling Tile-Address Dynamic Programs

Namespace note: concurrent mac-mini work owns HYP-2689/T925 for the ternary
Eisenstein inclusion-exclusion program.  This codex DP thread was renumbered
to HYP-2690/T926.

## Claim Being Tested

The two-clock tile address

```text
(beta,tau) = (a, a+b-1)
```

should support faster incremental computations for complement-even tournament
invariants as `n` grows.  THM-553 proves the local address layer:

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

## First DP Scaffold

`04-computation/half_tiling_address_dp_codex_s57.py` enumerates the
grid-symmetric/self-converse half cube directly by crossing-layer addresses.
For each half assignment, it reflects bits to the discarded side and computes
`c3`, score multiset, and Hamiltonian-path count on the full tournament.

Exact results through `n=8`:

```text
n  half_coords  assignments  Hmax
3            1            2     3
4            2            4     5
5            4           16    15
6            6           64    45
7            9          512   189
8           12         4096   661
```

The `n=8` half-cube leader has `c3=20` and score multiset
`(3,3,3,3,4,4,4,4)`.  This independently exercises the local-address DP
coordinates; KPS has separately pushed the global max-H half-cube verification
further to `n=9`.

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

OPEN.  THM-553 proves the exact local address calculus and the first
grid-symmetric DP now computes `c3` and `H` through `n=8`.  The next test is the
unordered-pair-state quotient for all complement-even tilings, plus a version
that imports the KPS half-tiling codec and n=9 maximizer workflow.

## kps-S20 realization (THM-554): the SCORE-determined half of this program is now a closed engine

The score generating function `Z_n = (prod_{v>=2} x_v) * prod_{tiles (a,b)} (x_a+x_b)` IS the
address DP this hypothesis asked for, for every **score-determined** invariant.  The beta-step
`Z_{n+1}=Z_n * x_{n+1} * prod_{b=1}^{n-1}(x_{n+1}+x_b)` is the incremental crossing-layer update;
the tau-clock complement reflection is the address quotient (2x fold for complement-invariant
observables).  It computes the EXACT c3-distribution to **n=10** (68.7e9 tilings, ~95s) without
enumerating tilings, and gives the PROVED closed form `E[c3]=(C(n,3)+(n-2))/4`.  Engine:
`04-computation/tile_address_score_gf_engine_kps.py`.  This settles the program's c3/score layer;
the OPEN remainder is exactly HYP-2690 step 2 — extend the state beyond scores (alpha_2 / disjoint
cycle pairs) to reach OCF/H, which `Z_n` alone cannot (the score->H / cut-space->cycle-space wall,
THM-554 Scope).  A kps-S20 application workflow is probing how far a richer-than-score address
state reaches toward alpha_2.
