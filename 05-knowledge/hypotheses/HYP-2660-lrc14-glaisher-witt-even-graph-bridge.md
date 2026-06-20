---
id: HYP-2660
title: LRC14 Glaisher-Witt-even graph bridge - parity and ghost closure as proof quotient
status: OPEN; executable invariant bridge
source: codex-2026-06-19-S37
depends_on:
  - HYP-2651
  - HYP-2656
  - HYP-2658
  - HYP-2659
  - HYP-2648
  - HYP-2652
related:
  - HYP-2655
  - HYP-2657
  - HYP-2653
  - HYP-2654
  - THM-541
  - THM-542
  - THM-543
  - THM-544
  - OPEN-Q-108
---

# HYP-2660 - LRC14 Glaisher-Witt-Even Graph Bridge

## Claim

The useful common structure behind the user's three equinumerosities is not just
"same count"; it is closure after choosing a gauge.

1. Euler/Glaisher: odd-part multiplicities close into dyadic distinct parts by
   binary expansion.
2. Witt/Laurent: Euler products close under the logarithmic ghost map
   `D = x d/dx`, turning product exponents into divisor sums.
3. Tournaments/even graphs: arbitrary tournament bits close to an even graph by
   adding one parity-root whose incident edges are forced.
4. LRC14: the proof quotient should close speed data into a Glaisher tower word,
   then keep CRT residue cells and endpoint owners before scalarizing.

For the one-hole AP-window layer, this exactly reframes KPS HYP-2656: the odd
speeds are the rigid skeleton, and the even speeds are dyadic tower refinements.
The drop-6 minimizer deletes the level-1 bit in the odd-3 tower while both `3`
and `12` survive.  The caveat is important: dyadic tower parity is not a total
order by itself.  Drop `8=2^3` is an even high-level outlier, so the live proof
object must be

```text
Glaisher tower word + CRT cell + endpoint-owner ledger,
```

not a raw even/odd scalar.

## Evidence

Script:

```text
04-computation/lrc14_glaisher_witt_even_graph_bridge_codex_s37.py
```

Stored output:

```text
05-knowledge/results/lrc14_glaisher_witt_even_graph_bridge_codex_s37.out
```

The script verifies `p_distinct(n)=p_odd(n)` for `n<=40`; selected values are

```text
n= 7: 5
n=14: 22
n=21: 76
n=28: 222
n=35: 585
n=40: 1113
```

The actual bijection is the Glaisher map: if an odd part `m` occurs with
multiplicity `a = sum epsilon_j 2^j`, replace those copies by the distinct parts
`epsilon_j 2^j m`.  This is exactly the speed-tower address used by the LRC14
single-deletion layer.

The logarithmic/Witt check finds no mismatch for `m<=40` in

```text
m [x^m] log prod_{n>=1}(1+x^n)
  = sum_{d|m} d (-1)^(m/d+1)
  = sigma_odd(m).
```

This is the Laurent-polynomial shadow of the Witt algebra: before invoking the
full derivation algebra of Laurent polynomials, the useful operator here is the
ghost/log-derivative mode `D=x d/dx`, which turns products into additive divisor
coordinates.

The LRC14 single-deletion Glaisher ledger reproduces the known fixed-observer
safe measures and adds tower addresses:

```text
drop  6: 7/858       odd=3 level=1 tower_after=(0,2)
drop 12: 426/35035   odd=3 level=2 tower_after=(0,1)
drop 10: 1520/63063  odd=5 level=1 tower_after=(0,)
drop  4: 97/4004     odd=1 level=2 tower_after=(0,1,3)
drop  2: 11/364      odd=1 level=1 tower_after=(0,2,3)
drop  8: 950/21021   odd=1 level=3 tower_after=(0,1,2)
```

The odd skeleton alone has safe measure

```text
75454/315315 = 0.2392972107.
```

Thus the tower coordinates explain why the dangerous rows delete dyadic
refinements while leaving the odd skeleton, but the high-level drop-8 outlier
shows why HYP-2658's endpoint-owner ledger, HYP-2659's odd-shell carry carrier,
and HYP-2651's exact components remain necessary.

The tournament/even-graph check verifies for `n<=6` that tournaments on `n`
vertices are in bijection with even graphs on `n+1` vertices:

```text
2^binom(n,2) = 2^(binom(n+1,2)-(n+1)+1).
```

Explicit bijection: mark each anti-order tournament edge as an internal graph
edge, then attach a root edge to exactly the vertices with odd internal degree.
The root degree is automatically even.  Inverse: delete the root and read the
internal parity bits.  The parallel with LRC14 is that raw pair bits are not the
right invariant until the parity closure has been imposed.

## Proof Route

The OPEN-Q-108 proof target suggested by this packet is:

1. Use the Glaisher tower word as a finite address for AP-window and AP-tail
   rows, especially rows with dyadic refinements of the deleted collar speed.
2. Use CRT residue cells to separate low-surviving-cell clean products from
   high-L outliers such as drop 8.
3. Use endpoint-owner ledgers from HYP-2651/HYP-2658 and the HYP-2659 odd-shell
   carry profile to decide the actual lower bound before passing to scalar
   excess.
4. Use the Witt ghost/Laurent derivative viewpoint to compare product-like
   wall functions by divisor coordinates rather than by raw product expansions.
5. Use the tournament/even-graph bijection as the model for proof bookkeeping:
   choose free anti-order bits, then add the forced parity-root closure.  In LRC
   terms, choose local wall/tower defects, then add the forced CRT/endpoint
   closure.

This does not prove LRC(14).  It gives a concrete invariant hierarchy for the
remaining far/state-template layer:

```text
glaisher_tower_word
> crt_residue_cell
> endpoint_owner_ledger
> even_graph_parity_closure
> witt_ghost_divisor_sum
> laurent_wall_polynomial
> raw_sumset_excess
> raw_speed
```

## Tournament Analysis

Vertices are proof quotients:

```text
glaisher_tower_word
crt_residue_cell
endpoint_owner_ledger
even_graph_parity_closure
witt_ghost_divisor_sum
laurent_wall_polynomial
raw_sumset_excess
raw_speed
```

Pairwise observable: which quotient preserves the LRC14 lower-bound predicate
while distinguishing the one-hole AP-window rows.

Switch/gauge: apply parity/ghost closure before scalarization.  Raw speed and raw
sumset excess are terminal quotients, not primary vertices.

Tie Hamiltonian path:

```text
glaisher_tower_word
> crt_residue_cell
> endpoint_owner_ledger
> even_graph_parity_closure
> witt_ghost_divisor_sum
> laurent_wall_polynomial
> raw_sumset_excess
> raw_speed
```

Fingerprint: transitive proof-priority tournament, score histogram
`{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}`, one Hamiltonian path, no directed
3-cycles, singleton SCCs.

## Assumption Challenge

Alternate vertex sets considered: runners, gaps, speed towers, fixed circle
sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier/Witt ghost modes, matroid circuits, even-graph parity bits, and proof
obligations.

Preserved predicate: the fixed-observer safe set `G_C` together with addressed
component/tower ownership.

Destroyed information: exact endpoint positions once one passes to raw partition
counts, raw tournament counts, or scalar excess.

Challenged assumption: tournament vertices should be runners or arcs.  Here the
live quotient is closer to a parity-closed tower word: Glaisher closure for
dyadic speeds, even-graph closure for tournaments, and Witt ghost closure for
Euler products.
