# The Collatz tree as a tournament built one node at a time, the sink trick with a larger prefix, the honeycomb rows, three tilings into five solids, and the port-move search

Session klein-2026-11-16-hex (2026-11-16), on the owner's pasted directive (Collatz as a hexagonal tiling of edges;
Collatz globally as an ordering by distance from 1 with an even/odd tie-break and locally as a comparison; the
tournament on infinite nodes, its Collatz-tree subset, built one node at a time, in other orders, in groups,
relabelled; the hexagonal, triangular and square lattices against the Platonic solids and the pattern of splitting 3
into 5; the binary port-move search or something better) and the owner's line "try the sink trick with a larger
prefix and a fresh sink before each digit". Canon: THM-5257 (new).

## Status

- **PROVED** ([COLLATZ_TOURNAMENT.md](COLLATZ_TOURNAMENT.md), THM-5257). The tree tournament `T(P, ℓ)` of a rooted
  tree built in the order `ℓ` (the new node beats everything but its ancestors) has as 3-cycles exactly the
  ancestor–descendant–interloper triples, and is transitive iff the build is depth-first; so `H = 1` iff depth-first.
  The Collatz distance order is breadth-first and gives the counts 1 (×8), 5, 11, 41, 131, 517, 2,107, …, 87,545,175
  (`N ≤ 20`, evens first) and 1 (×7), 3, 5, 15, 41, 157, 517, …, 105,295,407 (odds first); labels are irrelevant,
  inflation by 2-blocks multiplies by non-powers of two. Pre-registered P-D, P-E confirmed.
- **FINITE-EXACT, negative** ([SINK_DIGITS.md](SINK_DIGITS.md)). The sink trick with a larger prefix: over all rooted
  tournaments of order ≤ 6 (two added vertices), 7 (one), ≤ 5 (three), only count-5 prefixes carry three digits with a
  common continuation (P-A, P-C confirmed); the parity split is not the digit's resource (83 % against 80 % of bases
  have a `3H + 2` extension; P-B refuted).
- **FINITE-EXACT, negative** ([PORT_MOVES.md](PORT_MOVES.md)). The binary odometer search generalised from wave 45's
  two-root appends to arbitrary history-blind letters: no pair of letters survives the words of length two at orders
  3–6 (431,736 pairs at order 6).
- **Synthesis** ([HEX_PLATONIC.md](HEX_PLATONIC.md)). The honeycomb Collatz packets of 2026-11-02/04 and the rows of
  doubling chains are one drawing (rows = even arrows, Syracuse = odd jumps; `{7, 21}` one jump). Three tilings into
  five solids = one self-dual fixed object plus a dual pair splitting into two, on the Schläfli boundary; the same
  `fixed + 2·pairs` shape in the order-4 diamond 1, 3, 3, 5, the five-number window around the Collatz image, the rows
  `L_3`, `L_5` and the Busch floors; typed Platonic/lattice links; questions Q1–Q3.
- Pre-registration: [PREREGISTERED.md](PREREGISTERED.md) (P-A, P-C, P-D, P-E confirmed; P-B refuted).

## Files

- `scripts/`: `treetour.c` (tree tournament counts, orders and rules, inflation), `dfs_check.py` (random checks of
  THM-5257), `allsearch2.c` (sink-digit scans with 1–3 added vertices), `odometer.c` (generalised binary odometer).
- `receipts/`: `treetour_even/odd/dfs/numeric/*.out`, `dfs_check.out` / `.optimized.out` (identical), `as2_*.out`,
  `odometer_345.out`, `odometer_6.out`.
- `prediction_hex.txt`.
