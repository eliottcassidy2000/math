---
id: THM-5257
title: "Tree tournaments. For a rooted tree P and a linear extension l (a build order), let T(P, l) orient ancestor -> descendant and, for incomparable pairs, later-built -> earlier-built (the new node beats everything but its ancestors). (1) The 3-cycles of T(P, l) are exactly the triples x <_P y, z incomparable to x, l(x) < l(z) < l(y). (2) T(P, l) is transitive iff l is a preorder (depth-first) traversal; hence H(T(P, l)) = 1 iff the build is depth-first, and otherwise H >= 1 + 2 c_3(l). (3) The tournament depends only on (P, l), not on the labels. (4) FINITE-EXACT: for the Collatz inverse tree built breadth-first by distance from 1 (evens before odds, then value) the counts H(T_N), N = 1..20, are 1 (x8), 5, 11, 41, 131, 517, 2107, 10977, 58679, 334877, 2351809, 13778505, 87545175; with odds first 1 (x7), 3, 5, 15, 41, 157, 517, 2793, 12359, 70429, 332969, 2166681, 11604367, 105295407; the depth-first build gives 1 at every N."
status: >
  PROVED (klein-2026-11-16-hex; elementary case analysis on the number of comparable pairs in a triple; not independently
  audited). FINITE-EXACT: (4) by subset DP (scripts/treetour.c), brute-force checked to N = 9; the characterisation (2)
  and the 3-cycle description (1) checked on 406 random (tree, linear extension) pairs of 4-7 nodes (scripts/dfs_check.py).
  Pre-registered P-D, P-E confirmed (packet PREREGISTERED.md). CITED: THM-002 (odd-cycle formula), THM-4876 (path pairs).
session: klein-2026-11-16-hex
date: 2026-11-16
source: 05-knowledge/results/collatz_tree_tournament_sink_hex_20261116/COLLATZ_TOURNAMENT.md
scripts:
  - 05-knowledge/results/collatz_tree_tournament_sink_hex_20261116/scripts/treetour.c (+ receipts/treetour_*.out)
  - 05-knowledge/results/collatz_tree_tournament_sink_hex_20261116/scripts/dfs_check.py (+ receipts/dfs_check.out)
related:
  - 01-canon/theorems/THM-002-ocf.md
  - 01-canon/theorems/THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count.md
  - 05-knowledge/results/path_splice_rows_continuation_20261116/COLLATZ_ROWS.md
  - 05-knowledge/results/session140_collatz_hexagonal_geometry_20261102/README.md
---
# THM-5257 — tree tournaments built one node at a time

**Proof of (1).** Three vertices span 0, 1, 2 or 3 comparable pairs of `P`. With 0 pairs the orientation is by `ℓ`
alone, transitive. With 3 they form a chain, transitive. With 2, the two comparable pairs share a vertex which is an
ancestor of the other two (two ancestors of one vertex are comparable), a source of the triple. With exactly 1, say
`x <_P y` and `z` incomparable to both: `x → y`, and `y → z → x` iff `ℓ(y) > ℓ(z) > ℓ(x)`. ∎

**Proof of (2).** A tournament is transitive iff it has no 3-cycle. If `z` is incomparable to `x` it is incomparable
to every descendant `y` of `x` (an ancestor of `y` is comparable to `x`; a descendant of `y` is a descendant of `x`).
So the absence of the triples of (1) says that no node incomparable to `x` is built between `x` and any descendant of
`x`, for every `x`: the subtree of `x` is a contiguous block of `ℓ` starting at `x`, which is the preorder property.
`H = 1` iff transitive; the bound is THM-002 with one triangle. ∎

**(3)** is immediate from the definition. **(4)** receipts `treetour_even.out`, `treetour_odd.out`, `treetour_dfs.out`.

**What is new.** The characterisation and the count sequence. Known and used: the insertion argument, THM-002.
