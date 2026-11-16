# The charge-carrying port moves at orders 7 and 8 (FINITE-EXACT)

Session klein-2026-11-16-force. Tool `scripts/portmove.c` (all classes of an order, every ready root, every extension by
one or two vertices with arbitrary contacts, every new root; `seed … one` mode for iteration); receipts
`portmove_7.out`, `portmove_8one.out`, `iter_ports.out`, `iter_ports_summary.out`.

A port move at charge `a` is `r ↦ 2r + 1` at the same charge: `H′ = 2H + 1 + a`, `a′ = a`, `b′ = 0` (wave 45's
coordinates; ready means `b = O − E − g = 0`, i.e. `D = −g`).

| order (added vertices) | ready roots | extensions reaching `2H + 1 + a` | ready same-charge (extension, root) pairs | by parent count and charge |
|---|---|---|---|---|
| 7 (≤ 2) | 367 | 1,367 | 59 | `1` (a = 0): 45; `15` (4): 5; `27` (4): 3; `71` (6): 3; `81` (−2): 3 |
| 8 (1) | 3,571 | 522 | 10 | `1` (0): 3; `15` (4), `27` (4), `51` (8), `71` (6), `81` (−2), `115` (16), `255` (0): 1 each |

- **Classification.** Besides the trivial `1 → 3`, the carrying bases below order 9 are the counts 15, 27, 51, 71, 81,
  115, 255 at the charges 4, 4, 8, 6, −2, 16, 0, each sending `H ↦ 2H + 1 + a`: 15 → 35, 27 → 59, 51 → 111, 71 → 149,
  81 → 161, 115 → 247, 255 → 511. The charge-0 move 255 → 511 is hadamard's Mersenne step `F_8 → F_9` (balanced fan,
  one inserted loser), ready at the same root.
- **None iterates.** From every valid child (35, 59, 111, 149, 161, 247, 511, all ready at the parent's charge) no
  one-vertex extension even reaches the next target `2H′ + 1 + a` (`iter_ports_summary.out`); for 511 this is hadamard's
  observation that the plus sign at odd `k` needs an old-arc reversal. Pre-registered P-C (only charges 0 and 4 at order
  7) is refuted on the charges and confirmed on the absence of iteration.
- With the previous packet's negative results (orders ≤ 6; the 83 one-parent ready seeds) the picture to order 8 is:
  port moves exist sporadically at seven parent counts, never twice in a row.
