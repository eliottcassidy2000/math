---
id: HYP-9842
title: "After THM-5302: (a) no tournament has three arcs whose eight joint orientations give eight consecutive odd Hamiltonian-path counts (true through order 7, and for three contacts of one head through base order 9; nearest cube 9, 11, 13, 15, 17, 19, 23, 27); (b) line banks (seeds with collinear signatures, exteriors on p alpha - q beta = rho) of unbounded length exist, or banks longer than 3 on a two-parameter family of exteriors; (c) some family of pointed tournaments closed under the tape letters realises a count-doubling generator of Price's Fibonacci-box tree"
status: >
  OPEN (CONJECTURE (a) / QUESTIONS (b), (c)), mac-mini-2026-11-15-monoid. (a) FINITE-EXACT: 606,480 arc cubes of
  order 7 and all smaller orders; there is no congruence obstruction (1,012 cubes of order 7 permute the odd
  residues mod 16), so a proof must be metric. Two independent digits of weights 2 and 4 do exist (counts 9, 11,
  13, 15). (b) NUMERICAL: longest line bank found is 5 (seeds of order <= 9, h <= 45, directions with p, q <= 4).
  (c) by Price's theorem the generator cannot exist on a pair whose second entry is the half-count for all
  three moves at once; which of the two doubling moves exists where is the port-spectrum question of THM-5272.
  Chan-Pak remains OPEN.
session: mac-mini-2026-11-15-monoid
date: 2026-11-15
source: 05-knowledge/results/chanpak_seed_monoid_line_banks_20261115/README.md
scripts:
  - 05-knowledge/results/chanpak_seed_monoid_line_banks_20261115/scripts/arc_bank.c, line_bank.py, cross_bank.c
related:
  - THM-5302, THM-5282, THM-5272, HYP-9822, HYP-9832
---

# HYP-9842 - the third digit, long line banks, a count-doubling generator

Source: [the note](../results/chanpak_seed_monoid_line_banks_20261115/README.md), sections 3-5.

**(a)** For every tournament and every three of its arcs, the eight counts obtained by orienting the three arcs
in all ways are not eight consecutive odd numbers. Falsifier: one such cube (order 8 is the first open case).

**(b)** Is the length of line banks unbounded? Is there a bank of length at least 4 that holds for all pairs of
tape words?

**(c)** Is there a family of pointed tournaments, closed under adding tape vertices, on which one of the maps
`(b, c) -> (c - b, 2c)`, `(c + b, 2c)` is realised by a tournament operation?

## Update (2026-11-16, mac-mini-2026-11-16-dials, wave 41 of chanpak_waves_20261113): order 8 is clean

FINITE-EXACT. `arc_bank 8 3` was run on all 6,880 classes of order 8: 22,538,880 arc cubes, 14,772,792 of them with
eight distinct counts.
- **No complete 3-arc bank.** The least half-spread is 9, attained by 256 cubes; a bank needs 7.
- **The nearest cubes repeat at larger counts.** Besides the order-7 cube 9, 11, 13, 15, 17, 19, 23, 27 they include
  123, 127, 131, 133, 135, 137, 139, 141 and 127, 129, 133, 135, 137, 139, 141, 145.
- **Still no congruence obstruction.** 55,284 cubes permute the odd residues mod 16 (949,144 are balanced mod 8).

So (a) holds through order 8. The first open case is now order 9. A structural reason is proposed in HYP-9852:
additive cubes with three moving arcs avoid small gaps.
