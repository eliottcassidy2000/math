# Residues by ordinal sums, the hidden generators 35 and 49, the converse census of 33, the mod-3 tower and the pruned Collatz tree, and P107 on the ready seeds

Session klein-2026-11-16-monoid (2026-11-16), on the owner's directive (three tilings into five solids with the
corrected correspondence; 33 at order six as one self-converse class plus two converse pairs; odd numbers multiplied
or stepped by ±2 to reach the primes; prove every odd residue modulo every `2^a 3^b 5^c` realised and the pattern for
primes to 63 skipping 7; micro/macrocosm and the field with one element; the first 32 odd numbers split 3 | 1 | 28;
`{35, 49}` as precursors of an infinite hidden set; the root-moving port move on the ready seeds) and the pasted
second directive (multiplication modulo odd integers and Collatz as tricks; the mod-3 tower with period 12; the
accelerated chains and the shoelace; the pruned Collatz tree with edge lengths 1, 2, 3, 4). The hadamard packet
`one_apex_all_residues_20261116` was read first. Canon: THM-5258, THM-5259 (new).

## Status

- **PROVED** ([RESIDUES.md](RESIDUES.md), THM-5258 (A)). Every residue modulo every modulus was already proved
  realised (Session 25; hadamard's one-apex packet at every large order) — cited, not repeated. New: the first 32 odd
  numbers minus the holes generate every unit group modulo `2^a 3^b 5^c`, and the counts below 4000 with prime
  factors ≥ 67 generate every unit group with primes ≤ 63, by the exponent-free Frattini criterion; so every odd
  residue there is the count of an **ordinal sum of bounded blocks**, with the 7-lane supplied by the line `L_3`.
- **PROVED** ([MONOID.md](MONOID.md), THM-5258 (B)). The indecomposables of `odd ∖ {7, 21}` are the primes `≠ 7`, the
  `7p` (`p ≥ 5`), and 49, 63, 343: the hidden generators begin 35, 49 and are as many as the primes; `×7` is
  `t ↦ 7t + 3` in the address `t = (H − 1)/2`, the holes 7, 21 are the images of 1, 3, and the first 32 odd numbers
  split `3 | 1 | 28` at the hole `t = 3` with the block ending at `63 = 7·9`. The "microcosm" is hadamard's §5 (`3 + 7 = 10`
  in addresses); the field with one element is the transitive tournament (typed).
- **FINITE-EXACT** ([CONVERSE.md](CONVERSE.md)). 33 first appears at order 6 as 1 self-converse class + 2 converse
  pairs; full multiplicity tables (orders 4–7; 456 classes, 88 self-converse at order 7); the Platonic correspondence
  corrected as the owner states.
- **PROVED** ([TOWER_TREE.md](TOWER_TREE.md), THM-5259). The pruned accelerated Collatz tree is a complete binary tree
  with edge lengths in {1, 2, 3, 4} determined modulo 27, labels recovered from lengths and left/right; the mod-3
  tower's period-12 classes, shoelace halving arrows (+2 rows per +3) and three odd-arrow families (4 rows up, 2
  longer) — the owner's description made exact.
- **FINITE-EXACT** ([PORTMOVE.md](PORTMOVE.md)). P107 complemented: with arbitrary contacts and a moving root, only the
  trivial `1 → 3` and three one-off charge-4 moves `15 → 35` exist at orders ≤ 6 (none iterates); on the ready seeds
  13, 23, 41, 59, 75, 95 reconstructed as one-parent trees: none of the 83 representatives (every count and charge at orders 5–8) admits a ready same-charge move `r ↦ 2r + 1` with one or two added vertices (P-6, P-9 confirmed). The ready class is `D = −g`.
- Pre-registration: [PREREGISTERED.md](PREREGISTERED.md).

## Files

- `scripts/`: `converse_census.py`, `residues.py`, `frattini3.py`, `indecomposables.py`, `portmove.c` (all-orders and
  `seed` modes), `seeds.py`, `seeds2.py`, `onparent.py`, `collatz_tower_tree.py`.
- `receipts/`: `converse_4..7.out`, `residues.out`, `frattini3.out`, `savoid.out`, `indecomposables.out` /
  `.optimized.out`, `portmove_3..6.out`, `seeds.out`, `seeds2.out`, `onparent.out`, `iterate_children.out`,
  `seedruns.out`, `collatz_tower_tree.out` / `.optimized.out`, `tree_mod27.out`.
- `prediction_monoid.txt`.
