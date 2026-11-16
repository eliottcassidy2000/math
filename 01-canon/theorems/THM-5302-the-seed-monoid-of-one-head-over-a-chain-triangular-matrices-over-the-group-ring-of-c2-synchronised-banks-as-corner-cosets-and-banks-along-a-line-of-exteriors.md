---
id: THM-5302
title: "The seed monoid. (A) For a block C with head pattern X put M(C, X) = [[h, b, c], [0, G, a], [0, 0, h]] over R = Z[e]/(e^2 - 1): h = H(C), G, a, b, c the sums of H(C - Q) over the paths Q of C (all paths; ending at a non-contact; starting at a contact; from a contact to a non-contact), e marking an odd number of path vertices. Then M(C1 -> C2, X1 u X2) = M(C1, X1) M(C2, X2) and H(chain + head) = h + 2 [c]_even. (B) Tape letters: non-contact [[1,0,0],[0,1+e,e],[0,0,1]], contact [[1,e,0],[0,1+e,0],[0,0,1]]; under e = 1 they act on (b, c) by (2b, c + b) and (2b + 1, c), and H = 1 + c^+ + #(factors 10) for one-apex words. (C) Line banks: seeds with equal (h, G) and signatures (a, b) on a line of direction (p, -q) have counts differing by 2(c_i - c_0 + t_i rho) for every exterior with p alpha - q beta = rho. (D) Finite-exact: no synchronised bank longer than its block count through block order 9; line banks of length 5; no complete 3-arc bank through order 7"
status: >
  PROVED: (A) by splitting a path of an ordinal sum into its parts (same argument as THM-5282 (A), (D), of which
  it is the closed form); (B) by evaluation; (C) from THM-5282 (A) with the lowest upper vertex a non-contact and
  the highest lower vertex a contact. Checks: (A) 40 random products and 30 random chains against brute force;
  (C) one bank literally at four exteriors of its line.
  FINITE-EXACT: (D) 259,280 seeds on the 786 blocks of order 3..9 with h <= 45: longest synchronised bank by h is
  2 (h = 3), 3 (5), 4 (9), 3 (11), 5 (15), 4 (17), 5 (33), 6 (45), never longer than h; with two arbitrary tapes
  the longest is 3; line banks: length 5 for the directions (1,0), (0,1), (1,1), (1,2), (2,1), (1,3), (3,1),
  (1,4), (4,1), (3,4), (4,3), length 4 for (2,3), (3,2); arc cubes: complete 2-arc banks are the counts 9, 11,
  13, 15 at orders 5 and 6 (8 and 16 affine), no complete 3-arc bank among the 606,480 cubes of order 7 (nearest
  9, 11, 13, 15, 17, 19, 23, 27; 1,012 cubes permute the odd residues mod 16), no 4-arc bank at orders 5, 6.
  The comparison with Price's Fibonacci-box tree (arXiv:0809.4324; CITED) is an exact identity for one generator
  and an ANALOGY beyond it. NOT CLAIMED: Chan-Pak (OPEN); banks longer than their block count.
  Not independently audited (owner directive 2026-11-13: ideas over audits).
session: mac-mini-2026-11-15-monoid
date: 2026-11-15
source: 05-knowledge/results/chanpak_seed_monoid_line_banks_20261115/README.md
scripts:
  - 05-knowledge/results/chanpak_seed_monoid_line_banks_20261115/scripts/seed_monoid.py, cross_bank.c,
    line_bank.py, line_bank_check.py, arc_bank.c, master_formula_check.py
related:
  - HYP-9842, THM-5282, THM-5272, HYP-9832, HYP-9822
  - 05-knowledge/results/trees5_tournaments4_erdos_sos_20261005.md
  - 05-knowledge/results/recursive_interfaces_20261016b.md
---

# THM-5302 - the seed monoid

Source: [the note](../../05-knowledge/results/chanpak_seed_monoid_line_banks_20261115/README.md).

Notation of THM-5282. `R = Z[e]/(e^2 - 1)`; for a path `Q` of the block `C` the weight is `H(C - V(Q)) e^|Q|`.

**(A)** `M(C, X) = [[h, b, c], [0, G, a], [0, 0, h]]` with `h = H(C)`, `G` = sum over all paths (the empty one
included), `a` = sum over the paths ending at a non-contact, `b` = sum over the paths starting at a contact,
`c` = sum over the paths from a contact to a different non-contact. For an ordinal sum,
`M(C1 -> C2, X1 u X2) = M(C1, X1) M(C2, X2)`, and for every chain of blocks and tape vertices the count with the
head is `h + 2 [c]_even` of the product.

**(B)** One vertex: `M_0 = [[1,0,0],[0,1+e,e],[0,0,1]]` (non-contact), `M_1 = [[1,e,0],[0,1+e,0],[0,0,1]]`
(contact). At `e = 1`: `(b, c) -> (2b, c + b)`, `(2b + 1, c)`; at `e = -1` the corner counts the factors `10`.
For a one-apex word with binary value `b`: `H = 1 + sum over the zero bits k of floor(b / 2^(k+1)) + #(10)`.

**(C)** With `u1 = z1 = 0`: `H = h(1 + 2N_up + 2N_low) + 2[c + a alpha + b beta + 2G alpha beta]` (`a, b, G` the
values at `e = 1`). If seeds `i` share `(h, G)` and `a_i = a_0 + p t_i`, `b_i = b_0 - q t_i`, then for every
exterior with `p alpha - q beta = rho` the counts differ by `2(c_i - c_0 + t_i rho)`.

**(D)** See the status block and `out/` of the note.
