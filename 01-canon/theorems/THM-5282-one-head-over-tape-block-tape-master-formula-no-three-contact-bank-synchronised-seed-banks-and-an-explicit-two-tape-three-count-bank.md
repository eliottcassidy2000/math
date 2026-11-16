---
id: THM-5282
title: "One head over tape-block-tape. (A) Master formula: H = h(1 + 2N_up + 2N_low) + 2[cut_C(X) + a_o(alpha+u1) + a_e alpha + b_o(beta+z1) + b_e beta + G_e((alpha+u1)(beta+z1) + alpha beta) + G_o((alpha+u1) beta + alpha(beta+z1))]. (B) Seeds (block, pattern) with equal invariants (h, G_e, G_o, a_o, a_e, b_o, b_e) have counts differing by a constant for every pair of tape words. (C) For all A, B >= 0 the numbers K+90, K+92, K+94 are Hamiltonian-path counts, K = 432AB + 174A + 303B + 45 - 30(A mod 2) - 15(B mod 2) - 90(J(A) + J(B)). (D) Composition: for C1 -> C2, h = h1 h2, G = G1 (x) G2 by parity, cut = h2 cut1 + h1 cut2 + b_e a_e + b_o a_o. (E) Finite-exact: no complete three-contact bank for any block of order at most 8 and any tapes"
status: >
  PROVED: (A) from the cut form (THM-5272 (A)) by splitting a path into its runs in the upper tape, the block and
  the lower tape; (B) immediate from (A); (C) from (A) for three explicit seeds on 6-vertex blocks with
  h = 45, G_e = G_o = 216, a = (72, 57), b = (144, 114), cut = 45, 46, 47; (D) by the same splitting.
  Brute-force checks: (A) 60 random instances, (C) 40 random tape pairs (tournaments of order up to 15), (D) 40.
  FINITE-EXACT: (E) for the 7,410 isomorphism classes of blocks of order 3..8, every three contacts in the block,
  every block context and all tape data (exact linear systems in the tape parameters): no complete 3-bank, and
  complete 2-banks only at isolated tape data (blocks of order <= 6). Longest synchronised banks through block
  order 8 (1,823,696 seeds): 6 (upper tape, height 1 free), 4 (upper tape), 3 (two tapes, either convention).
  The bank (C) reaches 11.57 / 14.27 / 16.87 per cent of the odd numbers below 10^5 / 10^6 / 10^7 (3.9 per cent
  with the upper tape alone).
  NOT CLAIMED: Chan-Pak (OPEN); a selector for prescribed targets; density statements; three-contact banks on
  bases with several separated blocks. One-tape reader banks of three and five counts and the capacity bound are
  variable_reader_coverage_20261116 (another session); the two-tape banks and (D) in this form are added here.
  Not independently audited (owner directive 2026-11-13: ideas over audits).
session: mac-mini-2026-11-15-tapes
date: 2026-11-15
source: 05-knowledge/results/chanpak_two_tape_banks_jacobian_anatomy_20261115/README.md
scripts:
  - 05-knowledge/results/chanpak_two_tape_banks_jacobian_anatomy_20261115/scripts/master_formula_check.py,
    cube_bank_tapes.py, sync_bank.c, cross_bank.c, two_tape_bank_check.py, bank_coverage.py
related:
  - HYP-9832, THM-5272, HYP-9822, THM-5262
  - 05-knowledge/results/variable_reader_coverage_20261116/README.md
  - 05-knowledge/results/physical_radix_carry_20261116/README.md
---

# THM-5282 - one head over tape-block-tape

Source: [the note](../../05-knowledge/results/chanpak_two_tape_banks_jacobian_anatomy_20261115/README.md).

Base `TT_i -> C -> TT_j` (ordinal sum), head `x` with contacts `U` (upper tape), `X` (block), `L` (lower tape);
`Z` = lower non-contacts. Heights and depths are counted from the block.
`alpha = sum_{a in U, height >= 2} 2^(height-2)`, `u1 = [height 1 in U]`,
`beta = sum_{b in Z, depth >= 2} 2^(depth-2)`, `z1 = [depth 1 in Z]`.
For a directed path `Q` of `C` put `w(Q) = H(C - V(Q))`. `lam_par[s]` (resp. `kap_par[s]`) is the sum of `w(Q)`
over the paths ending (resp. starting) at `s` whose number of vertices has that parity; `G_par` the sum over all
paths, the empty one being even; `W(a,b)` the sum over the even paths from `a` to `b`;
`cut_C(X) = sum_{a in X, b notin X} W(a,b)`, `a_par(X) = sum_{s notin X} lam_par[s]`,
`b_par(X) = sum_{s in X} kap_par[s]`; `h = H(C)`; `N_up`, `N_low` the one-apex half-counts of the tape words.

**(A)** `H(B + x) = h(1 + 2N_up + 2N_low) + 2[cut_C(X) + a_o(alpha+u1) + a_e alpha + b_o(beta+z1) + b_e beta
+ G_e((alpha+u1)(beta+z1) + alpha beta) + G_o((alpha+u1) beta + alpha(beta+z1))]`.

**(B)** If two seeds have the same `(h, G_e, G_o, a_o, a_e, b_o, b_e)`, their counts differ by
`2(cut - cut')` for all tapes. With no lower tape only `(h, a_o, a_e)` matters.

**(C)** With the seeds `(22 44 24 49 34 5; 23)`, `(38 28 24 49 33 6; 39)`, `(38 28 24 49 33 6; 23)` (out-masks of
the block; head pattern), `A = 2 alpha + u1`, `B = 2 beta + z1` and `J(n) = sum over pairs of binary digits of n at
distance g of c(g)`, `c(1) = 1`, `c(g) = 2^(g-2)`: the three counts are `K + 90, K + 92, K + 94`,
`K = 432AB + 174A + 303B + 45 - 30(A mod 2) - 15(B mod 2) - 90(J(A) + J(B))`, for all `A, B >= 0`; `K = 3 (mod 6)`.

**(D)** For `C = C1 -> C2` and patterns `X1, X2`: `h = h1 h2`, `G_e = G1e G2e + G1o G2o`,
`G_o = G1e G2o + G1o G2e`, `cut_C(X1 u X2) = h2 cut_{C1}(X1) + h1 cut_{C2}(X2) + b_e(X1) a_e(X2) + b_o(X1) a_o(X2)`.

**(E)** See the status block and `out/cube_bank_tapes.txt`.
