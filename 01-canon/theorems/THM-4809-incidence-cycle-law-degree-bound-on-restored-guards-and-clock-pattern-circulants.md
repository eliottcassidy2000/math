---
id: THM-4809
title: "Three quantities of a finite response protocol, and clock patterns. (A) and (B) are session72's statements (published at 31002464ae minutes before this filing), restated with proofs and extended by the order-7 paired-cell row, the guard ranks at orders 6-7 and the cover-depth census; (C) is this session's. (A) Incidence-cycle law: if every state of a finite protocol has two terminal responses a(x) in A, b(x) in B, with E the attained pairs, V the attained responses and c the components of the bipartite graph (V, E), then a joint law on the paired cells is observed through its marginals with rank |V| - c, and exactly |E| - |V| + c directions, the cycle space, are invisible; for the anonymous reset cells ([x], [tau x]) of the marked carrier X_n this gives (E, V, c, rank, lost) = (7, 8, 1, 7, 0), (50, 24, 1, 23, 27), (732, 112, 1, 111, 621), (24654, 912, 1, 911, 23743) at n = 4, 5, 6, 7. (B) Degree bound: a named Hamilton-path guard on X_n is a subcube indicator fixing at most n - 1 of the M = (n-1)(n-2)/2 free pairs, so k-fold products are Boolean polynomials of degree <= k(n - 1); full recovery of the 2^M presentation-law coordinates by k-fold restored tests needs k >= ceil((n-2)/2); pairs span at orders 4, 5, 6 (ranks 8, 64, 1024) and cannot from order 7 (at most 32647 of 32768); single guards have rational rank 8, 44, 273 at n = 4, 5, 6. (C) Clock patterns: the pattern S_m of a sequence modulo m (its nonzero residues) defines the circulant Clk(Z/m, S_m), an oriented graph iff S_m and -S_m are disjoint and a tournament iff moreover |S_m| = (m-1)/2; the squares give a tournament iff p = 3 mod 4 (Paley), the powers a^k iff -1 is not in <a> and <a> has index 2."
status: >
  CITED-CONCURRENT, restated: (A) and (B) are session72's reset-action cycle theorem and Boolean-degree obstruction
  (05-knowledge/results/session72_even_graph_unifying_carriers_20261026/session72-response/REPORT.md, with its sharp-depth
  conjecture d_n = ceil((n-2)/2) verified through order 6); the proofs here are the classical ones (the marginal map is the
  transpose of a vertex-edge incidence map, of rank |V| - c with kernel the cycle space; monomial counting on the cube).
  PROVED: (C) THM-4532's circulant criterion with the order of -1 in <a>.
  FINITE-EXACT: guard_ranks_paired_cells_20261026.cpp (rational ranks modulo two primes at n <= 6, exact over Q with
  fractions at n <= 5 in guard_rank_restored_pairs_20261026.py, GF(2) ranks, the paired-cell graphs at n <= 7, all
  connected), restore_depth_20261026.cpp (cover depths: least cover depth among coverable states equals the degree
  bound 1, 2, 2, 3 at n = 4..7 and 3 on 10,000 sampled states at n = 8; 1, 13, 68, 492 states with a dark free pair are
  never a common cube), clock_arrows_circulant_20261026.py (census to modulus 120: squares at the primes 3 mod 4,
  powers of 2 at 7, 23, 47, 71, 79, 103, powers of 3 at 3, 11, 23, 47, 59, 71, 83, 107, triangular numbers at 3 and 5;
  cubes, Fibonacci numbers, primes and odd Collatz successors never oriented).
  CITED: the unnamed-carrier packet (U2 reset protocol, the 50 cells at n = 5, the three-layer proposal, the 44/5
  ledger), session72 (the announced incidence-cycle reading, the positive cone), session33 and boundary_support
  Proposition 0 (the same rank law in converse-merged coordinates, r - c - b and 2r - c - b), session70 (A2-A4),
  THM-4806 / THM-4807 (dark arcs), THM-4532 (clocks), THM-4644 and THM-126 (circulant tournaments), THM-1600
  (derangements 9, 44, 265: the guard ranks 8, 44, 273 are not derangement counts), session71 (the 44 square and 5
  triangle cut types; the angle charge; the field obstruction), the concurrent branch clock_even_guarded_composition
  (the index / successor / iteration arrow trichotomy; not on main at filing).
  NUMEROLOGY: every 44 / 5 coincidence except the Gaussian ray (2 + i)^3 = 2 + 11i, (2 + 11i)^2 = -117 + 44i.
  AUDITED: pending (note, section 12).
  OPEN: HYP-9419.
session: opus-2026-10-26-S62
source: 05-knowledge/results/three_quantities_fortyfour_five_dudeney_fermat_clock_arrows_20261026.md
depends_on:
  - 01-canon/definitions.md
  - 01-canon/theorems/THM-4532-tournament-clocks-two-sheet-clocks-all-odd-tournaments-and-the-syracuse-clock.md
  - 05-knowledge/results/unnamed_carrier_dissection_20261026/CONNECTIONS.md
  - 05-knowledge/results/session72_even_graph_unifying_carriers_20261026/README.md
  - 01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md
related:
  - 05-knowledge/hypotheses/HYP-9419-unifying-candidates-three-layer-composition-rank-laws-dissection-obstructions-pattern-clocks.md
  - 05-knowledge/results/session71_dissection_norm_alignment_20261026/README.md
  - 05-knowledge/results/session70_attachment_guard_alignment_20261026/README.md
  - 01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md
  - 01-canon/theorems/THM-4806-named-legality-language-single-pair-lemma-goldbach-free-shapes-and-row-companions-of-even-perfect-numbers.md
---

# THM-4809: three quantities, and clock patterns

## Theorem A (incidence-cycle law; session72, restated)

Let a finite protocol attach to every state `x` terminal responses `a(x) in A` and `b(x) in B`; let `E` be the set of attained pairs `(a(x), b(x))`, `V` the attained responses, `c` the number of components of the bipartite graph `(V, E)`. A joint law `mu` on the paired cells (a vector indexed by `E`) is observed through its marginals `(mu_A, mu_B)`; this linear map has rank `|V| - c` and kernel the cycle space of `(V, E)`, of dimension `|E| - |V| + c`.

*Proof.* The marginal map is the transpose of the vertex-edge incidence matrix of the bipartite graph; its rank is `|V| - c` and its kernel on edge space is the cycle space. QED

For the anonymous reset cells `([x], [tau x])` of `X_n` (`tau` complements every free pair; `[.]` the isomorphism class): orders 4, 5, 6, 7 give `E = 7, 50, 732, 24654`, `V = 8, 24, 112, 912` (twice the number of tournament classes), `c = 1`, rank `7, 23, 111, 911`, lost directions `0, 27, 621, 23743`. The order-5 row is the pasted "50 / 23 / 27".

## Theorem B (degree bound on restored named guards; session72's lower bound, restated)

On `X_n` a named guard `g_p` (`p in S_n`) is the indicator of the subcube of states whose free pairs along `p`, at most `n - 1` of them, are oriented as `p` requires. A product `g_(p_1) ... g_(p_k)` is therefore the indicator of a subcube fixing at most `k(n - 1)` free pairs, i.e. a Boolean polynomial of degree `<= k(n - 1)` in the `M` coordinates, and the span of all `k`-fold products has dimension at most `sum_(j <= k(n-1)) binom(M, j)`. Full recovery of the `2^M` presentation-law coordinates by `k`-fold restored tests needs `k(n - 1) >= M`, i.e. `k >= ceil((n - 2)/2)`.

*Proof.* Monomial counting: `[x_i = 0] = 1 - x_i`, so an indicator fixing `d` coordinates has degree `d`. QED

FINITE-EXACT: single guards have rational rank 8, 44, 273 at `n = 4, 5, 6` (47 distinct guards among 120 permutations at order 5; GF(2) ranks 8, 43, 272, and 1934 at order 7); pair products span the whole space at orders 4, 5, 6 (ranks 8, 64, 1024) and cannot at order 7 (at most `32768 - 105 - 15 - 1 = 32647`). The point indicator of a state is itself a product of its guards iff some `k` of its Hamiltonian paths cover all free pairs; states with a dark free pair (`1, 13, 68, 492` at orders 4-7) never are, and the least cover depth among the others equals the bound `1, 2, 2, 3` (order 8: 3 on a sample of 10,000). So the full rank of pair products at orders 5 and 6 is a linear-combination phenomenon (only 21 of the 64 order-5 states are two-path cubes).

## Theorem C (clock patterns)

For a sequence `(a_k)` and a modulus `m`, let `S_m` be the set of nonzero residues attained. The circulant `Clk(Z/m, S_m)` (arcs `i -> i + s`, `s in S_m`) is an oriented graph iff `S_m` and `-S_m` are disjoint, and a tournament iff moreover `|S_m| = (m - 1)/2` (THM-4532's criterion). For an odd prime `p`: the squares give a tournament iff `p = 3 mod 4`, the Paley tournament; the powers `a^k` give an oriented graph iff `-1` is not in `<a>` and a tournament iff moreover `<a>` has index 2, i.e. `a` is a quadratic residue of order `(p - 1)/2` (for `a = 2`: `p = 7 mod 8` with `ord_p(2) = (p - 1)/2`).

*Proof.* `-1 in <a>` iff the pattern meets its negative; index 2 iff half-set; squares: `-1` is a square iff `p = 1 mod 4`. QED

Census to modulus 120: note, section 7.
