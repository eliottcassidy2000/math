---
id: THM-4741
title: "Presentation counts of metagraph nodes. (A1) The presentation count tc = H/|Aut| of a class is the number of Aut-orbits of Hamiltonian paths, i.e. of distinct tilings (views). (A2) tc is multiplicative over the strong-component stack (H and |Aut| are). (A3) tc = 1 iff the tournament is a stack over {K1, C3}; the number of such classes of order n is the number of compositions of n into parts 1 and 3 (A000930: 4, 6, 9, 13, 19 at n = 5..9), each with H = |Aut| = 3^(number of C3 blocks). (A4) Rotation inequality: for strong T and any Hamiltonian cycle C, the n rotations of C realise exactly n/|Z_C| tilings (Z_C = cyclic shifts of C that are automorphisms), all closed, so tc(T) >= n/|Z_C|; tc(T) < n forces a nontrivial cyclic automorphism along every Hamiltonian cycle, and a rigid strong tournament has tc = H >= n. (A5) The pendant nodes of the fundamental black graph B (THM-4731) are exactly the non-palindromic {K1, C3}-stacks, a(n) - p(n) of them (14 at n = 9), and the palindromic ones are pure-blue (p(n) = a(n/2) or a((n-1)/2) + a((n-3)/2)). (B) No tournament has exactly 7 or 21 Hamiltonian paths (stack product + insertion monotonicity H(T) >= H(T - v) + Moon's strong subtournament theorem + the census minima h(3..7) = 3, 5, 9, 15, 25), so a rigid node never has fibre size 7 or 21; at n <= 9 no node of any kind does. (C) A strong class with tc = 3 is circulant with respect to every Hamiltonian cycle or has |Z_C| = n/3; at n <= 9 the only strong class with tc = 3 is T5, and the tc = 3 classes are exactly the stacks over {K1, C3, T5} with one T5 (1, 2, 3, 6, 11 at n = 5..9). (D) H is not a descent potential on the black graph: 4, 19, 99, 867, 13258 non-pure-blue classes at n = 5..9 have no black neighbour of smaller H"
status: >
  PROVED: (A1)-(A5) (with THM-4730 (B) for the black path in every odd circulant of order >= 7, THM-2444's rigidity
  lemma for the self-converse single-presentation classes, THM-4731 (3) for the degrees of B, Camion and LEM-003);
  (B) with Moon 1966 (every strong tournament on n >= 4 vertices has a strong subtournament on n - 1 vertices) CITED
  and the census minima FINITE-EXACT; (C) the reduction lemma.
  CITED (prior art): the exclusion of 7 and 21 and the realisation of every other odd number up to 80555 is
  MathOverflow 232751 (bof, Gordon Royle) as reported in D. Grinberg's Rédei-Berge notes (KTH 2024), which record as
  open whether 7 and 21 are the only exceptions; A. H. Busch, Electron. J. Combin. 13 (2006) N3 (a strong tournament
  of order n has at least 5^((n-1)/3) Hamiltonian paths, sharp for n = 1 mod 3; the census minima 3, 5, 9, 15, 25,
  45, 75 are 5^((n-1)/3), 3 . 5^((n-3)/3), 9 . 5^((n-5)/3) by residue, an independent route to h(n) > 21 for n >= 8);
  THM-849 section 1 (observer paths); THM-643 (tc odd); THM-4729 L2
  (closed presentations > 0 iff strong); Rédei.
  FINITE-EXACT: the class files of tiling_line_metagraph_n9_20261020.cpp (closed and minHblack columns, S48) at
  n = 5..9: tc = 1 counts 4, 6, 9, 13, 19 with H = 3^j; pendant nodes 14 at n = 9; tc = 3 counts 1, 2, 3, 6, 11; the
  strong presentation-count spectra (minima 3, 9, 9, 17, 25 at n = 5..9) and strong Hamiltonian-path minima
  9, 15, 25, 45, 75 (178,133 strong classes at n = 9, THM-4729's value); 7 and 21 absent from the tc and H spectra at every n <= 9; the sink counts.
  REFUTED: the heuristic (D).
  NUMERICAL: the n = 11 random sample of the three-move law (HYP-9352 (3)).
  CONCURRENT: mac-mini-2026-10-21-electric reserved THM-4739/4740 / HYP-9351 at 20:18 and published THM-4739 / HYP-9351
  at 20:44 (Hodge decomposition of the arc flow with the scores as potentials and the cyclic-triangle 2-form as curl,
  Kirchhoff's law of the merged metagraph for exact-flux invariants, second-order determination, consecutive-start-port
  ledgers; THM-4740 ceded); the Codex packet metagraph_network_invariants_20261021 (20:47) proved componentwise laws of
  the reflection-paired network B + 2Q (spanning-tree polynomials, critical group, resistances) on the n <= 7 atlas.
  Lane split by message; nothing electrical is claimed here and neither packet is modified.
  OPEN: tc >= 5 for every strong tournament other than K1, C3, T5 (HYP-9352 (1)); fibre sizes 7, 21 for non-rigid
  classes (HYP-9352 (2)).
session: opus-2026-10-20-S48
source: 05-knowledge/results/presentation_counts_single_view_classes_rotation_inequality_20261020.md
depends_on:
  - 01-canon/definitions.md
  - 01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md
  - 01-canon/theorems/THM-643-blue-black-line-parity-structure.md
  - 01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md
related:
  - 01-canon/theorems/THM-2444-pure-blue-count-refutation-and-rigid-class-law.md
  - 01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md
  - 01-canon/theorems/THM-849-self-line-census-n8-klein-orbits.md
  - 05-knowledge/hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md
scripts:
  - 04-computation/experiments/tiling_line_metagraph_n9_20261020.cpp (S48 columns; + _classes_n5..n8.tsv, _classes_n9.tsv.gz)
  - 04-computation/experiments/presentation_counts_20261020.py (+ .out)
  - 04-computation/experiments/blue_line_black_detours_20261020.cpp (random-sample mode; + _n11_rand.out)
---

# THM-4741: presentation counts, single-view classes, the rotation inequality, and the forbidden fibre sizes

Conventions: THM-643's tilings (presentations) of a class, `tc = H/|Aut|`; stacks `(B_1, ..., B_m)` of strong components with all arcs downward; THM-4731's fundamental black graph `B` (degree `tc` at a merged non-self-converse node, `(tc - bm)/2` at a self-converse node). Proofs in the [session note](../../05-knowledge/results/presentation_counts_single_view_classes_rotation_inequality_20261020.md), sections 2-5.

**(A1)** Two Hamiltonian paths have the same tiling iff an automorphism carries one to the other (the positional map preserves all arcs); `Aut` acts freely (LEM-003), so `tc` counts the orbits, i.e. the distinct views.

**(A2)** A Hamiltonian path exhausts each strong component before entering the next, so `H = prod H(B_i)`; automorphisms preserve the condensation order, so `|Aut| = prod |Aut(B_i)|`; hence `tc = prod tc(B_i)`.

**(A3)** If a strong `T` has one tiling, the rotations of a Hamiltonian cycle `C` force `v_a -> v_b` iff `b - a in S = {d : v_0 -> v_d}`, so `T` is a circulant of odd order; `n = 3` gives `C3`, `n = 5` gives `T5` (`tc = 3`), and for odd `n >= 7` THM-4730 (B)'s black path has a tiling different from the blue rotation tiling. With (A2), `tc = 1` iff every block is `K1` or `C3`; the count is the number of compositions of `n` into parts 1 and 3.

**(A4)** Rotations `i`, `j` of `C` share a tiling iff the shift by `j - i` is an automorphism; the rotations are closed (apex = cycle arc).

**(A5)** Degree 1 in `B` with `tc = 1` needs a non-self-converse class; a self-converse one has `bm = 1 = tc` (THM-2444's rigidity lemma: pure-blue); a `{K1, C3}`-stack is self-converse iff palindromic.

**(B)** `H(T) >= H(T - v)` (insert `v` at the front, at the end, or at a sign change); Moon: a strong `T` on `n >= 4` vertices has a strong `(n-1)`-subtournament, so the strong minimum `h(n)` is non-decreasing; `h(3..7) = 3, 5, 9, 15, 25` and the strong counts at `n = 5, 6, 7` avoid 21 (census); `H = 7` or `21` would need a strong block with `H = 7` or `21`, or blocks `3 x 7`.

**(C)** By (A4) `n/|Z_C| in {1, 3}` when `tc = 3` (`|Z_C|` is odd).

| n | tc = 1 (A000930) | pendant nodes of B | tc = 3 (one T5) | strong classes | min strong tc | min strong H | odd values absent from the tc / H spectra |
|---|---|---|---|---|---|---|---|
| 5 | 4 | 2 | 1 | 6 | 3 | 9 | 7 (spectrum ends at 13) / 7, 17, ... |
| 6 | 6 | 4 | 2 | 35 | 9 | 15 | 7, 21, 35, ... / 7, 21, 35, ... |
| 7 | 9 | 6 | 3 | 353 | 9 | 25 | 7, 21 (and 63 for H) |
| 8 | 13 | 10 | 6 | 6008 | 17 | 45 | 7, 21 |
| 9 | 19 | 14 | 11 | 178133 | 25 | 75 | 7, 21 |
