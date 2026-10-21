---
id: THM-4741
title: "Presentation counts of metagraph nodes: the rotation inequality, the pendant nodes of the fundamental black graph, the strong presentation-count spectrum, and the H-descent refutation. (A4) Rotation inequality: for a strong tournament T and any Hamiltonian cycle C, the n rotations of C realise exactly n/|Z_C| tilings (Z_C = the cyclic shifts of C that are automorphisms), all closed, so tc(T) = H/|Aut| >= n/|Z_C| (a corollary of THM-4729 L2's closure bound); tc(T) < n forces a nontrivial cyclic automorphism along every Hamiltonian cycle, and a rigid strong tournament has tc = H >= n. (A5) The pendant nodes of THM-4731's converse-merged fundamental black graph B are exactly the converse pairs of non-palindromic stacks over {K1, C3} ((a(n) - p(n))/2 of them, a = A000930, p = pal13) together with the self-converse classes with tc - bm = 2, which at n <= 9 are the two S4-centred palindromic stacks (1, S4, 1) and (1, 1, S4, 1, 1): B has 1, 3, 3, 6, 7 pendant nodes at n = 5..9, and its degree-0 nodes are exactly the pure-blue classes; this corrects HYP-9347 (4) and THM-4732 section 1 (MISTAKE-648). (A4') A second route to THM-2477's all-order theorem: a strong tournament with one tiling is a circulant with respect to any Hamiltonian cycle (the rotations share the tiling), and the audited uniform prefix-rerouting theorem supplies a black path in every circulant of odd order >= 7 against the blue rotation tiling; n = 5 gives T5 with tc = 3. (C) A strong class with tc = 3 has, along every Hamiltonian cycle, n/|Z_C| in {1, 2, 3}: a circulant, or |Z_C| = n/3, or n = 2 mod 4 with |Z_C| = n/2; at n <= 9 the only strong class with tc = 3 is T5 and the tc = 3 classes are exactly the stacks over {K1, C3, T5} with one T5 (1, 2, 3, 6, 11). (S) The strong presentation-count spectrum at n = 3..9 has minima 1, 5, 3, 9, 9, 17, 25 and the strong Hamiltonian-path minima 3, 5, 9, 15, 25, 45, 75 are Busch's values 5^((n-1)/3), 3 . 5^((n-3)/3), 9 . 5^((n-5)/3) by residue; the fibre sizes 7 and 21 (THM-1370's gaps) are absent for every class, rigid or not, at n <= 9. (D) H is not a descent potential on the black graph (unlike c3, black_descent_20261020c): 4, 19, 99, 867, 13258 non-pure-blue classes at n = 5..9 have no black neighbour of smaller H. (N) 50,190 random closed blue lines at n = 11 (all layers) have three-move detours"
status: >
  PROVED: (A4) (corollary of THM-4729 L2, also directly), (A5) (from THM-2477 all-order + THM-2444's rigidity lemma +
  THM-4731 (3); the identification of the self-converse pendant classes at n <= 9 is FINITE-EXACT), (A4') (second
  route; uses pure_blue_atom_classification_20261020's uniform prefix-rerouting theorem, audited and independent of path
  rigidity; THM-4730 (B)'s explicit paths are FINITE-EXACT to n = 23 only), (C) the index lemma.
  CITED (prior art, found by grep after the pre-audit checkpoint, MISTAKE-647): THM-2477 (as superseded on 2026-10-20 by
  merged_metagraph_finite_pure_atoms_20261020 section 3: the only strong tournament with H = |Aut| is C3 at every order,
  by the one-orbit closing property against Benhocine-Wojda's Hamiltonian-bypass theorem; the rigid classes are the
  ordered {1, C3} stacks with generating function 1/(1 - z - z^3); its successor-automorphism step is shared by (A4'));
  THM-2453 (palindromic Narayana law pal13(n)); THM-2450 section 2 (stack multiplicativity); THM-1370 (the
  Hamiltonian-path spectrum omits 7 and 21 for all n: multiplicativity, the insertion lemma and the Busch floor, the
  same proof this session re-derived; attainment of every other odd number is a conjecture there); HYP-9347 (4) (the
  first, erroneous pendant statement); black_descent_20261020c (the cyclic-triangle descent theorem); Busch 2006
  (Theorem 1 and Corollary 1: h(n) = 3 beta^(n-3), beta^(n-1), 9 beta^(n-5) for n = 0, 1, 2 mod 3, beta = 5^(1/3));
  MathOverflow 232751 via Grinberg's Math 530 notes (Lecture 12, Remark 1.1.11: 7 and 21 excluded, every other odd
  number to 80555 attained, higher unknown); THM-643, THM-849 section 1, THM-4729 L2, LEM-003, Camion.
  FINITE-EXACT: the class files of tiling_line_metagraph_n9_20261020.cpp (closed and minHblack columns, S48) at
  n = 5..9 (presentation_counts_20261020.out): tc = 1 counts 4, 6, 9, 13, 19 (= A000930) with H = 3^j; non-self-converse
  tc = 1 classes 2, 4, 6, 10, 14 (converse pairs: 1, 2, 3, 5, 7 nodes) and the two tc - bm = 2 self-converse classes;
  tc = 3 counts 1, 2, 3, 6, 11; strong spectra and minima (178,133 strong classes at n = 9 = THM-4729's value); 7 and 21
  absent from the tc and H spectra, every other odd value below 2883 present at n = 9; sink counts.
  REFUTED: (D).
  NUMERICAL: (N), blue_line_black_detours_20261020_n11_rand.out (random-sample mode, with replacement); the full
  order-11 census was running at close-out.
  CONCURRENT: mac-mini-2026-10-21-electric (THM-4739 / HYP-9351, Kirchhoff/Hodge laws, published 20:44; THM-4740
  ceded) and the Codex packet metagraph_network_invariants_20261021 (20:47, reflection-paired network B + 2Q,
  spanning-tree polynomials, critical group); lane split by message; nothing electrical is claimed here.
  AUDITED: independent reproduction of every data column at n <= 7 by a full canonical-form census, of every spectrum
  number at n <= 9, of the rotation statements on all Hamiltonian cycles at n <= 7, of the sampler and the stride
  partition; eight majors corrected before publication (note section 10; MISTAKE-648).
  OPEN: tc >= 5 for every strong tournament other than K1, C3, T5 (HYP-9352 (1)); fibre sizes 7, 21 for non-rigid
  classes (HYP-9352 (2)); the self-converse pendant classes beyond n = 9 (HYP-9352 (4)).
session: opus-2026-10-20-S48
source: 05-knowledge/results/presentation_counts_single_view_classes_rotation_inequality_20261020.md
depends_on:
  - 01-canon/definitions.md
  - 01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md
  - 01-canon/theorems/THM-643-blue-black-line-parity-structure.md
  - 01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md
related:
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 01-canon/theorems/THM-2444-pure-blue-count-refutation-and-rigid-class-law.md
  - 01-canon/theorems/THM-2450-rigid-self-converse-classes-are-cyclic-ternary-towers.md
  - 01-canon/theorems/THM-2453-palindromic-narayana-law-for-the-rigid-stratum.md
  - 01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md
  - 01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md
  - 05-knowledge/results/merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md
  - 05-knowledge/results/pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md
  - 05-knowledge/results/black_descent_20261020c/cycles.md
  - 05-knowledge/hypotheses/HYP-9347-three-move-law-layer-theorems-fundamental-black-graph-and-the-klein-dictionary-beyond-order-4.md
  - 05-knowledge/hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md
scripts:
  - 04-computation/experiments/tiling_line_metagraph_n9_20261020.cpp (S48 columns; + _classes_n5..n8.tsv, _classes_n9.tsv.gz, _n9_s48.out)
  - 04-computation/experiments/presentation_counts_20261020.py (+ .out)
  - 04-computation/experiments/blue_line_black_detours_20261020.cpp (random-sample and stride modes; + _n11_rand.out)
---

# THM-4741: presentation counts of metagraph nodes

Conventions: THM-643's tilings (presentations) of a class, `tc = H/|Aut|`, the number of `Aut`-orbits of Hamiltonian paths (LEM-003) and so of distinct views (THM-849 section 1); stacks `(B_1, ..., B_m)` of strong components with all arcs downward, over which `H`, `|Aut|` and `tc` multiply (THM-2450 section 2); THM-4731's converse-merged fundamental black graph `B` (degree `tc` at a merged non-self-converse node, `(tc - bm)/2` at a self-converse class); `S4` the strong tournament of order 4. Proofs and tables in the [session note](../../05-knowledge/results/presentation_counts_single_view_classes_rotation_inequality_20261020.md).

**Prior art used.** THM-2477 (updated 2026-10-20): the only strong tournament with `H = |Aut|` is `C3`; hence `tc = 1` iff the tournament is an ordered stack over `{K1, C3}`, counted by `1/(1 - z - z^3)` (A000930; `4, 6, 9, 13, 19` at `n = 5..9`, each with `H = |Aut| = 3^j`). THM-2453: the self-converse ones are the palindromic stacks, `p(n)` of them. THM-1370: no tournament has exactly 7 or 21 Hamiltonian paths. `black_descent_20261020c`: the cyclic-triangle count descends along a black line from every strong black-active node.

**(A4) Rotation inequality.** Let `T` be strong, `C = (v_0, ..., v_{n-1})` a Hamiltonian cycle (Camion), `Z_C` the group of cyclic shifts of `C` that are automorphisms of `T`. Rotations `i` and `j` of `C` have the same tiling iff the positional map `v_{i+a} -> v_{j+a}` is an automorphism, i.e. the shift by `j - i` lies in `Z_C`; so the `n` rotations realise exactly `n/|Z_C|` tilings, all closed (the apex pair of a rotation carries a cycle arc). Hence `tc(T) >= n/|Z_C|`; if `Aut(T)` is trivial the rotations are `n` distinct tilings and `tc = H >= n`. Equivalently, by THM-4729 L2 the closed presentations number `n HC/|Aut| <= tc`, and the `Aut`-orbit of `C` consists of `|Aut|/|Z_C|` Hamiltonian cycles.

**(A4') Second route to THM-2477's theorem.** If a strong `T` has a single tiling, all rotations of `C` share it, so `v_a -> v_b` iff `b - a in S := {d : v_0 -> v_d}` (compare the rotations at `0` and at `-a`): `T` is the circulant with symbol `S`, of odd order (this successor-automorphism step is also the first step of THM-2477's superseded proof). `n = 3` is `C3`; `n = 5` is `T5`, with `tc = 3`; for odd `n >= 7` the rotation tiling is blue (`sigma`-fixed, the tile orientation depending only on `b - a`) while the uniform prefix-rerouting theorem of `pure_blue_atom_classification_20261020` exhibits a black Hamiltonian path, a second tiling. (THM-2477's own last step is Benhocine-Wojda's bypass theorem.)

**(A5) Pendant nodes.** A merged node of `B` has degree 1 iff it is a non-self-converse class with `tc = 1` (by THM-2477 a non-palindromic `{K1, C3}`-stack; its converse is the reversed stack and the pair is one node) or a self-converse class with `tc - bm = 2`; degree 0 iff `tc = bm`, i.e. pure-blue (the palindromic `{K1, C3}`-stacks are pure-blue by THM-2444's rigidity lemma). The census identifies the `tc - bm = 2` self-converse classes at `n <= 9` as `(1, S4, 1)` and `(1, 1, S4, 1, 1)`, `(H, |Aut|, tc, bm) = (5, 1, 5, 3)`; so the pendant counts are `(a(n) - p(n))/2 + [n in {6, 8}] = 1, 3, 3, 6, 7` at `n = 5..9`.

**(C) `tc = 3`.** By (A4), a strong class with `tc = 3` has `n/|Z_C| <= 3` dividing `n` with `|Z_C|` odd, so `n/|Z_C| in {1, 2, 3}` with `2` only for `n = 2 mod 4` and `|Z_C| = n/2` (the index-2 case occurs at `n = 6` for the strong self-converse classes with `(H, |Aut|, tc) = (45, 3, 15)`). The census shows `tc = 3` for a strong class only at `n = 5` (`T5`), so by multiplicativity the `tc = 3` classes at `n <= 9` are exactly the stacks over `{K1, C3, T5}` with one `T5` (`1, 2, 3, 6, 11`; the `(H, |Aut|) = (45, 15)` members at `n = 8, 9` contain a `C3` as well).

**(S) Strong spectrum.** Strong iff `closed > 0` (THM-4729 L2):

| n | strong classes | smallest strong tc (count) | min strong H | tc = 1 / pendant nodes / tc = 3 classes |
|---|---|---|---|---|
| 5 | 6 | 3 (1), 5 (1), 9 (2), 11 (1), 13 (1) | 9 | 4 / 1 / 1 |
| 6 | 35 | 9 (2), 11 (2), 15 (4), 17 (2), 19 (2) | 15 | 6 / 3 / 2 |
| 7 | 353 | 9 (1), 11 (1), 15 (2), 17 (3), 19 (2) | 25 | 9 / 3 / 3 |
| 8 | 6008 | 17 (1), 23 (4), 25 (2), 27 (4), 29 (2) | 45 | 13 / 6 / 6 |
| 9 | 178133 | 25 (1), 33 (5), 37 (2), 39 (1), 45 (10) | 75 | 19 / 7 / 11 |

The strong `H` minima `3, 5, 9, 15, 25, 45, 75` (`n = 3..9`) are Busch's `5^((n-1)/3)`, `3 . 5^((n-3)/3)`, `9 . 5^((n-5)/3)` by residue. Among the odd values below 2883, exactly 7 and 21 are missing from the presentation-count spectrum at `n = 9` (the spectrum is sparse between 2883 and its maximum 3299; at `n = 7, 8` the same two gaps are the only ones below 63 and 611), and 7 and 21 are missing from the `H` spectrum at every `n = 5..9` (THM-1370); for a rigid node `tc = H`, so rigid nodes never have fibre size 7 or 21; for non-rigid nodes this is HYP-9352 (2).

**(D)** With `minHblack` the minimum `H` over the black neighbours of a class, `4, 19, 99, 867, 13258` non-pure-blue classes at `n = 5..9` have no black neighbour of smaller `H` (at `n = 9`, 16 of them have a neighbour of equal `H` and the rest only larger ones), with sink `H` values spread through the range; `H` carries no gradient structure on the black graph, in contrast with `c3`.

**(N)** Eight seeded samples of 12,500 sigma-orbit codes at `n = 11` (with replacement, about 80 duplicate draws) gave 50,190 closed blue lines across all nine layers (`181, 1555, 5446, 11060, 13719, 10960, 5446, 1628, 195`), every one with a literal three-black-move detour, none pure-blue-touching (59 s per sample, 9.4 ms per line).
