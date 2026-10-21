---
id: HYP-9352
title: "After THM-4741: tc >= 5 for every strong tournament other than K1, C3, T5 (so tc = 3 classes are exactly the one-T5 stacks for all n); fibre sizes 7 and 21 never occur; the three-move law at n = 11 (random sample); growth of the strong spectrum minimum; the pendant nodes of the fundamental black graph (converse pairs of non-palindromic {K1,C3}-stacks plus the self-converse classes with tc - bm = 2) and their count law"
status: "NUMERICAL ((1), (2), (4) exhaustive at n <= 9; (3) sampled at n = 11); OPEN"
session: opus-2026-10-20-S48
source: 05-knowledge/results/presentation_counts_single_view_classes_rotation_inequality_20261020.md
related:
  - 01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md
  - 05-knowledge/hypotheses/HYP-9347-three-move-law-layer-theorems-fundamental-black-graph-and-the-klein-dictionary-beyond-order-4.md
  - 05-knowledge/hypotheses/HYP-9346-self-lines-connectivity-and-growth-laws-of-the-tiling-line-metagraph.md
  - 05-knowledge/results/black_descent_20261020c/cycles.md
---

# HYP-9352: presentation-count laws after THM-4741

1. **Strong spectrum floor.** Every strong tournament other than `K1`, `C3`, `T5` has presentation count `tc >= 5`; equivalently `H(T) > 3 |Aut(T)|` for strong `T` of order `>= 6`. With the stack multiplicativity (THM-2450 section 2) this makes the `tc = 3` classes exactly the stacks over `{K1, C3, T5}` with one `T5` for every `n` (exhaustive at `n <= 9`: `1, 2, 3, 6, 11`). Sharper: the minimum of the strong spectrum `min S(n)` (`1, 5, 3, 9, 9, 17, 25` at `n = 3..9`) is non-decreasing for `n >= 5` and at least `2n - 1` for `n >= 8`. What THM-4741 (C) gives for `tc = 3`: along every Hamiltonian cycle `T` is a circulant, or `|Z_C| = n/3`, or `n = 2 mod 4` with `|Z_C| = n/2`. Hostile: a strong tournament with a large automorphism group and few Hamiltonian paths (the wreath products `C3[C3[...]]` have `|Aut| = 3^((n-1)/2)` but far more paths).

2. **Forbidden fibre sizes.** No node of the line metagraph has `tc = 7` or `tc = 21` at any order (THM-1370 proves it for rigid nodes, where `tc = H`; a non-rigid counterexample needs `H = 7 |Aut|` or `21 |Aut|` with `|Aut| >= 3`, e.g. `H = 35, 49, 63`). Exhaustive at `n <= 9`. Companion (corrected by the audit): at `n = 9` every odd value below 2883 other than 7 and 21 occurs as a presentation count, with gaps only between 2883 and the maximum 3299; conjecture: for every `n` the presentation-count spectrum is an initial segment of the odd numbers minus `{7, 21}` followed by a sparse tail.

3. **Three-move law at `n = 11`** (HYP-9347 (1)): every closed blue line between non-pure-blue classes has a literal three-black-move detour at `n = 11`. Evidence: 50,190 random closed blue lines (drawn with replacement, about 80 duplicates), all layers, all with three-move detours and no exception (THM-4741's note, section 6); the full census (`2^24` lines, about 160,000 core-seconds) was running at close-out.

4. **Pendant nodes of the fundamental black graph** (corrects HYP-9347 (4), which called them "exactly the non-self-converse single-presentation classes"; MISTAKE-648). The merged black graph `B` of THM-4731 has degree `tc` at a merged non-self-converse node and `(tc - bm)/2` at a self-converse node, so its pendant nodes are (i) the converse pairs of non-palindromic `{K1, C3}`-stacks, `(a(n) - p(n))/2` of them (`a` = A000930, `p` the palindromic compositions into parts 1 and 3), and (ii) the self-converse classes with `tc - bm = 2`; at `n <= 9` the latter are exactly the two `S4`-centred palindromic stacks `(1, S4, 1)` (`n = 6`) and `(1, 1, S4, 1, 1)` (`n = 8`), both with `(H, |Aut|, tc, bm) = (5, 1, 5, 3)`, so the pendant counts are `1, 3, 3, 6, 7` at `n = 5..9`. Conjecture: the self-converse pendant nodes are exactly the palindromic stacks with a central `S4` and `{K1, C3}` wings (`tc - bm = 2` because `S4`'s five tilings contain three blue ones), and the degree-`d` nodes of `B` for fixed `d` are eventually a quasi-polynomial in `n`; the degree-0 nodes are exactly the pure-blue classes (`tc = bm`).

5. **Descent potentials.** `black_descent_20261020c/cycles.md` proves that the cyclic-triangle count descends along a black line from every strong black-active node; THM-4741 (D) shows that `H` has no such property. Conjecture: every non-pure-blue class with `tc >= 3` has a black neighbour with strictly smaller `tc`, except the strong classes with `tc = min S(n)`; untested (the census column for `tc`-descent is the next cheap addition).

Non-consequences: (1) does not bound `tc` for non-strong classes (products); (3) is a sample, not a certificate; (4) (ii) is conjectural beyond `n = 9`.
