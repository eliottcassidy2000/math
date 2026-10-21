---
id: HYP-9352
title: "After THM-4741: tc >= 5 for every strong tournament other than K1, C3, T5 (so tc = 3 classes are exactly the one-T5 stacks for all n); fibre sizes 7 and 21 never occur; the three-move law at n = 11 (random sample); growth of the strong spectrum minimum; the pendant-node count law a(n) - p(n) at n = 10"
status: "NUMERICAL ((1), (2), (4) exhaustive at n <= 9; (3) sampled at n = 11); OPEN"
session: opus-2026-10-20-S48
source: 05-knowledge/results/presentation_counts_single_view_classes_rotation_inequality_20261020.md
related:
  - 01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md
  - 05-knowledge/hypotheses/HYP-9347-three-move-law-layer-theorems-fundamental-black-graph-and-the-klein-dictionary-beyond-order-4.md
  - 05-knowledge/hypotheses/HYP-9346-self-lines-connectivity-and-growth-laws-of-the-tiling-line-metagraph.md
---

# HYP-9352: presentation-count laws after THM-4741

1. **Strong spectrum floor.** Every strong tournament other than `K1`, `C3`, `T5` has presentation count `tc >= 5`; equivalently `H(T) > 3 |Aut(T)|` for strong `T` of order `>= 6`. With THM-4741 (A2) this makes the `tc = 3` classes exactly the stacks over `{K1, C3, T5}` with one `T5` for every `n` (exhaustive at `n <= 9`: `1, 2, 3, 6, 11`). Sharper: the minimum of the strong spectrum `min S(n)` (`1, 5, 3, 9, 9, 17, 25` at `n = 3..9`) is non-decreasing for `n >= 5` and at least `2n - 1` for `n >= 8`. Hostile: a strong tournament with a large automorphism group and few Hamiltonian paths (the wreath products `C3[C3[...]]` have `|Aut| = 3^((n-1)/2)` but far more paths).

2. **Forbidden fibre sizes.** No node of the line metagraph has `tc = 7` or `tc = 21` at any order (THM-4741 (B) proves it for rigid nodes; a non-rigid counterexample needs `H = 7 |Aut|` or `21 |Aut|` with `|Aut| >= 3`, e.g. `H = 35, 49, 63`). Exhaustive at `n <= 9`. Companion: the tc-spectrum of all classes at order `n` is the set of odd numbers up to its maximum except 7 and 21, for every `n >= 9`.

3. **Three-move law at `n = 11`** (HYP-9347 (1)): every closed blue line between non-pure-blue classes has a literal three-black-move detour at `n = 11`. Evidence: 50,190 random closed blue lines (all layers) with three-move detours and no exception (THM-4741's note, section 6); the full census (`2^24` lines, about 160,000 core-seconds) was running at close-out.

4. **Pendant-node law.** The pendant nodes of the fundamental black graph number `a(n) - p(n)` (`a = A000930`, `p` the palindromic compositions into parts 1 and 3): `2, 4, 6, 10, 14` at `n = 5..9` (THM-4741 (A5), proved), and the degree-2 nodes of `B` are exactly the self-converse classes with `(tc, bm) = (5, 1)` (one at `n = 9`, `(H, |Aut|) = (15, 3)`); conjecture: the number of degree-`d` nodes of `B` for fixed `d` is eventually a quasi-polynomial in `n`.

5. **H is not a gradient, but `tc` might be.** Conjecture: every non-pure-blue class with `tc >= 3` has a black neighbour with strictly smaller `tc`, except the strong classes with `tc = min S(n)`; untested (the census column for `tc`-descent is the next cheap addition).

Non-consequences: (1) does not bound `tc` for non-strong classes (products); (3) is a sample, not a certificate.
